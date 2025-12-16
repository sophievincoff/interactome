import os
import pandas as pd
import logging
import os
import multiprocessing
from functools import partial
from math import ceil
from datetime import datetime
from hydra.core.hydra_config import HydraConfig
from omegaconf import DictConfig
from pathlib import Path
from lxml import etree as LXML_ET
from typing import List
import xml.etree.ElementTree as ET

import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)
ROOT_PATH = Path(root).resolve()  # reuse your existing 'root'

logger = logging.getLogger(__name__)

# --- XML utilities (namespace-safe) -----------------------------------------
def _S(tag: str) -> str:
    """Strip XML namespace: "{ns}tag" -> "tag"."""
    return tag.split("}", 1)[1] if tag.startswith("{") else tag

def _xml_text(elem, default=None):
    return (elem.text if elem is not None and elem.text is not None else default)

def _xml_attr(elem, name, default=None):
    return (elem.get(name) if elem is not None and elem.get(name) is not None else default)

def _iter_interaction_nodes_xml(root):
    # MIF300: interactions under entrySet/entry/interactionList/interaction
    for inter in root.iter():
        if _S(inter.tag) == "interaction":
            yield inter

def _safe_label(obj):
    """Return obj.label or obj.name, safely."""
    try:
        return getattr(obj, "label", None) or getattr(obj, "name", None)
    except Exception:
        return None
    
def _has_type_label(interactor, expected: str | List[str] ="protein"):
    """Returns True iff interactor.type has label==expected or label in expected; never raises."""
    if type(expected)==str:
        expected = [expected]
    try:
        t = interactor.type  # may explode internally
        return (_safe_label(t) in expected)
    except Exception:
        return False
    
def _safe_get(obj, attr, default=None):
        try:
            return getattr(obj, attr, default)
        except Exception:
            return default

def _build_interactor_map_xml(root):
    """
    Build dict: interactor @id -> rich info including:
      - db/acc from primaryRef
      - lists from secondaryRef: ensp/ensg/enst/interpro/rcsbpdb/go/intact_ids
      - names.shortLabel and gene_symbol from alias(typeAc=MI:0301 or type="gene name")
      - type label, sequence, and organism (optional)
    """
    inter_map = {}

    def _xml_text(elem, default=None):
        return (elem.text if elem is not None and elem.text is not None else default)

    def _xml_attr(elem, name, default=None):
        return (elem.get(name) if elem is not None and elem.get(name) is not None else default)

    def _push_lists(dbv, acv, out_lists):
        if not acv or not dbv:
            return
        dbn = dbv.strip().lower()

        if dbn == "uniprotkb":
            out_lists["uniprot_ids"].append(acv)

        # Ensembl
        elif dbn == "ensembl":
            if   acv.startswith("ENSP"): out_lists["ensp"].append(acv)
            elif acv.startswith("ENSG"): out_lists["ensg"].append(acv)
            elif acv.startswith("ENST"): out_lists["enst"].append(acv)

        elif dbn == "interpro":
            out_lists["interpro"].append(acv)

        elif dbn in ("rcsb pdb", "rscb pdb", "pdb", "protein databank"):
            out_lists["rcsbpdb"].append(acv)

        elif dbn in ("go", "gene ontology"):
            out_lists["go"].append(acv)

        elif dbn == "intact":
            out_lists["intact_ids"].append(acv)
            
        elif dbn == "reactome":
            out_lists["reactome"].append(acv)
            
        elif dbn == "dip":
            out_lists["dip"].append(acv)

    for inter in root.iter():
        if _S(inter.tag) != "interactor":
            continue

        iid = _xml_attr(inter, "id")
        if not iid:
            continue

        # names / shortLabel
        names = next((c for c in inter if _S(c.tag) == "names"), None)
        sl = next((c for c in names if _S(c.tag) == "shortLabel"), None) if names is not None else None
        short = _xml_text(sl)

        # alias → gene_symbol
        gene_symbol = None
        if names is not None:
            for al in names.findall(".//{*}alias"):
                t_ac  = _xml_attr(al, "typeAc")
                t_txt = (_xml_attr(al, "type") or "").strip().lower()
                if (t_ac == "MI:0301") or (t_txt == "gene name"):
                    gene_symbol = _xml_text(al)
                    if gene_symbol:
                        break

        # xref primary + secondary
        xref = next((c for c in inter if _S(c.tag) == "xref"), None)
        pref = None
        if xref is not None:
            for c in xref:
                if _S(c.tag) == "primaryRef":
                    pref = c
                    break
        db  = _xml_attr(pref, "db")
        acc = _xml_attr(pref, "id") or _xml_attr(pref, "ac")

        lists = {
            "ensp": [], "ensg": [], "enst": [],
            "interpro": [], "rcsbpdb": [], "reactome": [], "go": [], "intact_ids": [],
            "uniprot_ids": [], "dip": []
        }
        # include primary in lists (some files stash useful info there)
        _push_lists(db, acc, lists)

        if xref is not None:
            for c in xref:
                if _S(c.tag) != "secondaryRef":
                    continue
                _push_lists(_xml_attr(c, "db"), (_xml_attr(c, "id") or _xml_attr(c, "ac")), lists)

        # type label
        itype = next((c for c in inter if _S(c.tag) == "interactorType"), None)
        it_names = next((c for c in itype if _S(c.tag) == "names"), None) if itype is not None else None
        it_sl = next((c for c in it_names if _S(c.tag) == "shortLabel"), None) if it_names is not None else None
        type_label = _xml_text(it_sl)

        # organism (species), not experimental host
        org_el = next((c for c in inter if _S(c.tag) == "organism"), None)
        species_taxid = _xml_attr(org_el, "ncbiTaxId") if org_el is not None else None
        species_label = None
        if org_el is not None:
            on = next((c for c in org_el if _S(c.tag) == "names"), None)
            osl = next((c for c in on if _S(c.tag) == "shortLabel"), None) if on is not None else None
            species_label = _xml_text(osl)

        # sequence (optional)
        seq_el = next((c for c in inter if _S(c.tag) == "sequence"), None)
        seq = _xml_text(seq_el)
        
        # chain-seq-* attributes on the interactor (like IntAct bulk chains)
        chain_seq_start = None
        chain_seq_end = None
        # Look for any <attribute> named "chain-seq-start"/"chain-seq-end"
        # anywhere under this <interactor> (not just direct children).
        for attr in inter.iter():
            if _S(attr.tag) != "attribute":
                continue
            name = _xml_attr(attr, "name")
            if name == "chain-seq-start" and _xml_text(attr):
                chain_seq_start = _xml_text(attr).strip()
            elif name == "chain-seq-end" and _xml_text(attr):
                chain_seq_end = _xml_text(attr).strip()
                
        inter_map[iid] = {
            "db": db, "acc": acc, "short": short, "type": type_label, "sequence": seq,
            "gene_symbol": gene_symbol, "reactome": lists["reactome"],
            "ensp": lists["ensp"], "ensg": lists["ensg"], "enst": lists["enst"],
            "interpro": lists["interpro"], "rscbpdb": lists["rcsbpdb"], "go": lists["go"],
            "intact_ids": lists["intact_ids"],
            "uniprot_ids": lists["uniprot_ids"],
            "dip": lists["dip"],
            "species_taxid": species_taxid, "species_label": species_label,
            "chain_seq_start": chain_seq_start,
            "chain_seq_end": chain_seq_end,
        }

    return inter_map


def _build_experiment_map_xml(root):
    """
    Map experimentDescription @id -> node, so we can resolve <experimentRef>.
    """
    exp_map = {}
    for el in root.iter():
        if _S(el.tag) == "experimentDescription":
            eid = _xml_attr(el, "id")
            if eid:
                exp_map[eid] = el
    return exp_map

def _parse_experiment_node_xml(exp_node):
    """
    Mirror parse_psi30_experiment() return fields from raw XML experimentDescription.

    Returns dict like:
      {
        "pubmed": ...,
        "method": ...,
        "partmethod": ...,
        "hosts": [
          {
            "taxid": "9606",
            "label": "human",           # shortLabel (for convenience)
            "label_short": "human",     # explicit shortLabel
            "label_full": "Homo sapiens",  # fullName, if present
            "cell_type": "HeLa",
            "compartment": "cytoplasm",
            "tissue": "liver",
          },
          ...
        ]
      }
    """
    new_experiment = {}
    # pubmed
    bibref = next((c for c in exp_node if _S(c.tag) == 'bibref'), None)
    xref = next((c for c in bibref if _S(c.tag) == 'xref'), None) if bibref is not None else None
    pref = None
    if xref is not None:
        for c in xref:
            if _S(c.tag) == 'primaryRef':
                pref = c
                break
    if pref is not None and _xml_attr(pref, 'db') == 'pubmed':
        new_experiment["pubmed"] = _xml_attr(pref, 'id') or _xml_attr(pref, 'ac')

    # interactionDetectionMethod
    method = next((c for c in exp_node if _S(c.tag) == 'interactionDetectionMethod'), None)
    if method is not None:
        nms = next((c for c in method if _S(c.tag) == 'names'), None)
        sl  = next((c for c in nms if _S(c.tag) == 'shortLabel'), None) if nms is not None else None
        if _xml_text(sl):
            new_experiment["method"] = _xml_text(sl)

    # participantIdentificationMethod (partmethod)
    pim = next((c for c in exp_node if _S(c.tag) == 'participantIdentificationMethod'), None)
    if pim is not None:
        nms = next((c for c in pim if _S(c.tag) == 'names'), None)
        sl  = next((c for c in nms if _S(c.tag) == 'shortLabel'), None) if nms is not None else None
        if _xml_text(sl):
            new_experiment["partmethod"] = _xml_text(sl)

    # hosts
    host_list = next((c for c in exp_node if _S(c.tag) == 'hostOrganismList'), None)
    if host_list is not None:
        hosts = []
        for host in host_list:
            if _S(host.tag) != 'hostOrganism':
                continue

            h = {}
            taxid = _xml_attr(host, 'ncbiTaxId') or _xml_attr(host, 'taxid')
            if taxid is not None:
                h["taxid"] = taxid

            # names → shortLabel + fullName
            nms = next((c for c in host if _S(c.tag) == 'names'), None)
            if nms is not None:
                sl = next((c for c in nms if _S(c.tag) == 'shortLabel'), None)
                fl = next((c for c in nms if _S(c.tag) == 'fullName'), None)

                short_label = _xml_text(sl)
                full_label  = _xml_text(fl)

                if short_label:
                    h["label_short"] = short_label
                    # keep old 'label' for convenience / backward compat
                    h["label"] = short_label
                elif full_label:
                    # fallback if only fullName exists
                    h["label"] = full_label

                if full_label:
                    h["label_full"] = full_label

            # cellType, compartment, tissue (optional child controlled vocab terms)
            for subtag, key in (("cellType", "cell_type"),
                                ("compartment", "compartment"),
                                ("tissue", "tissue")):
                sub = next((c for c in host if _S(c.tag) == subtag), None)
                if sub is not None:
                    sn = next((c for c in sub if _S(c.tag) == 'names'), None)
                    ssl = next((c for c in sn if _S(c.tag) == 'shortLabel'), None) if sn is not None else None
                    lbl = _xml_text(ssl)
                    if lbl:
                        h[key] = lbl

            hosts.append(h)

        if hosts:
            new_experiment["hosts"] = hosts

    return new_experiment

def _xml_to_interactor_info(x):
    # Collect all uniprot IDs: primaryRef + all collected 'uniprot_ids'
    uniprot_list = []
    # put primary acc first
    if x.get("db") == "uniprotkb" and x.get("acc"):
        uniprot_list.append(x["acc"])
    uniprot_list.extend(x.get("uniprot_ids", []))

    # de-duplicate while preserving order
    seen = set()
    uniprot_dedup = []
    for u in uniprot_list:
        if not u:
            continue
        if u not in seen:
            seen.add(u)
            uniprot_dedup.append(u)

    uniprotkb = ",".join(uniprot_dedup) or None
    out = {
        "mol_type": x.get("type"),
        "gene_name": x.get("short"),
        "gene_symbol": x.get("gene_symbol"),
        "length": len(x.get("sequence")) if isinstance(x.get("sequence"), str) else None,
        "protein": x.get("sequence"),
        "uniprotkb": uniprotkb,
        "ensp": ",".join(x.get("ensp", [])) or None,
        "ensg": ",".join(x.get("ensg", [])) or None,
        "enst": ",".join(x.get("enst", [])) or None,
        "interpro": ",".join(x.get("interpro", [])) or None,
        "reactome": ",".join(x.get("reactome", [])) or None,
        "rscbpdb": ",".join(x.get("rscbpdb", [])) or None,
        "intactid": ",".join(x.get("intact_ids", [])) or None,
        "dip": ",".join(x.get("dip", [])) or None,
        "primaryref_db": x.get("db"),
        "primaryref_id": x.get("acc"),
        "go": ",".join(x.get("go", [])) or None,
        # (Optional) you can add species fields if you want them:
        "species_taxid": x.get("species_taxid"),
        "species_label": x.get("species_label"),
        "chain_seq_start": x.get("chain_seq_start"),
        "chain_seq_end": x.get("chain_seq_end"),
    }

    # Optional: feature fields that may have been attached at participant level
    for k in [
        "mutation_mi", "mutation_name", "mutation_begin", "mutation_end",
        "mutation_orig", "mutation_new", "mutation_short",
        "binding_mi", "binding_name", "binding_begin", "binding_end", "binding_short",
        "ptm_mi", "ptm_name", "ptm_begin", "ptm_end", "ptm_orig", "ptm_new", "ptm_short",
    ]:
        out[k] = x.get(k)

    return out
    
def _xml_interaction_label(inter_el):
    # Prefer fullName, fallback to shortLabel
    full  = inter_el.findtext(".//{*}interactionType/{*}names/{*}fullName")
    short = inter_el.findtext(".//{*}interactionType/{*}names/{*}shortLabel")
    return (full or short or "").strip()

def _xml_interaction_mi(inter_el):
    pref = inter_el.find(".//{*}interactionType/{*}xref/{*}primaryRef")
    if pref is None: 
        return None
    # Some MIFs use `id`, others `ac`; prefer `id`, fallback to `ac`
    return pref.get("id") or pref.get("ac")
    
def _extract_features_from_xml(inter_node, mutation_mi_ok, bindsite_mi_ok, ptm_mi_ok):
    """
    Backwards-compatible wrapper that returns *_1 / *_2 keys
    but delegates the real work to _extract_features_from_participant.
    """
    out = {
        "mutation_mi_1": None, "mutation_name_1": None,
        "mutation_begin_1": None, "mutation_end_1": None,
        "mutation_orig_1": None, "mutation_new_1": None,
        "binding_mi_1": None, "binding_name_1": None,
        "binding_begin_1": None, "binding_end_1": None,
        "chain_seq_start_1": None, "chain_seq_end_1": None,
        "ptm_mi_1": None, "ptm_name_1": None,
        "ptm_begin_1": None, "ptm_end_1": None,
        "ptm_orig_1": None, "ptm_new_1": None,

        "mutation_mi_2": None, "mutation_name_2": None,
        "mutation_begin_2": None, "mutation_end_2": None,
        "mutation_orig_2": None, "mutation_new_2": None,
        "binding_mi_2": None, "binding_name_2": None,
        "binding_begin_2": None, "binding_end_2": None,
        "chain_seq_start_2": None, "chain_seq_end_2": None,
        "ptm_mi_2": None, "ptm_name_2": None,
        "ptm_begin_2": None, "ptm_end_2": None,
        "ptm_orig_2": None, "ptm_new_2": None,
    }

    plist = next((c for c in inter_node if _S(c.tag) == "participantList"), None)
    parent = plist if plist is not None else inter_node
    parts = [p for p in parent if _S(p.tag) == "participant"]
    if not parts:
        return out

    if len(parts) >= 1:
        p1 = _extract_features_from_participant(parts[0], mutation_mi_ok, bindsite_mi_ok, ptm_mi_ok)
        for k, v in p1.items():
            out[f"{k}_1"] = v

    if len(parts) >= 2:
        p2 = _extract_features_from_participant(parts[1], mutation_mi_ok, bindsite_mi_ok, ptm_mi_ok)
        for k, v in p2.items():
            out[f"{k}_2"] = v

    return out

def _xml_interaction_intact_id(inter_el):
    """
    Given an <interaction> element, return its IntAct interaction ID
    (xref/primaryRef db="intact") or None.
    """
    xref = next((c for c in inter_el if _S(c.tag) == "xref"), None)
    if xref is None:
        return None
    pref = next((c for c in xref if _S(c.tag) == "primaryRef"), None)
    if pref is None:
        return None
    if _xml_attr(pref, "db") == "intact":
        return _xml_attr(pref, "id") or _xml_attr(pref, "ac")
    return None

def _xml_interaction_all_intact_ids(inter_el):
    """
    Given an <interaction> element, return *all* IntAct interaction IDs
    in its <xref> (primaryRef + secondaryRef with db="intact"), as a list.
    """
    ids = []
    xref = next((c for c in inter_el if _S(c.tag) == "xref"), None)
    if xref is None:
        return ids

    # primaryRef
    pref = next((c for c in xref if _S(c.tag) == "primaryRef"), None)
    if pref is not None and _xml_attr(pref, "db") == "intact":
        val = _xml_attr(pref, "id") or _xml_attr(pref, "ac")
        if val:
            ids.append(val)

    # all secondaryRef
    for c in xref:
        if _S(c.tag) != "secondaryRef":
            continue
        if _xml_attr(c, "db") == "intact":
            val = _xml_attr(c, "id") or _xml_attr(c, "ac")
            if val:
                ids.append(val)

    # de-dup while preserving order
    seen = set()
    out = []
    for v in ids:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out

def _build_interaction_xml_index(root):
    """
    Build:
      - list of all <interaction> nodes (for legacy index-based usage)
      - dict: IntAct interaction id -> <interaction> node
    """
    nodes = []
    by_intact = {}
    for inter in _iter_interaction_nodes_xml(root):
        nodes.append(inter)
        intact_id = _xml_interaction_intact_id(inter)
        if intact_id:
            by_intact[intact_id] = inter
    return nodes, by_intact
def _extract_features_from_participant(
    part,
    mutation_mi_ok,
    bindsite_mi_ok,
    ptm_mi_ok,
):
    """
    Extract mutation / binding-site / PTM / chain-seq-* features
    for a single <participant> node.

    Returns a dict with UNSUFFIXED keys:

      # mutation features: one pipe per <feature>
      mutation_mi, mutation_name, mutation_short,
      mutation_begin, mutation_end,
      mutation_orig, mutation_new,

      # binding-site features: one pipe per <feature>
      binding_mi, binding_name, binding_short,
      binding_begin, binding_end,

      # PTM features: one pipe per <feature>
      ptm_mi, ptm_name, ptm_short,
      ptm_begin, ptm_end,
      ptm_orig, ptm_new,

      # plus chain-seq-* on the participant
      chain_seq_start, chain_seq_end
    """

    # Collect full feature records so pipes stay aligned
    mutation_features = []  # each: {"mi","name","short","begin","end","orig","new"}
    binding_features  = []  # each: {"mi","name","short","begin","end"}
    ptm_features      = []  # each: {"mi","name","short","begin","end","orig","new"}

    ptm_mi_ok = set(ptm_mi_ok or [])

    # --- walk all <feature> under this participant ---
    for feat in part.iter():
        if _S(feat.tag) != "feature":
            continue

        # --- feature-level shortLabel (e.g., O60566:p.[Leu128Ala;Leu131Ala]) ---
        feat_names = feat.find(".//{*}names")
        feat_sl = feat_names.find(".//{*}shortLabel") if feat_names is not None else None
        feat_short = _xml_text(feat_sl)

        # --- featureType MI + type name (e.g., mutation disrupting interaction) ---
        pref = feat.find(".//{*}featureType/{*}xref/{*}primaryRef")
        ftype_mi = None
        psi_mod_id = None
        if pref is not None:
            dbv = pref.get("db")
            if dbv == "psi-mi":
                ftype_mi = pref.get("id") or pref.get("ac")
            elif dbv == "psi-mod":
                # e.g. MOD:00394 monoacetylated residue
                psi_mod_id = pref.get("id") or pref.get("ac")

        ftype_name = (
            feat.findtext(".//{*}featureType/{*}names/{*}fullName")
            or feat.findtext(".//{*}featureType/{*}names/{*}shortLabel")
        )

        # --- featureRole MI + name (e.g., observed-ptm) ---
        role_pref = feat.find(".//{*}featureRole/{*}xref/{*}primaryRef")
        role_mi = None
        if role_pref is not None and (role_pref.get("db") == "psi-mi"):
            role_mi = role_pref.get("id") or role_pref.get("ac")

        role_name = (
            feat.findtext(".//{*}featureRole/{*}names/{*}fullName")
            or feat.findtext(".//{*}featureRole/{*}names/{*}shortLabel")
        )

        is_mut  = ftype_mi in mutation_mi_ok
        is_bind = ftype_mi in bindsite_mi_ok

        # --- PTM flag via attributeList, featureRole, featureType, or psi-mod ---
        is_ptm = False
        ptm_attrs_mi = []
        ptm_attrs_name = []

        # (A) Original attributeList-based PTM detection
        attr_list = feat.find(".//{*}attributeList")
        if attr_list is not None:
            for attr in attr_list:
                if _S(attr.tag) != "attribute":
                    continue
                nameAc = attr.get("nameAc")
                nameTxt = attr.get("name")
                if nameAc and nameAc in ptm_mi_ok:
                    is_ptm = True
                    ptm_attrs_mi.append(nameAc)
                    if nameTxt:
                        ptm_attrs_name.append(nameTxt)

        # (B) PTM via featureRole MI term (e.g., MI:0925 observed-ptm subtree)
        if (not is_ptm) and role_mi and role_mi in ptm_mi_ok:
            is_ptm = True
            ptm_attrs_mi.append(role_mi)
            if role_name:
                ptm_attrs_name.append(role_name)

        # (C) PTM via featureType MI term in ptm_mi_ok (if your CV includes these)
        if (not is_ptm) and ftype_mi and ftype_mi in ptm_mi_ok:
            is_ptm = True
            ptm_attrs_mi.append(ftype_mi)
            if ftype_name:
                ptm_attrs_name.append(ftype_name)

        # (D) PTM via psi-mod db (e.g., MOD:00394 monoacetylated residue)
        if (not is_ptm) and psi_mod_id:
            is_ptm = True
            ptm_attrs_mi.append(psi_mod_id)
            if ftype_name:
                ptm_attrs_name.append(ftype_name)

        # --- per-feature ranges + residues ---
        begins_this_feat, ends_this_feat = [], []
        mut_origs_this_feat, mut_news_this_feat = [], []
        ptm_origs_this_feat, ptm_news_this_feat = [], []

        frl = feat.find(".//{*}featureRangeList")
        if frl is not None:
            for fr in frl:
                if _S(fr.tag) != "featureRange":
                    continue
                begin_el = fr.find(".//{*}begin")
                end_el   = fr.find(".//{*}end")
                bpos = begin_el.get("position") if (begin_el is not None and begin_el.get("position")) else None
                epos = end_el.get("position") if (end_el is not None and end_el.get("position")) else None
                begins_this_feat.append(bpos)
                ends_this_feat.append(epos)

                res = fr.find(".//{*}resultingSequence")
                if res is not None:
                    orig = res.findtext(".//{*}originalSequence")
                    new  = res.findtext(".//{*}newSequence")
                else:
                    orig = None
                    new  = None

                if is_mut:
                    if orig is not None:
                        mut_origs_this_feat.append(orig)
                    if new is not None:
                        mut_news_this_feat.append(new)

                if is_ptm:
                    if orig is not None:
                        ptm_origs_this_feat.append(orig)
                    if new is not None:
                        ptm_news_this_feat.append(new)

        def _join_or_none(xs):
            xs = [x for x in xs if x is not None]
            return ",".join(xs) if xs else None

        begin_str = _join_or_none(begins_this_feat)
        end_str   = _join_or_none(ends_this_feat)

        # --- stash per-feature records ---
        if is_mut:
            mutation_features.append({
                "mi":    ftype_mi,
                "name":  ftype_name,
                "short": feat_short,
                "begin": begin_str,
                "end":   end_str,
                "orig":  _join_or_none(mut_origs_this_feat),
                "new":   _join_or_none(mut_news_this_feat),
            })

        if is_bind:
            binding_features.append({
                "mi":    ftype_mi,
                "name":  ftype_name,
                "short": feat_short,
                "begin": begin_str,
                "end":   end_str,
            })

        if is_ptm:
            # --- build PTM names: attributes + role + featureType fullName ---
            def _dedup(seq):
                out = []
                seen = set()
                for v in seq:
                    if not v:
                        continue
                    if v not in seen:
                        seen.add(v)
                        out.append(v)
                return out

            names_for_this_ptm = list(ptm_attrs_name)
            if role_name:
                names_for_this_ptm.append(role_name)
            if ftype_name:
                names_for_this_ptm.append(ftype_name)

            ptm_mi_val = _join_or_none(_dedup(ptm_attrs_mi))
            ptm_name_val = _join_or_none(_dedup(names_for_this_ptm))

            ptm_features.append({
                "mi":    ptm_mi_val,
                "name":  ptm_name_val,
                "short": feat_short,
                "begin": begin_str,
                "end":   end_str,
                "orig":  _join_or_none(ptm_origs_this_feat),
                "new":   _join_or_none(ptm_news_this_feat),
            })

    # --- chain-seq-* attributes on participant ---
    cstart = None
    cend = None
    for attr in part.iter():
        if _S(attr.tag) != "attribute":
            continue
        name = attr.get("name")
        if name == "chain-seq-start" and attr.text:
            cstart = attr.text.strip()
        elif name == "chain-seq-end" and attr.text:
            cend = attr.text.strip()

    # --- helper: join field across features WITHOUT losing alignment ---
    def _join_feature_field(feats, key):
        if not feats:
            return None
        # keep empty slots so number of pipes stays equal to len(feats) - 1
        vals = []
        for f in feats:
            v = f.get(key)
            vals.append("" if v is None else str(v))
        return "|".join(vals)

    return {
        # mutation features
        "mutation_mi":    _join_feature_field(mutation_features, "mi"),
        "mutation_name":  _join_feature_field(mutation_features, "name"),
        "mutation_short": _join_feature_field(mutation_features, "short"),
        "mutation_begin": _join_feature_field(mutation_features, "begin"),
        "mutation_end":   _join_feature_field(mutation_features, "end"),
        "mutation_orig":  _join_feature_field(mutation_features, "orig"),
        "mutation_new":   _join_feature_field(mutation_features, "new"),

        # binding-site features
        "binding_mi":     _join_feature_field(binding_features, "mi"),
        "binding_name":   _join_feature_field(binding_features, "name"),
        "binding_short":  _join_feature_field(binding_features, "short"),
        "binding_begin":  _join_feature_field(binding_features, "begin"),
        "binding_end":    _join_feature_field(binding_features, "end"),

        # PTM features
        "ptm_mi":     _join_feature_field(ptm_features, "mi"),
        "ptm_name":   _join_feature_field(ptm_features, "name"),
        "ptm_short":  _join_feature_field(ptm_features, "short"),
        "ptm_begin":  _join_feature_field(ptm_features, "begin"),
        "ptm_end":    _join_feature_field(ptm_features, "end"),
        "ptm_orig":   _join_feature_field(ptm_features, "orig"),
        "ptm_new":    _join_feature_field(ptm_features, "new"),

        # participant-level chain info
        "chain_seq_start": cstart,
        "chain_seq_end":   cend,
    }

def organize_cv(
    path="/scratch/pranamlab/sophie/interactome/interactome/data_files/raw/intact/cv/intact.obo",
    starting_id="MI:0190",
    out_dir="/scratch/pranamlab/sophie/interactome/interactome/data_files/processed/intact/cv",
):
    """
    Organize IntAct PSI-MI OBO into a CSV subtree rooted at `starting_id`.

    Output columns:
      - id:              term id
      - label:           term name
      - parent_id:       the direct parent along the traversed edge (one per row)
      - parent_ids_all:  comma-separated *all* is_a parent ids for this term
      - parent_names_all:comma-separated *all* is_a parent names for this term
    """

    # -------- Parse OBO into terms dict --------
    # terms[id] = {"name": str, "is_a": [parent_ids]}
    terms = {}
    current = None

    def _commit(curr):
        if not curr:
            return
        if "id" in curr and "name" in curr:
            terms[curr["id"]] = {
                "name": curr["name"],
                "is_a": curr.get("is_a", []),
            }

    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if line == "[Term]":
                # Commit previous term before starting a new one
                _commit(current)
                current = {}
                continue
            if not line:
                # Blank line ends a stanza
                _commit(current)
                current = None
                continue

            if current is None:
                # Skip header or [Typedef] sections, etc.
                continue

            if line.startswith("id:"):
                current["id"] = line[3:].strip()
            elif line.startswith("name:"):
                current["name"] = line[5:].strip()
            elif line.startswith("is_a:"):
                # Keep only the ID portion before "!"
                parent_id = line[5:].split("!")[0].strip()
                current.setdefault("is_a", []).append(parent_id)

    # Commit last term if file doesn"t end with a blank line
    _commit(current)

    # -------- Build reverse index for traversal --------
    children_map = {}
    for term_id, meta in terms.items():
        for p in meta.get("is_a", []):
            children_map.setdefault(p, []).append(term_id)

    # -------- Traverse subtree from starting_id --------
    results = []

    def traverse(term_id, parent_id=None):
        if term_id not in terms:
            return
        meta = terms[term_id]
        all_parents = meta.get("is_a", [])
        # Map parent IDs → names if present
        parent_names = [terms[p]["name"] for p in all_parents if p in terms]

        results.append({
            "label": meta["name"],
            "id": term_id,
            "parent_id": parent_id,  # edge used to reach this node
            "parent_ids_all": ", ".join(all_parents) if all_parents else "",
            "parent_names_all": ", ".join(parent_names) if parent_names else "",
        })

        for child in children_map.get(term_id, []):
            traverse(child, parent_id=term_id)

    traverse(starting_id)

    df = pd.DataFrame(results).drop_duplicates()

    os.makedirs(out_dir, exist_ok=True)
    out_path = f"{out_dir}/{starting_id.lower().replace(':','_')}_subtree.csv"
    df.to_csv(out_path, index=False)
    return out_path
        
def run_debug_mode(
    debug_file: str,
    interaction_milabel_ok: list,
    mutation_mi_ok,
    bindsite_mi_ok,
    ptm_mi_ok,
    output_dir: Path,
    ok_part_types: str | List[str] = ("protein", "peptide"),
):
    logger.info(f"--- Running XML DEBUG mode on: {debug_file} ---")
    debug_file = str(debug_file)
    example_fname = Path(debug_file).stem

    # normalize ok_part_types to a set of strings
    if isinstance(ok_part_types, str):
        ok_part_types = [ok_part_types]
    ok_part_types = set(ok_part_types)

    try:
        tree = ET.parse(debug_file)
        root = tree.getroot()
    except Exception as e:
        logger.error(f"Cannot parse debug XML {debug_file}: {e}")
        return

    xml_inter_nodes, xml_inter_by_intact = _build_interaction_xml_index(root)
    xml_inter_map   = _build_interactor_map_xml(root)
    xml_exp_map     = _build_experiment_map_xml(root)

    interaction_mi_norm = {
        str(x).strip().lower()
        for x in interaction_milabel_ok
        if pd.notna(x)
    }

    logger.info(f"File has {len(xml_inter_nodes)} <interaction> nodes (XML).")

    # Log the first few interactions that pass filters
    max_to_log = 5
    logged = 0

    for i, node in enumerate(xml_inter_nodes):
        if logged >= max_to_log:
            break

        label_xml = _xml_interaction_label(node)
        if (label_xml or "").strip().lower() not in interaction_mi_norm:
            continue

        plist = next((c for c in node if _S(c.tag) == "participantList"), None)
        part_parent = plist if plist is not None else node

        xml_parts = []
        for p in part_parent:
            if _S(p.tag) != "participant":
                continue

            info = None

            # 1) Try referenced interactor
            ir = p.find(".//{*}interactorRef")
            if ir is not None and ir.text:
                info = xml_inter_map.get(ir.text)

            # 2) Inline interactor
            if info is None:
                inter_el = p.find(".//{*}interactor")
                if inter_el is not None:
                    iid = _xml_attr(inter_el, "id")
                    if iid and iid in xml_inter_map:
                        info = xml_inter_map[iid]
                    else:
                        # synthesize a dict for inline (db/acc/short/type/sequence)
                        names = inter_el.find(".//{*}names")
                        sl = names.find(".//{*}shortLabel") if names is not None else None
                        short = _xml_text(sl)

                        # gene_symbol from alias(typeAc=MI:0301) or type="gene name"
                        gene_symbol = None
                        if names is not None:
                            for al in names.findall(".//{*}alias"):
                                t_ac  = _xml_attr(al, "typeAc")
                                t_txt = _xml_attr(al, "type")
                                if (t_ac == "MI:0301") or (t_txt and t_txt.strip().lower() == "gene name"):
                                    gene_symbol = _xml_text(al)
                                    break

                        xref = inter_el.find(".//{*}xref")
                        pref = xref.find(".//{*}primaryRef") if xref is not None else None
                        db  = _xml_attr(pref, "db")
                        acc = _xml_attr(pref, "id") or _xml_attr(pref, "ac")

                        ensp = []; ensg = []; enst = []
                        interpro = []; rscbpdb = []; go = []; reactome = []; intact_ids = []
                        uniprot_ids = []; dip = []

                        def _push(dbv, acv):
                            if not acv or not dbv:
                                return
                            dbn = dbv.strip().lower()
                            if dbn == "uniprotkb":
                                uniprot_ids.append(acv)
                            elif dbn == "ensembl":
                                if   acv.startswith("ENSP"): ensp.append(acv)
                                elif acv.startswith("ENSG"): ensg.append(acv)
                                elif acv.startswith("ENST"): enst.append(acv)
                            elif dbn == "interpro":
                                interpro.append(acv)
                            elif dbn in ("rcsb pdb", "rscb pdb", "pdb", "protein databank"):
                                rscbpdb.append(acv)
                            elif dbn in ("go", "gene ontology"):
                                go.append(acv)
                            elif dbn == "intact":
                                intact_ids.append(acv)
                            elif dbn == "reactome":
                                    intact_ids.append(acv)
                            elif dbn == "dip":
                                    dip.append(acv)

                        _push(db, acc)
                        if xref is not None:
                            for c in xref:
                                if _S(c.tag) != "secondaryRef":
                                    continue
                                _push(
                                    _xml_attr(c, "db"),
                                    _xml_attr(c, "id") or _xml_attr(c, "ac"),
                                )

                        # organism (species), not experimental host
                        org_el = inter_el.find(".//{*}organism")
                        species_taxid = _xml_attr(org_el, "ncbiTaxId") if org_el is not None else None
                        species_label = None
                        if org_el is not None:
                            on = org_el.find(".//{*}names")
                            osl = on.find(".//{*}shortLabel") if on is not None else None
                            species_label = _xml_text(osl)

                        itype2 = inter_el.find(".//{*}interactorType")
                        it_names = itype2.find(".//{*}names") if itype2 is not None else None
                        it_sl = it_names.find(".//{*}shortLabel") if it_names is not None else None
                        type_label = _xml_text(it_sl)

                        seq_el = inter_el.find(".//{*}sequence")
                        seq = _xml_text(seq_el)

                        info = {
                            "db": db, "acc": acc, "short": short, "type": type_label,
                            "sequence": seq, 
                            "gene_symbol": gene_symbol,
                            "ensp": ensp, "ensg": ensg, "enst": enst, "reactome": reactome,
                            "interpro": interpro, "rscbpdb": rscbpdb, "go": go,
                            "intact_ids": intact_ids,
                            "uniprot_ids": uniprot_ids, "dip": dip,
                            "species_taxid": species_taxid, "species_label": species_label,
                            
                        }

            if info is None:
                continue

            # 3) Attach participant-specific feature info
            feat = _extract_features_from_participant(
                p, mutation_mi_ok, bindsite_mi_ok, ptm_mi_ok
            )
            info = dict(info)  # avoid mutating shared map entries
            for k in ("chain_seq_start", "chain_seq_end"):
                    if feat.get(k) is None:
                        feat.pop(k, None)
            info.update(feat)

            xml_parts.append(info)

        # Require exactly two participants of allowed types
        if len(xml_parts) != 2:
            continue
        if not (
            (xml_parts[0].get("type") in ok_part_types) and
            (xml_parts[1].get("type") in ok_part_types)
        ):
            continue

        logged += 1

        # Normalize via _xml_to_interactor_info so we see exactly
        # what parse_psi30 will emit into the final table
        info_1 = _xml_to_interactor_info(xml_parts[0])
        info_2 = _xml_to_interactor_info(xml_parts[1])

        logger.info(f"\nInteraction {i} [xml]")
        logger.info(f"  label: {label_xml}")
        logger.info(
            f"  participant types: {info_1.get('mol_type')}, {info_2.get('mol_type')}"
        )
        logger.info(
            f"  participant 1: gene_symbol={info_1.get('gene_symbol')} "
            f"uniprot={info_1.get('uniprotkb')}"
        )
        logger.info(
            f"  participant 2: gene_symbol={info_2.get('gene_symbol')} "
            f"uniprot={info_2.get('uniprotkb')}"
        )

        logger.info(
            "  features p1: "
            f"mut={info_1.get('mutation_mi')} "
            f"bind={info_1.get('binding_mi')} "
            f"ptm={info_1.get('ptm_mi')} "
            f"chain=({info_1.get('chain_seq_start')}, {info_1.get('chain_seq_end')})"
        )
        logger.info(
            "  features p2: "
            f"mut={info_2.get('mutation_mi')} "
            f"bind={info_2.get('binding_mi')} "
            f"ptm={info_2.get('ptm_mi')} "
            f"chain=({info_2.get('chain_seq_start')}, {info_2.get('chain_seq_end')})"
        )

    # Also run full parse_psi30 and dump CSV for inspection
    logger.info("\n--- Running parse_psi30(XML-only) on debug file ---")
    test_info, test_skip_info = parse_psi30(
        [debug_file],
        interaction_milabel_ok,
        mutation_mi_ok,
        bindsite_mi_ok,
        ptm_mi_ok,
        ok_part_types=list(ok_part_types),
    )
    df = pd.DataFrame.from_dict(test_info).rename(columns={
        "gene_name_1": "protein_1",
        "gene_name_2": "protein_2",
        "protein_1": "aa_1",
        "protein_2": "aa_2",
    })
    output_dir.mkdir(parents=True, exist_ok=True)
    debug_path = output_dir / f"DEBUG_{example_fname}_parsed.csv"
    df.to_csv(debug_path, index=False)
    logger.info(f"Saved debug output to {debug_path}")
    
    test_skip_info = test_skip_info["skipped_interactions"]
    if len(test_skip_info)>0:
        skipped_df = pd.DataFrame(test_skip_info)
        skipped_path = os.path.join(output_dir, f"DEBUG_{example_fname}_skipped_interactions.csv")
        skipped_df.to_csv(skipped_path, index=False)    
        logger.info(f"Saved debug output SKIPPED interactions to {debug_path}") 

def parse_psi30_interactor(interactor):
    """
    Parse a single interactor from PSI-MI 3.0 format, resilient to missing
    Pymex internals (e.g., names["alias"] -> KeyError).
    """
    if interactor is None:
        return {
            "gene_name": None, "gene_symbol": None, "length": None, "protein": None,
            "uniprotkb": None, "ensp": None, "ensg": None, "enst": None, "interpro": None,
            "rscbpdb": None, "intactid": None, "primaryref_db": None, "primaryref_id": None,
            "go": None,
        }

    # SAFELY access all Pymex properties (some can raise internally)
    mol_type    = _safe_label(interactor.type)
    gene_name   = _safe_get(interactor, "label")
    gene_symbol = _safe_get(interactor, "gene")          # <-- was raising KeyError: "alias"
    protein     = _safe_get(interactor, "sequence")
    length      = len(protein) if isinstance(protein, str) else None

    go = []; uniprotkb = []; ensp = []; ensg = []; enst = []; interpro = []; rscbpdb = []; intactid = []; dip = []
    primaryref_db = None; primaryref_id = None

    # primaryRef (safe)
    pref = _safe_get(interactor, "primaryRef")
    if pref is not None:
        try:
            database = _safe_get(pref, "db")
            primaryref_db = database
            primaryref_id = _safe_get(pref, "ac")
            if database == "uniprotkb":
                if primaryref_id: uniprotkb.append(primaryref_id)
            elif database == "ensembl" and primaryref_id:
                if   primaryref_id.startswith("ENSP"): ensp.append(primaryref_id)
                elif primaryref_id.startswith("ENSG"): ensg.append(primaryref_id)
                elif primaryref_id.startswith("ENST"): enst.append(primaryref_id)
            elif database == "interpro":    interpro.append(primaryref_id)
            elif database == "rscb pdb":    rscbpdb.append(primaryref_id)
            elif database == "go":          go.append(primaryref_id)
            elif database == "intact":      intactid.append(primaryref_id)
            elif database == "dip":         dip.append(primaryref_id)
        except Exception:
            pass

    # secondaryRef (safe iterate)
    secs = _safe_get(interactor, "secondaryRef") or []
    try:
        for sec in secs:
            database = _safe_get(sec, "db")
            ac = _safe_get(sec, "ac")
            if database == "uniprotkb": uniprotkb.append(ac)
            elif database == "ensembl" and ac:
                if   ac.startswith("ENSP"): ensp.append(ac)
                elif ac.startswith("ENSG"): ensg.append(ac)
                elif ac.startswith("ENST"): enst.append(ac)
            elif database == "interpro": interpro.append(ac)
            elif database == "rscb pdb": rscbpdb.append(ac)
            elif database == "go":       go.append(ac)
            elif database == "intact":   intactid.append(ac)
            elif database == "dip":      dip.append(ac)
    except Exception:
        pass

    return {
        "gene_name": gene_name, "gene_symbol": gene_symbol, "length": length, "protein": protein, "mol_type": mol_type,
        "uniprotkb": ",".join([x for x in uniprotkb if x]) or None,
        "ensp":      ",".join([x for x in ensp if x]) or None,
        "ensg":      ",".join([x for x in ensg if x]) or None,
        "enst":      ",".join([x for x in enst if x]) or None,
        "interpro":  ",".join([x for x in interpro if x]) or None,
        "rscbpdb":   ",".join([x for x in rscbpdb if x]) or None,
        "intactid":  ",".join([x for x in intactid if x]) or None,
        "dip":       ",".join([x for x in dip if x]) or None,
        "primaryref_db": primaryref_db, "primaryref_id": primaryref_id,
        "go": ",".join([x for x in go if x]) or None,
    }


def parse_psi30_experiment(experiment):
    """
    Parse a single PSI-30 experiment
    """
    new_experiment = {}
    bib = getattr(experiment, "bibref", None)
    pref = getattr(bib, "primaryRef", None) if bib else None
    try:
        if pref is not None:
            # PubMed ID
            new_experiment["pubmed"] = pref.ac if getattr(pref, "db", None) == "pubmed" else None

        # Experiment info
        method = getattr(experiment, "method", None)
        if method and getattr(method, "name", None) is not None:
            new_experiment["method"] = method.name

        partmethod = getattr(experiment, "partmethod", None)
        if partmethod and getattr(partmethod, "name", None) is not None:
            new_experiment["partmethod"] = partmethod.name

        hosts = getattr(experiment, "hosts", []) or []
        if hosts:
            new_experiment["hosts"] = []
            for host in hosts:
                new_host = {}
                taxid = getattr(host, "taxid", None)
                if taxid is not None:
                    new_host["taxid"] = taxid
                    if taxid != "-1":
                        # label, cell type, compartment, tissue
                        if getattr(host, "label", None) is not None:
                            new_host["label"] = host.label
                        if getattr(host, "cellType", None) is not None:
                            new_host["cell_type"] = host.cellType
                        if getattr(host, "compartment", None) is not None:
                            new_host["compartment"] = host.compartment
                        tissue = getattr(host, "tissue", None)
                        if tissue and getattr(tissue, "label", None) is not None:
                            new_host["tissue"] = tissue.label
                        new_experiment["hosts"].append(new_host)
                    else:
                        logger.warning("Host taxid is -1, not full host")
    except Exception as e:
        logger.error(f"Error processing experiment: {e}")
    return new_experiment

def parse_psi30(files: List[str | Path], interaction_milabel_ok: List[str], mutation_mi_ok: List[str], bindsite_mi_ok: List[str], ptm_mi_ok: List[str], ok_part_types: str | List[str]):
    """
    Parse PSI-MI 3.0 formatted files to get interaction information.
    """

    interaction_label = []
    interaction_mi = []
    interaction_xml_id = [] # from <interaction id="...">" line in XML file
    # There can be multiple experiments per interaction 
    # Each experiment entry will be a list of dictionaries
    experiments = []
    year = []
    interaction_intactid = []
    process_method = []
    
    gene_name_1 = []; gene_symbol_1 = []; length_1 = []; protein_1 = []; uniprotkb_1 = []; species_label_1 = []; species_taxid_1 = []
    ensp_1 = []; ensg_1 = []; enst_1 = []; interpro_1 = []; reactome_1 = []; rscbpdb_1 = []; intactid_1 = []; dip_1 = []
    primaryref_db_1 = []; primaryref_id_1 = []; go_1 = []; host_taxid_1 = []; host_cell_type_1 = []
    host_compartment_1 = []; host_tissue_1 = []; host_label_short_1 = []; host_label_full_1 = []
    mol_type_1 = []

    gene_name_2 = []; gene_symbol_2 = []; length_2 = []; protein_2 = []; uniprotkb_2 = []; species_label_2 = []; species_taxid_2 = []
    ensp_2 = []; ensg_2 = []; enst_2 = []; interpro_2 = []; reactome_2 = []; rscbpdb_2 = []; intactid_2 = []; dip_2 = []
    primaryref_db_2 = []; primaryref_id_2 = []; go_2 = []; host_taxid_2 = []; host_cell_type_2 = []
    host_compartment_2 = []; host_tissue_2 = []; host_label_short_2 = []; host_label_full_2 = []
    mol_type_2 = []
    
    # features (per participant)
    mutation_mi_1 = []; mutation_name_1 = []; mutation_short_1 = []
    mutation_begin_1 = []; mutation_end_1 = []
    mutation_orig_1 = []; mutation_new_1 = []
    binding_mi_1 = [];  binding_name_1 = []; binding_short_1 = []
    binding_begin_1 = []; binding_end_1 = []
    chain_seq_start_1 = []; chain_seq_end_1 = []

    mutation_mi_2 = []; mutation_name_2 = []; mutation_short_2 = []
    mutation_begin_2 = []; mutation_end_2 = []
    mutation_orig_2 = []; mutation_new_2 = []
    binding_mi_2 = [];  binding_name_2 = []; binding_short_2 = []
    binding_begin_2 = []; binding_end_2 = []
    chain_seq_start_2 = []; chain_seq_end_2 = []
    
    ptm_mi_1 = []; ptm_name_1 = []; ptm_begin_1 = []; ptm_end_1 = []; ptm_orig_1 = []; ptm_new_1 = []; ptm_short_1 = []
    ptm_mi_2 = []; ptm_name_2 = []; ptm_begin_2 = []; ptm_end_2 = []; ptm_orig_2 = []; ptm_new_2 = []; ptm_short_2 = []
    
    # SKIPPED INTERACTION LOG below - track things we don't include and why
    skipped_interactions = []
    
    # normalize acceptable labels for robust matching
    interaction_mi_norm = {str(x).strip().lower() for x in interaction_milabel_ok if pd.notna(x)}
    
    def _log_skip(file_path, node, reason):
        """Record an interaction we decided not to fully process."""
        p = Path(file_path).resolve()
        try:
            # Strip off the project root prefix; this will NEVER produce "../"
            rel_file = p.relative_to(ROOT_PATH)
        except ValueError:
            # If the file isn’t under root for some reason, fall back gracefully
            rel_file = p.name  # or just 'p' if you prefer the full path

        xml_id = _xml_attr(node, "id")
        intact_ids = _xml_interaction_all_intact_ids(node)  # may be None
        skipped_interactions.append({
            "file": rel_file,
            "interaction_xml_id": xml_id,
            "interaction_intactid": ",".join(intact_ids) if intact_ids else None,
            "reason": reason,
        })
        
    for file in files:
        # get the year. Example file path ends in psi30/pmid/2003/14609943.xml
        try:
            fyear = int(file.split("psi30/pmid/")[-1].split("/")[0])
        except:
            fyear = None
        
        # Build XML mirrors once per file (for fallback)
        try:
            tree = ET.parse(file)
            root = tree.getroot()
            xml_inter_nodes, xml_inter_by_intact = _build_interaction_xml_index(root)
            xml_inter_map   = _build_interactor_map_xml(root)
            xml_exp_map     = _build_experiment_map_xml(root)
        except Exception as e:
            logger.warning(f"XML parse fallback disabled for {file}: {e}")
            xml_inter_nodes = []
            xml_inter_map   = {}
            xml_exp_map     = {}

        n = len(xml_inter_nodes)
        
        for i in range(n):            
            # ---------- XML path (binary protein-protein) ----------
            if i >= len(xml_inter_nodes):
                # no XML node to fallback to
                continue

            node = xml_inter_nodes[i]
            xml_id = _xml_attr(node, "id")

            # interaction type label + MI code
            label_xml = _xml_interaction_label(node)
            mi_code = _xml_interaction_mi(node)
                    
            # gate by acceptable labels
            if (label_xml or "").strip().lower() not in interaction_mi_norm:
                continue

            # participants (resolve refs or inline)
            plist = next((c for c in node if _S(c.tag) == "participantList"), None)
            part_parent = plist if plist is not None else node

            xml_parts = []
            for p in part_parent:
                if _S(p.tag) != "participant":
                    continue

                # ref
                info = None
                
                ir = p.find(".//{*}interactorRef")
                if ir is not None and ir.text:
                    info = xml_inter_map.get(ir.text)
                    
                if info is None:
                    inter_el = p.find(".//{*}interactor")
                    if inter_el is not None:
                        iid = _xml_attr(inter_el, "id")
                        if iid and iid in xml_inter_map:
                            info = xml_inter_map[iid]
                        else:
                            # # synthesize a dict for inline (db/acc/short/type/sequence)
                            names = inter_el.find(".//{*}names")
                            sl = names.find(".//{*}shortLabel") if names is not None else None
                            short = _xml_text(sl)
                            
                            # gene_symbol from alias(typeAc=MI:0301) or type="gene name"
                            gene_symbol = None
                            if names is not None:
                                for al in names.findall(".//{*}alias"):
                                    t_ac  = _xml_attr(al, "typeAc")
                                    t_txt = _xml_attr(al, "type")
                                    if (t_ac == "MI:0301") or (t_txt and t_txt.strip().lower() == "gene name"):
                                        gene_symbol = _xml_text(al)
                                        break
                
                            # xref
                            xref = inter_el.find(".//{*}xref")
                            pref = xref.find(".//{*}primaryRef") if xref is not None else None
                            db  = _xml_attr(pref, "db")
                            acc = _xml_attr(pref, "id") or _xml_attr(pref, "ac")
                            
                            ensp = []; ensg = []; enst = []; interpro = []; rscbpdb = []; go = []; intact_ids = []; dip = []
                            def _push(dbv, acv):
                                if not acv: return
                                if dbv == "ensembl":
                                    if   acv.startswith("ENSP"): ensp.append(acv)
                                    elif acv.startswith("ENSG"): ensg.append(acv)
                                    elif acv.startswith("ENST"): enst.append(acv)
                                elif dbv == "interpro":     interpro.append(acv)
                                elif dbv == "rscb pdb":     rscbpdb.append(acv)
                                elif dbv == "go":           go.append(acv)
                                elif dbv == "intact":       intact_ids.append(acv)
                                elif dbv == "dip":          dip.append(acv)

                            # include primaryRef too
                            _push(db, acc)
                            # include all secondaryRef
                            if xref is not None:
                                for c in xref:
                                    if _S(c.tag) != "secondaryRef": continue
                                    _push(_xml_attr(c, "db"), _xml_attr(c, "id") or _xml_attr(c, "ac"))
                
                            acc = _xml_attr(pref, "id") or _xml_attr(pref, "ac")
                            
                            # organism (species), not experimental host
                            org_el = inter_el.find(".//{*}organism")
                            species_taxid = _xml_attr(org_el, "ncbiTaxId") if org_el is not None else None
                            species_label = None
                            if org_el is not None:
                                on = org_el.find(".//{*}names")
                                osl = on.find(".//{*}shortLabel") if on is not None else None
                                species_label = _xml_text(osl)

                            # type
                            itype2 = inter_el.find(".//{*}interactorType")
                            it_names = itype2.find(".//{*}names") if itype2 is not None else None
                            it_sl = it_names.find(".//{*}shortLabel") if it_names is not None else None
                            type_label = _xml_text(it_sl)
                            # sequence
                            seq_el = inter_el.find(".//{*}sequence")
                            seq = _xml_text(seq_el)
                            
                            intact_ids = []
                            if xref is not None:
                                for c in xref:
                                    if _S(c.tag) == "secondaryRef" and _xml_attr(c, "db") == "intact":
                                        intact_ids.append(_xml_attr(c, "id") or _xml_attr(c, "ac"))
                            
                            # attributeList on INLINE interactor (same pattern as in _build_interactor_map_xml)
                            chain_seq_start = None
                            chain_seq_end = None
                            for attr in inter_el.iter():
                                if _S(attr.tag) != "attribute":
                                    continue
                                name = _xml_attr(attr, "name")
                                if name == "chain-seq-start" and _xml_text(attr):
                                    chain_seq_start = _xml_text(attr).strip()
                                elif name == "chain-seq-end" and _xml_text(attr):
                                    chain_seq_end = _xml_text(attr).strip()

                            info = {
                                "db": db, "acc": acc, "short": short, "type": type_label, "sequence": seq,
                                "gene_symbol": gene_symbol,
                                "ensp": ensp, "ensg": ensg, "enst": enst,
                                "interpro": interpro, "rscbpdb": rscbpdb, "go": go,
                                "intact_ids": intact_ids, "dip": dip,
                                "species_taxid": species_taxid, "species_label": species_label,
                                "chain_seq_start": chain_seq_start,
                                "chain_seq_end": chain_seq_end,
                            }
                if info is None:
                    continue
                    
                feat = _extract_features_from_participant(p, mutation_mi_ok, bindsite_mi_ok, ptm_mi_ok)
                info = dict(info)  # copy, just to be safe
                # Do NOT overwrite existing chain_seq_* with None from participant
                for k in ("chain_seq_start", "chain_seq_end"):
                    if feat.get(k) is None:
                        feat.pop(k, None)
                info.update(feat)

                xml_parts.append(info)
                
            # require exactly two protein OR peptide participants
            if len(xml_parts) != 2:
                _log_skip(file, node, "non_binary_participant_count")
                continue
            if not ( (xml_parts[0].get("type") in ok_part_types) and (xml_parts[1].get("type") in ok_part_types) ):
                _log_skip(file, node, "non_protein_peptide_participant_type")
                continue
            
            # OK: append XML-derived row identical to pymex path
            interaction_label.append(label_xml)
            interaction_mi.append(mi_code)
            interaction_xml_id.append(xml_id) 
            year.append(fyear)
            
            # IntAct IDs on interaction xref
            cur_interaction_intactid = []
            xref_i = next((c for c in node if _S(c.tag) == "xref"), None)
            if xref_i is not None:
                pref_i = next((c for c in xref_i if _S(c.tag) == "primaryRef"), None)
                if pref_i is not None and _xml_attr(pref_i, "db") == "intact":
                    cur_interaction_intactid.append(_xml_attr(pref_i, "id") or _xml_attr(pref_i, "ac"))
                for c in xref_i:
                    if _S(c.tag) == "secondaryRef" and _xml_attr(c, "db") == "intact":
                        cur_interaction_intactid.append(_xml_attr(c, "id") or _xml_attr(c, "ac"))
            interaction_intactid.append(",".join([x for x in cur_interaction_intactid if x]) if cur_interaction_intactid else None)

            # interactor infos
            info_1 = _xml_to_interactor_info(xml_parts[0])
            info_2 = _xml_to_interactor_info(xml_parts[1])
            
            # interactor fields
            gene_name_1.append(info_1["gene_name"]);   gene_symbol_1.append(info_1["gene_symbol"])
            length_1.append(info_1["length"]);         protein_1.append(info_1["protein"])
            species_label_1.append(info_1["species_label"]);    species_taxid_1.append(info_1["species_taxid"])
            uniprotkb_1.append(info_1["uniprotkb"]);   ensp_1.append(info_1["ensp"])
            ensg_1.append(info_1["ensg"]);             enst_1.append(info_1["enst"])
            interpro_1.append(info_1["interpro"]);      rscbpdb_1.append(info_1["rscbpdb"])
            reactome_1.append(info_1["reactome"])
            intactid_1.append(info_1["intactid"]);      primaryref_db_1.append(info_1["primaryref_db"])
            dip_1.append(info_1["dip"])
            primaryref_id_1.append(info_1["primaryref_id"]); go_1.append(info_1["go"])
            mol_type_1.append(info_1["mol_type"])

            gene_name_2.append(info_2["gene_name"]);   gene_symbol_2.append(info_2["gene_symbol"])
            length_2.append(info_2["length"]);         protein_2.append(info_2["protein"])
            species_label_2.append(info_2["species_label"]);    species_taxid_2.append(info_2["species_taxid"])
            uniprotkb_2.append(info_2["uniprotkb"]);   ensp_2.append(info_2["ensp"])
            ensg_2.append(info_2["ensg"]);             enst_2.append(info_2["enst"])
            interpro_2.append(info_2["interpro"]);      rscbpdb_2.append(info_2["rscbpdb"])
            reactome_2.append(info_1["reactome"])
            intactid_2.append(info_2["intactid"]);      primaryref_db_2.append(info_2["primaryref_db"])
            dip_2.append(info_2["dip"])
            primaryref_id_2.append(info_2["primaryref_id"]); go_2.append(info_2["go"])
            mol_type_2.append(info_2["mol_type"])

            # feature fields – now *directly* from info_1/info_2
            mutation_mi_1.append(info_1["mutation_mi"])
            mutation_name_1.append(info_1["mutation_name"])
            mutation_begin_1.append(info_1["mutation_begin"])
            mutation_end_1.append(info_1["mutation_end"])
            mutation_orig_1.append(info_1["mutation_orig"])
            mutation_new_1.append(info_1["mutation_new"])
            mutation_short_1.append(info_1.get("mutation_short"))

            binding_mi_1.append(info_1["binding_mi"])
            binding_name_1.append(info_1["binding_name"])
            binding_begin_1.append(info_1["binding_begin"])
            binding_end_1.append(info_1["binding_end"])
            binding_short_1.append(info_1.get("binding_short"))

            chain_seq_start_1.append(info_1["chain_seq_start"])
            chain_seq_end_1.append(info_1["chain_seq_end"])

            ptm_mi_1.append(info_1["ptm_mi"])
            ptm_name_1.append(info_1["ptm_name"])
            ptm_begin_1.append(info_1["ptm_begin"])
            ptm_end_1.append(info_1["ptm_end"])
            ptm_orig_1.append(info_1["ptm_orig"])
            ptm_new_1.append(info_1["ptm_new"])
            ptm_short_1.append(info_1.get("ptm_short"))

            mutation_mi_2.append(info_2["mutation_mi"])
            mutation_name_2.append(info_2["mutation_name"])
            mutation_begin_2.append(info_2["mutation_begin"])
            mutation_end_2.append(info_2["mutation_end"])
            mutation_orig_2.append(info_2["mutation_orig"])
            mutation_new_2.append(info_2["mutation_new"])
            mutation_short_2.append(info_2.get("mutation_short"))

            binding_mi_2.append(info_2["binding_mi"])
            binding_name_2.append(info_2["binding_name"])
            binding_begin_2.append(info_2["binding_begin"])
            binding_end_2.append(info_2["binding_end"])
            binding_short_2.append(info_2.get("binding_short"))

            chain_seq_start_2.append(info_2["chain_seq_start"])
            chain_seq_end_2.append(info_2["chain_seq_end"])

            ptm_mi_2.append(info_2["ptm_mi"])
            ptm_name_2.append(info_2["ptm_name"])
            ptm_begin_2.append(info_2["ptm_begin"])
            ptm_end_2.append(info_2["ptm_end"])
            ptm_orig_2.append(info_2["ptm_orig"])
            ptm_new_2.append(info_2["ptm_new"])
            ptm_short_2.append(info_2.get("ptm_short"))

            process_method.append("xml")
            
            # experiments (XML): either inline <experimentDescription> or <experimentRef>
            cur_experiments = None
            exp_list = next((c for c in node if _S(c.tag) == "experimentList"), None)
            if exp_list is not None:
                cur = []
                for e in exp_list:
                    tag = _S(e.tag)
                    if tag == "experimentDescription":
                        cur.append(_parse_experiment_node_xml(e))
                    elif tag == "experimentRef":
                        refid = _xml_text(e)
                        if refid and refid in xml_exp_map:
                            cur.append(_parse_experiment_node_xml(xml_exp_map[refid]))
                if cur:
                    cur_experiments = cur
            experiments.append(cur_experiments)
            
            host_taxids = []
            host_cell_types = []
            host_compartments = []
            host_tissues = []
            host_labels_short = []
            host_labels_full = []

            if cur_experiments:
                for exp in cur_experiments:
                    for h in exp.get("hosts", []):
                        taxid = h.get("taxid")
                        if taxid:
                            host_taxids.append(taxid)

                        ct = h.get("cell_type")
                        if ct:
                            host_cell_types.append(ct)

                        comp = h.get("compartment")
                        if comp:
                            host_compartments.append(comp)

                        tis = h.get("tissue")
                        if tis:
                            host_tissues.append(tis)

                        ls = h.get("label_short")
                        if ls:
                            host_labels_short.append(ls)

                        lf = h.get("label_full")
                        if lf:
                            host_labels_full.append(lf)

            # de-duplicate while preserving order
            def _dedup(seq):
                out = []
                seen = set()
                for v in seq:
                    if not v:
                        continue
                    if v not in seen:
                        seen.add(v)
                        out.append(v)
                return out

            host_taxid_str       = ",".join(_dedup(host_taxids))       or None
            host_cell_type_str   = ",".join(_dedup(host_cell_types))   or None
            host_compartment_str = ",".join(_dedup(host_compartments)) or None
            host_tissue_str      = ",".join(_dedup(host_tissues))      or None
            host_label_short_str = ",".join(_dedup(host_labels_short)) or None
            host_label_full_str  = ",".join(_dedup(host_labels_full))  or None

            host_taxid_1.append(host_taxid_str)
            host_taxid_2.append(host_taxid_str)

            host_cell_type_1.append(host_cell_type_str)
            host_cell_type_2.append(host_cell_type_str)

            host_compartment_1.append(host_compartment_str)
            host_compartment_2.append(host_compartment_str)

            host_tissue_1.append(host_tissue_str)
            host_tissue_2.append(host_tissue_str)

            host_label_short_1.append(host_label_short_str)
            host_label_short_2.append(host_label_short_str)

            host_label_full_1.append(host_label_full_str)
            host_label_full_2.append(host_label_full_str)
            
    ret_dict = {
        ## Full interaction data
        "interaction_label": interaction_label,
        "interaction_mi": interaction_mi,
        "interaction_intactid": interaction_intactid,
        "interaction_xml_id": interaction_xml_id,
        "experiments": experiments,
        "year": year,
        "process_method": process_method,
        ## Interactor 1
        "gene_name_1": gene_name_1,
        "gene_symbol_1": gene_symbol_1,
        "mol_type_1": mol_type_1,
        "species_label_1": species_label_1,
        "species_taxid_1": species_taxid_1,
        "length_1": length_1,
        "protein_1": protein_1,
        "chain_seq_start_1": chain_seq_start_1,
        "chain_seq_end_1": chain_seq_end_1,
        "uniprotkb_1": uniprotkb_1,
        "ensp_1": ensp_1,
        "ensg_1": ensg_1,
        "enst_1": enst_1,
        "interpro_1": interpro_1,
        "reactome_1": reactome_1,
        "rscbpdb_1": rscbpdb_1,
        "intactid_1": intactid_1,
        "dip_1": dip_1,
        "primaryref_db_1": primaryref_db_1,
        "primaryref_id_1": primaryref_id_1,
        "go_1": go_1,
        "host_taxid_1": host_taxid_1,
        "host_label_short_1": host_label_short_1,
        "host_label_full_1": host_label_full_1,
        "host_cell_type_1": host_cell_type_1,
        "host_compartment_1": host_compartment_1,
        "host_tissue_1": host_tissue_1,
        # --- feature columns ---
        "mutation_mi_1": mutation_mi_1,
        "mutation_name_1": mutation_name_1,
        "mutation_short_1": mutation_short_1,
        "mutation_begin_1": mutation_begin_1,
        "mutation_end_1": mutation_end_1,
        "mutation_orig_1": mutation_orig_1,
        "mutation_new_1": mutation_new_1,
        "binding_mi_1": binding_mi_1,
        "binding_name_1": binding_name_1,
        "binding_short_1": binding_short_1,
        "binding_begin_1": binding_begin_1,
        "binding_end_1": binding_end_1,
        "chain_seq_start_1": chain_seq_start_1,
        "chain_seq_end_1": chain_seq_end_1,
        "ptm_mi_1": ptm_mi_1,
        "ptm_name_1": ptm_name_1,
        "ptm_short_1": ptm_short_1,
        "ptm_begin_1": ptm_begin_1,
        "ptm_end_1": ptm_end_1,
        "ptm_orig_1": ptm_orig_1,
        "ptm_new_1": ptm_new_1,
        ## Interactor 2
        "gene_name_2": gene_name_2,
        "gene_symbol_2": gene_symbol_2,
        "mol_type_2": mol_type_2,
        "species_label_2": species_label_2,
        "species_taxid_2": species_taxid_2,
        "length_2": length_2,
        "protein_2": protein_2,
        "chain_seq_start_2": chain_seq_start_2,
        "chain_seq_end_2": chain_seq_end_2,
        "uniprotkb_2": uniprotkb_2,
        "ensp_2": ensp_2,
        "ensg_2": ensg_2,
        "enst_2": enst_2,
        "interpro_2": interpro_2,
        "reactome_2": reactome_2,
        "rscbpdb_2": rscbpdb_2,
        "intactid_2": intactid_2,
        "dip_2": dip_2,
        "primaryref_db_2": primaryref_db_2,
        "primaryref_id_2": primaryref_id_2,
        "go_2": go_2,
        "host_taxid_2": host_taxid_2,
        "host_label_short_2": host_label_short_1,
        "host_label_full_2": host_label_full_1,
        "host_cell_type_2": host_cell_type_2,
        "host_compartment_2": host_compartment_2,
        "host_tissue_2": host_tissue_2,
        "mutation_mi_2": mutation_mi_2,
        "mutation_name_2": mutation_name_2,
        "mutation_short_2": mutation_short_2,
        "mutation_begin_2": mutation_begin_2,
        "mutation_end_2": mutation_end_2,
        "mutation_orig_2": mutation_orig_2,
        "mutation_new_2": mutation_new_2,
        "binding_mi_2": binding_mi_2,
        "binding_name_2": binding_name_2,
        "binding_short_2": binding_short_2,
        "binding_begin_2": binding_begin_2,
        "binding_end_2": binding_end_2,
        "ptm_mi_2": ptm_mi_2,
        "ptm_name_2": ptm_name_2,
        "ptm_short_2": ptm_short_2,
        "ptm_begin_2": ptm_begin_2,
        "ptm_end_2": ptm_end_2,
        "ptm_orig_2": ptm_orig_2,
        "ptm_new_2": ptm_new_2,
    }
    
    skip_dict = {"skipped_interactions": skipped_interactions}
    
    return ret_dict, skip_dict


def get_pos_neg_files(folder):
    """
    Get positive and negative files from the given path.

    Args:
        path (str): Path to the directory containing files.

    Returns:
        tuple: Two lists, one for positive files and one for negative files.
    """
    posfiles = []
    negfiles = []

    for path, subdirs, files in os.walk(folder):
        for name in files:
            file_name = os.path.join(path, name)
            if "negative" not in file_name:
                posfiles.append(file_name)
            else:
                negfiles.append(file_name)

    return posfiles, negfiles

def parse_and_save(file_idx, file_batch, interaction_milabel_ok, mutation_mi_ok, bindsite_mi_ok, ptm_mi_ok, output_dir, ok_part_types):
    """
    Parses a batch of PSI-MI 3.0 files and saves output to a CSV file.
    """
    psi30_info, skip_info = parse_psi30(file_batch, interaction_milabel_ok, mutation_mi_ok, bindsite_mi_ok, ptm_mi_ok, ok_part_types)

    df = pd.DataFrame.from_dict(psi30_info)
    df = df.rename(columns={
        "gene_name_1": "protein_1",
        "gene_name_2": "protein_2",
        "protein_1": "aa_1",
        "protein_2": "aa_2",
    })

    out_path = os.path.join(output_dir, f"interactome_part_{file_idx}.csv")
    df.to_csv(out_path, index=False)
    
    skip_info = skip_info["skipped_interactions"]
    if len(skip_info)>0:
        skipped_df = pd.DataFrame(skip_info)
        skipped_path = os.path.join(output_dir, f"skipped_interactions_part_{file_idx}.csv")
        skipped_df.to_csv(skipped_path, index=False)             
                              
    return out_path  # Return path for tracking

def split_into_chunks(lst, n_chunks):
    """
    Splits a list into approximately equal-sized chunks.
    """
    chunk_size = ceil(len(lst) / n_chunks)
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]

def run_parallel_parsing(posfiles, interaction_milabel_ok, mutation_mi_ok, bindsite_mi_ok, ptm_mi_ok, output_dir="../data_files/processed/intact/output_chunks", ok_part_types: str | List[str]=["protein","peptide"]):
    os.makedirs(output_dir, exist_ok=True)

    n_cores = multiprocessing.cpu_count()
    # for safety, leave one out
    n_cores = n_cores-1
    file_chunks = split_into_chunks(posfiles, n_cores)

    with multiprocessing.Pool(processes=n_cores) as pool:
        tasks = [
            (i, chunk)
            for i, chunk in enumerate(file_chunks)
        ]
        results = pool.starmap(
            partial(parse_and_save, 
                    interaction_milabel_ok=interaction_milabel_ok, 
                    mutation_mi_ok=mutation_mi_ok,
                    bindsite_mi_ok=bindsite_mi_ok,
                    ptm_mi_ok=ptm_mi_ok,
                    output_dir=output_dir,
                    ok_part_types=ok_part_types),
            tasks
        )
    return results

def merge_outputs(intermediate_output_dir: Path, interaction_type: str, final_output_dir: Path):
    dtypes = {
        # amino acid sequences
        "aa_1": "string",
        "aa_2": "string",

        # feature positions (stored as strings like "10|20" or "10,20|30")
        "binding_begin_1": "string",
        "binding_begin_2": "string",
        "binding_end_1": "string",
        "binding_end_2": "string",
        "mutation_begin_1": "string",
        "mutation_begin_2": "string",
        "mutation_end_1": "string",
        "mutation_end_2": "string",
        "chain_seq_start_1": "string",
        "chain_seq_start_2": "string",
        "chain_seq_end_1": "string",
        "chain_seq_end_2": "string",
        "ptm_begin_1": "string",
        "ptm_begin_2": "string",
        "ptm_end_1": "string",
        "ptm_end_2": "string",

        # feature MI terms / names / residues (all pipe-joined strings)
        "binding_mi_1": "string",
        "binding_mi_2": "string",
        "binding_name_1": "string",
        "binding_name_2": "string",
        "mutation_mi_1": "string",
        "mutation_mi_2": "string",
        "mutation_name_1": "string",
        "mutation_name_2": "string",
        "mutation_new_1": "string",
        "mutation_new_2": "string",
        "mutation_orig_1": "string",
        "mutation_orig_2": "string",
        "ptm_mi_1": "string",
        "ptm_mi_2": "string",
        "ptm_name_1": "string",
        "ptm_name_2": "string",
        "ptm_new_1": "string",
        "ptm_new_2": "string",
        "ptm_orig_1": "string",
        "ptm_orig_2": "string",
        
        # short names of features
        "mutation_short_1": "string",
        "mutation_short_2": "string",
        "binding_short_1": "string",
        "binding_short_2": "string",
        "ptm_short_1": "string",
        "ptm_short_2": "string",

        # gene symbols, types, etc.
        "gene_symbol_1": "string",
        "gene_symbol_2": "string",
        "mol_type_1": "string",
        "mol_type_2": "string",

        # sequences (same as aa_1/aa_2 before renaming)
        "protein_1": "string",
        "protein_2": "string",

        # lengths (true numeric, but nullable)
        "length_1": "Int64",
        "length_2": "Int64",

        # ID lists (mostly comma-separated)
        "ensg_1": "string",
        "ensg_2": "string",
        "ensp_1": "string",
        "ensp_2": "string",
        "enst_1": "string",
        "enst_2": "string",
        "go_1": "string",
        "go_2": "string",
        "dip_1": "string",
        "dip_2": "string",
        "interpro_1": "string",
        "interpro_2": "string",
        "intactid_1": "string",
        "intactid_2": "string",
        "interaction_intactid": "string",
        "rscbpdb_1": "string",
        "rscbpdb_2": "string",
        "uniprotkb_1": "string",
        "uniprotkb_2": "string",
        "reactome_1": "string",
        "reactome_2": "string",
        "species_taxid_1": "string",
        "species_taxid_2": "string",

        # species / host labels
        "species_label_1": "string",
        "species_label_2": "string",
        "host_taxid_1": "string",
        "host_taxid_2": "string",
        "host_cell_type_1": "string",
        "host_cell_type_2": "string",
        "host_compartment_1": "string",
        "host_compartment_2": "string",
        "host_tissue_1": "string",
        "host_tissue_2": "string",
        "host_label_full_1": "string",
        "host_label_full_2": "string",
        "host_label_short_1": "string",
        "host_label_short_2": "string",

        # interaction-level info
        "interaction_label": "string",
        "interaction_mi": "string",
        "process_method": "string",

        # primary ref info
        "primaryref_db_1": "string",
        "primaryref_db_2": "string",
        "primaryref_id_1": "string",
        "primaryref_id_2": "string",

        # experiment blobs (list-of-dicts as repr)
        "experiments": "string",

        # miscellaneous
        "go_1": "string",
        "go_2": "string",

        # year (numeric but may be missing)
        "year": "Int64",
    }
    dfs = []

    intermediate_output_dir = Path(intermediate_output_dir)
    final_output_dir = Path(final_output_dir)

    for file in sorted(intermediate_output_dir.glob("interactome_part_*.csv")):
        df = pd.read_csv(file, dtype=dtypes)
        dfs.append(df)

    if not dfs:
        logger.warning(f"No CSV files found in {intermediate_output_dir}. Skipping merge.")
        return

    merged = pd.concat(dfs, ignore_index=True)

    today = datetime.now().strftime("%Y-%m-%d")
    merged_filename = final_output_dir / f"intact_processed_{interaction_type}PPIs_{today}.csv"
    final_output_dir.mkdir(parents=True, exist_ok=True)

    merged.to_csv(merged_filename, index=False)
    logger.info(f"Merged output written to {merged_filename}")
    
def merge_skipped_outputs(intermediate_output_dir: Path, interaction_type: str, final_output_dir: Path):
    intermediate_output_dir = Path(intermediate_output_dir)
    final_output_dir = Path(final_output_dir)

    dfs = []
    for file in sorted(intermediate_output_dir.glob("skipped_interactions_part_*.csv")):
        df = pd.read_csv(file, dtype="string")
        dfs.append(df)

    if not dfs:
        logger.info(f"No skipped interaction files found in {intermediate_output_dir}.")
        return

    merged = pd.concat(dfs, ignore_index=True)
    today = datetime.now().strftime("%Y-%m-%d")
    out = final_output_dir / f"intact_skipped_{interaction_type}PPIs_{today}.csv"
    final_output_dir.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, index=False)
    logger.info(f"Merged skipped-interaction log written to {out}")


def main(cfg: DictConfig):
    """
    Main method. Process interactome data.
    """
    mode = cfg.process.mode
    logger.info(f"Running IntAct processor in mode: {mode}")

    intact_folder = Path(root) / cfg.process.input_dir
    mypath_pmid = intact_folder / "psi30/pmid"
    mypaths_terms = {}
    for labeltype in cfg.process.label_term_files:
        curpath = (Path(root) / cfg.process.label_term_files[labeltype]).resolve()
        fname = str(curpath.name if isinstance(curpath, Path) else os.path.basename(curpath))
        if "mi_" in fname and "_subtree" in fname:
            curmi = "MI:" + fname.split("mi_")[1].split("_subtree")[0]
            if not curpath.exists():
                logger.info("Generating controlled vocabulary terms...")
                organize_cv(starting_id=curmi)
        else:
            # optional: raise/log
            logger.warning(f"Could not infer MI root from {fname}")
        mypaths_terms[labeltype] = curpath
            
    pos_output_dir = Path(root) / cfg.process.pos_output_dir
    neg_output_dir = Path(root) / cfg.process.neg_output_dir
    example_output_dir = Path(root) / cfg.process.example_output_dir

    posfiles, negfiles = get_pos_neg_files(mypath_pmid)
    logger.info(
        f"Total positive files: {len(posfiles)}\nTotal negative files: {len(negfiles)}"
    )

    interaction_milabel_ok = pd.read_csv(mypaths_terms["interaction"])["label"].tolist()
    mutation_mi_ok = pd.read_csv(mypaths_terms["mutation"])["id"].tolist()
    bindsite_mi_ok = pd.read_csv(mypaths_terms["binding_site"])["id"].tolist()
    ptm_mi_ok = pd.read_csv(mypaths_terms["observed_ptm"])["id"].tolist()
    
    # If running in debug mode, just do the debug 
    if cfg.process.debug:
        run_debug_mode(
            debug_file=str(Path(root) / cfg.process.debug_file),
            interaction_milabel_ok=interaction_milabel_ok,
            mutation_mi_ok=mutation_mi_ok,
            bindsite_mi_ok=bindsite_mi_ok,
            ptm_mi_ok=ptm_mi_ok,
            output_dir=Path(HydraConfig.get().runtime.output_dir)
        )
        return

    ok_part_types = cfg.process.ok_participant_types
    if isinstance(ok_part_types, str):
        ok_part_types = [ok_part_types]
    ok_part_types = list(set(ok_part_types))
    
    if mode in ["positive", "all"]:
        logger.info("Running sample positive file")
        example_fname = Path(posfiles[0]).stem
        test_info, test_skip_info = parse_psi30([posfiles[0]], interaction_milabel_ok, mutation_mi_ok, bindsite_mi_ok, ptm_mi_ok, ok_part_types)
        df = pd.DataFrame.from_dict(test_info).rename(
            columns={
                "gene_name_1": "protein_1",
                "gene_name_2": "protein_2",
                "protein_1": "aa_1",
                "protein_2": "aa_2",
            }
        )
        example_output_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(example_output_dir / f"positive_{example_fname}_interactome.csv", index=False)
        
        test_skip_info = test_skip_info["skipped_interactions"]
        if len(test_skip_info)>0:
            skipped_df = pd.DataFrame(test_skip_info)
            skipped_path = example_output_dir / f"positive_{example_fname}_skipped_interactions.csv"
            skipped_df.to_csv(skipped_path, index=False)

        logger.info("Running parallel processing for positives")
        run_parallel_parsing(posfiles, 
                             interaction_milabel_ok, 
                             mutation_mi_ok,
                             bindsite_mi_ok,
                             ptm_mi_ok,
                             output_dir=pos_output_dir,
                             ok_part_types=ok_part_types)
        merge_outputs(
            intermediate_output_dir=pos_output_dir,
            interaction_type="positive",
            final_output_dir=Path(root) / cfg.process.final_output_dir
        )
        
        merge_skipped_outputs(
            intermediate_output_dir=pos_output_dir,
            interaction_type="positive",
            final_output_dir=Path(root) / cfg.process.final_output_dir,
        )

    if mode in ["negative", "all"]:
        logger.info("Running parallel processing for negatives")
        run_parallel_parsing(negfiles, 
                             interaction_milabel_ok, 
                             mutation_mi_ok,
                             bindsite_mi_ok,
                             ptm_mi_ok, ok_part_types=ok_part_types, output_dir=neg_output_dir)
        merge_outputs(
            intermediate_output_dir=neg_output_dir,
            interaction_type="negative",
            final_output_dir=Path(root) / cfg.process.final_output_dir
        )

if __name__ == "__main__":
    main()
