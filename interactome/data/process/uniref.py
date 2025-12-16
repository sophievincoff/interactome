import gzip
import csv
import xml.etree.ElementTree as ET
from pathlib import Path
import logging
import os
from omegaconf import DictConfig
import rootutils

root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

# If you're not using Hydra logging for this script, uncomment this to be sure logs show up:
# logging.basicConfig(
#     level=logging.INFO,
#     format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
# )

logger = logging.getLogger(__name__)

def stream_uniref_mapping(
    xml_gz_path: str | Path,
    out_tsv_gz_path: str | Path,
    cluster_prefix: str = "",
    max_clusters: int | None = None,
):
    """
    Stream UniRef XML (.xml.gz) and write mapping rows to a gzipped TSV.

    Columns:
      - uniprot_acc    : UniProtKB accession (primary or secondary)
      - uniprot_id     : UniProtKB ID (e.g. FOXO1_HUMAN)
      - is_primary     : 1 if this accession is the primary accession for that ID, else 0
      - uniref_cluster : UniRef cluster ID (e.g. UniRef50_Q9R1E0)

    Logic:
      - For each <entry>, get entry/@id as the cluster ID.
      - For each <dbReference type="UniProtKB ID"> under that entry:
          * Get dbReference/@id        → uniprot_id
          * Get all <property type="UniProtKB accession" value="..."> → list of accessions
          * The *first* accession in that list is treated as primary (is_primary=1),
            all others as secondary (is_primary=0).
      - Emits one row per accession.
      - Includes both representativeMember and member (they share the same structure).
    """

    xml_gz_path = Path(xml_gz_path)
    out_tsv_gz_path = Path(out_tsv_gz_path)
    out_tsv_gz_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Parsing {xml_gz_path} → {out_tsv_gz_path}")

    n_entries = 0
    n_rows = 0

    with gzip.open(xml_gz_path, "rt") as fh, \
         gzip.open(out_tsv_gz_path, "wt", newline="") as out_f:

        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["uniprot_acc", "uniprot_id", "is_primary", "uniref_cluster"])

        for event, elem in ET.iterparse(fh, events=("end",)):
            tag = elem.tag.rsplit("}", 1)[-1]  # strip namespace

            if tag != "entry":
                continue

            n_entries += 1

            # Cluster ID, e.g. "UniRef50_Q9R1E0"
            cluster_id = elem.attrib["id"]
            if cluster_prefix:
                if not cluster_id.startswith(cluster_prefix):
                    suffix = cluster_id.split("_", 1)[-1]
                    cluster_id = f"{cluster_prefix}{suffix}"

            seen_acc = set()

            # Look at all UniProtKB entries in this cluster
            for dbref in elem.findall(".//{*}dbReference[@type='UniProtKB ID']"):
                uniprot_id = dbref.attrib.get("id")  # e.g. FOXO1_HUMAN

                # All accessions for this ID (primary + secondary)
                acc_props = dbref.findall(".//{*}property[@type='UniProtKB accession']")
                accessions = [p.attrib.get("value") for p in acc_props if p.attrib.get("value")]

                if not accessions:
                    continue

                # Convention: first accession listed is primary
                for i, acc in enumerate(accessions):
                    if acc in seen_acc:
                        continue
                    is_primary = 1 if i == 0 else 0
                    writer.writerow([acc, uniprot_id, is_primary, cluster_id])
                    seen_acc.add(acc)
                    n_rows += 1

            elem.clear()

            if n_entries % 10_000 == 0:
                logger.info(f"Processed {n_entries:,} clusters, {n_rows:,} mappings...")

            if max_clusters is not None and n_entries >= max_clusters:
                logger.info(
                    f"Reached max_clusters={max_clusters}. "
                    f"Stopping early with {n_rows:,} mappings."
                )
                break

    logger.info(f"Done. Total clusters: {n_entries:,}, total mappings: {n_rows:,}")
    logger.info(f"To see file preview enter: zcat {out_tsv_gz_path} | head -n 20")

def main(cfg: DictConfig):
    """
    Main method. Process UniRef data.
    """

    uniref50_input_path = Path(root) / cfg.process.uniref_50_input_path
    uniref90_input_path = Path(root) / cfg.process.uniref_90_input_path
    uniref50_output_dir = Path(root) / cfg.process.uniref_50_output_dir
    uniref90_output_dir = Path(root) / cfg.process.uniref_90_output_dir

    os.makedirs(uniref50_output_dir, exist_ok=True)
    os.makedirs(uniref90_output_dir, exist_ok=True)

    # ---- DEBUG MODE SETTINGS ----
    # Expecting something like cfg.process.debug: bool in your Hydra config
    debug = getattr(cfg.process, "debug", False)
    max_clusters = 10_000 if debug else None
    uniref50_output_dir = uniref50_output_dir/"debug" if debug else uniref50_output_dir
    uniref90_output_dir = uniref90_output_dir/"debug" if debug else uniref90_output_dir
    suffix = "_debug" if debug else ""

    # Process UniRef90
    stream_uniref_mapping(
        uniref90_input_path,
        uniref90_output_dir / f"uniref90_mapping{suffix}.tsv.gz",
        cluster_prefix="UniRef90_",
        max_clusters=max_clusters,
    )

    # Process UniRef50
    stream_uniref_mapping(
        uniref50_input_path,
        uniref50_output_dir / f"uniref50_mapping{suffix}.tsv.gz",
        cluster_prefix="UniRef50_",
        max_clusters=max_clusters,
    )
