import pymex.mif

rec = pymex.mif.Record()
file = "/scratch/pranamlab/sophie/interactome/interactome/2025-10-02-20-14.xml"

# get the year. Example file path ends in psi30/pmid/2003/14609943.xml
rec.parseMif(file, "mif300")
intact_ids = []
# try:
for i, interaction in enumerate(rec.interactions):
    #print(interaction.xrefs)
    print("PRIMARY")
    print(interaction.primaryRef)
    print("SECONDARY BELOW")
    print(interaction.secondaryRef)
    cur_interact_ids = []
    try:
        if interaction.primaryRef and interaction.primaryRef.db =="intact":
            cur_interact_ids.append(interaction.primaryRef.ac)
        if interaction.secondaryRef:
            for sec in interaction.secondaryRef:
                if sec.db=="intact":
                    cur_interact_ids.append(interaction.secondaryRef.ac)
    except Exception as e:
        print(f"Error getting IntAct ID: {e}")
    
    if len(cur_interact_ids)>0:
        intact_ids.append(cur_interact_ids)
    else:
        intact_ids.append(None)
    
    
print(intact_ids)