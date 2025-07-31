# required files:
# JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt
# jaspar2022_motifs/
# JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar_nameFIX.txt
# peakRData/
# jaspar2022_motifs_to_clusters.txt
# jaspar2022_dndbd.txt
# USF_list.csv

########## TOBIAS single file ########################
#bash
wget https://jaspar2022.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt



########## HINT individual files ######################
#bash
python ~/rgtdata/motifs/createPwm.py -i \
JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
-f jaspar-2016 \
-o jaspar2022_motifs



######### PRINT single file #######################
#PRINT requires single file motif but can strangely chop off motif header so:
#python
import re
with open("JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt", "r") as f:
    lines = f.readlines()

new_lines=[]
for i, line in enumerate(lines):
    if line.startswith(">"):
        n = line.split()
        motif_id = re.sub(">" ,"", n[0])
        motif_name = n[1]
        new_str = f">{motif_name}_{motif_id}\t{motif_name}_{motif_id}\n"
        new_lines.append(new_str)
    else:
        new_lines.append(line) 

with open("JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar_nameFIX.txt", "w") as f:
    f.writelines(new_lines)


######### PRINT input RData #######################
see /peakRData/print_peaksub.sh



######### Motif cluster info #########################
wget https://jaspar2022.genereg.net/static/clustering/2022/vertebrates/CORE/interactive_trees/JASPAR_2022_matrix_clustering_vertebrates_archive.zip

#navigate to 
/interactive_trees/JASPAR_2022_matrix_clustering_vertebrates_CORE_tables/alignment_table.tab

awk '{print $1, $2, $3}' alignment_table.tab | awk -F"_" '{print $6"."$7_$8}' | awk '{print $2"_"$1, $3}' > /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/jaspar2022_motifs_to_clusters.txt



########## Motif DNA Binding Domain annoations #######################
#python

import requests
# Function to get all motifs from the JASPAR vertebrate core set with pagination
def get_all_motifs():
    url_core_set = "http://jaspar.genereg.net/api/v1/matrix/?tax_group=vertebrates&collection=CORE&format=json"
    motifs = []
    while url_core_set:
        response_core_set = requests.get(url_core_set)
        if response_core_set.status_code == 200:
            core_set_data = response_core_set.json()
            motifs.extend(core_set_data['results'])
            url_core_set = core_set_data['next']  # URL for the next page of results
        else:
            print(f"Failed to retrieve core set data: {response_core_set.status_code}")
            break
    return motifs

# Retrieve all motifs
motifs = get_all_motifs()

# Step 2: Fetch details for each motif and extract the DNA binding domain
res = []

for motif in motifs:
    #print(motif)
    jaspar_id = motif['matrix_id']
    url_motif = f"http://jaspar.genereg.net/api/v1/matrix/{jaspar_id}/"
    response_motif = requests.get(url_motif)
    
    if response_motif.status_code == 200:
        motif_data = response_motif.json()

        name = motif_data.get('name', 'NA')
        motifid = motif_data.get('matrix_id', 'NA')
        dbd = motif_data.get('class', 'NA')
        family = motif_data.get('family', 'NA')
        res.append((motifid, name, dbd, family))
    else:
        print(f"Failed to retrieve data for {jaspar_id}: {response_motif.status_code}")

print(res)

with open('JASPAR_dbd.txt', 'w') as file:
    for a,b,c,d in res:
        file.write(f"{a},{b},{c},{d}\n")

