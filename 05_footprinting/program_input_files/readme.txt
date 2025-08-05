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



########## Motif DNA Binding Domain annoations #######################

#Rscript


```{r, warning=FALSE}
suppressPackageStartupMessages({
  library(TFBSTools)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(gridExtra)
  library(knitr)
  library(RColorBrewer)
  library(monaLisa)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(seqLogo)
})


jaspar_pcm <- TFBSTools::readJASPARMatrix("../../05_footprinting/program_input_files/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt", matrixClass=c("PFM"))
names(jaspar_pcm) = sapply(jaspar_pcm, function(x) paste0(toupper(x@name),"_", x@ID))

jaspar_ppm <- TFBSTools::toPWM(jaspar_pcm, pseudocounts=0.8, type = "prob")

jaspar_pwm <- TFBSTools::toPWM(jaspar_pcm, pseudocounts=0.8, type = "log2probratio")

jaspar_icm <- toICM(jaspar_pcm, pseudocounts=0.8, schneider=TRUE)



mylist = list()

for (tf in names(jaspar_icm)){
  icm_mat = jaspar_icm[[tf]]@profileMatrix
  pcm_mat = jaspar_pcm[[tf]]@profileMatrix
  ppm_mat = jaspar_ppm[[tf]]@profileMatrix
  
  n_observed = colSums(pcm_mat)[1]
  n_len = ncol(icm_mat)
    
  h = sum(ppm_mat * log2((ppm_mat/0.25)))
  if (n_observed <= 30){
    e = (2-TFBSTools:::schneider_Hnb_precomputed(n_observed)) * n_len
  }else{
    e = (2-TFBSTools:::schneider_Hnb_approx(n_observed, 2)) * n_len
  }
  
  mylist[[tf]] = c(sum(icm_mat), #ic 
                   n_len, #len
                   n_observed, #obs
                   h-e, #ic by hand. for no correction remove e.
                   mean(ppm_mat["G",] + ppm_mat["C",]),
                   max(mean(ppm_mat["A",]), mean(ppm_mat["G",]), 
                       mean(ppm_mat["C",]), mean(ppm_mat["T",])))
}


ic_df = as.data.frame(do.call("rbind",mylist))
colnames(ic_df) = c("IC","len","obs","IC_byHand","GC_content", "Repeat_content")
ic_df$ic_len_norm = ic_df$IC/ic_df$len #check these are below 2. 
ic_df$TF = rownames(ic_df)


### Calculate genome prevalance of chr10. 
seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, "chr10")

res <- findMotifHits(query = jaspar_pwm,
                     subject = seqs,
                     min.score = 6.0, 
                     method = "matchPWM",
                     BPPARAM = BiocParallel::MulticoreParam(workers=10)
                     #BPPARAM = BiocParallel::SerialParam()
                     )

res$fullname = paste0(toupper(res$pwmname), "_", res$pwmid)
tmp_df = data.frame(table(res$fullname))
colnames(tmp_df) = c("TF", "chr10_prevalence")

merged_df <- merge(ic_df, tmp_df, by="TF", all=TRUE) 

write.csv(merged_df[, c("TF", "IC", "len", "obs", "GC_content", "Repeat_content", "chr10_prevalence")], file="jaspar2022_ICcontent.csv", quote=F, row.names = F)
```



########## RNA expression #######################
#bash part
curl -o HEPG2_RNAExpression.tsv "https://www.encodeproject.org/report.tsv?type=RNAExpression&file.assay_title=polyA+plus+RNA-seq&file.biosample_ontology.classification=cell+line&file.assembly=GRCh38&dataset.replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&file.biosample_ontology.term_name=HepG2"
curl -o GM12878_RNAExpression.tsv "https://www.encodeproject.org/report.tsv?type=RNAExpression&file.assay_title=polyA+plus+RNA-seq&file.biosample_ontology.classification=cell+line&file.assembly=GRCh38&dataset.replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&file.biosample_ontology.term_name=GM12878"
curl -o K562_RNAExpression.tsv "https://www.encodeproject.org/report.tsv?type=RNAExpression&file.assay_title=polyA+plus+RNA-seq&file.biosample_ontology.classification=cell+line&file.assembly=GRCh38&dataset.replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&file.biosample_ontology.term_name=K562"
curl -o SKNSH_RNAExpression.tsv "https://www.encodeproject.org/report.tsv?file.biosample_ontology.term_name=SK-N-SH&type=RNAExpression&file.assay_title=polyA+plus+RNA-seq&dataset.replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&file.assembly=GRCh38"
curl -o MCF7_RNAExpression.tsv "https://www.encodeproject.org/report.tsv?type=RNAExpression&file.assembly=GRCh38&dataset.replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&file.biosample_ontology.term_name=MCF-7&file.assay_title=polyA+plus+RNA-seq"


#python part. not totally needed, but easiest. 
df_list = []
for cl in ["HEPG2", "K562", "GM12878", "MCF7", "SKNSH"]:
    infile = f"{cl}_RNAExpression.tsv"
    df = pd.read_csv(infile, sep="\t") #, usecols = [1, 2, 3, 4])
    df = df[df["Assay title"] == "polyA plus RNA-seq"]
    df = df[df['Gene symbol'].notna()]
    tmp = df[["Gene symbol", "TPM"]].groupby('Gene symbol').mean()
    tmp.reset_index(inplace=True)
    tmp.columns=["GeneSymbol", cl]
    print(tmp.shape); print(cl)
    df_list.append(tmp)

df_merged = reduce(lambda  left, right: pd.merge(left, right, on="GeneSymbol", how='outer'), df_list)
df_merged.to_csv(f"/pollard/data/projects/aseveritt/encode_snatacseq/09_simPWMs/ENCODE_bulkRNA/ENCODE_polyAplus.csv", index=False)
