#PART 1
#remove reads that didn't make it into the fragment file. 
#input: *_rep*_unfiltered.bam & *_fragments_sinto.tsv_readnames.txt
#output: *_sinto.bam
######################################################################
samtools view -N 01_sinto/CALU3_rep1_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/CALU3_sinto.bam 00_unfiltered_bams/CALU3_rep1_unfiltered.bam

samtools view -N 01_sinto/IMR90_rep1_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/IMR90_sinto.bam 00_unfiltered_bams/IMR90_rep1_unfiltered.bam

samtools view -N 01_sinto/MCF10A_rep1_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/MCF10A_sinto.bam 00_unfiltered_bams/MCF10A_rep1_unfiltered.bam

samtools view -N 01_sinto/OCILY7_rep1_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/OCILY7_sinto.bam 00_unfiltered_bams/OCILY7_rep1_unfiltered.bam

samtools view -N 01_sinto/PC9_rep1_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/PC9_sinto.bam 00_unfiltered_bams/PC9_rep1_unfiltered.bam

samtools view -N 01_sinto/A549_rep1_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/A549_rep1_sinto.bam 00_unfiltered_bams/A549_rep1_unfiltered.bam
samtools view -N 01_sinto/A549_rep2_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/A549_rep2_sinto.bam 00_unfiltered_bams/A549_rep2_unfiltered.bam
samtools merge 03_filtered_bams/A549_sinto.bam 03_filtered_bams/A549_rep1_sinto.bam 03_filtered_bams/A549_rep2_sinto.bam

samtools view -@ 3 -N 01_sinto/GM12878_rep1_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/GM12878_rep1_sinto.bam 00_unfiltered_bams/GM12878_rep1_unfiltered.bam
samtools view -@ 3 -N 01_sinto/GM12878_rep2_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/GM12878_rep2_sinto.bam 00_unfiltered_bams/GM12878_rep2_unfiltered.bam
samtools view -@ 3 -N 01_sinto/GM12878_rep3_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/GM12878_rep3_sinto.bam 00_unfiltered_bams/GM12878_rep3_unfiltered.bam
samtools merge 03_filtered_bams/GM12878_sinto.bam 03_filtered_bams/GM12878_rep1_sinto.bam 03_filtered_bams/GM12878_rep2_sinto.bam 03_filtered_bams/GM12878_rep3_sinto.bam

samtools view -@ 5 -N 01_sinto/HEPG2_rep1_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/HEPG2_rep1_sinto.bam 00_unfiltered_bams/HEPG2_rep1_unfiltered.bam
samtools view -@ 5 -N 01_sinto/HEPG2_rep2_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/HEPG2_rep2_sinto.bam 00_unfiltered_bams/HEPG2_rep2_unfiltered.bam
samtools view -@ 5 -N 01_sinto/HEPG2_rep3_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/HEPG2_rep3_sinto.bam 00_unfiltered_bams/HEPG2_rep3_unfiltered.bam
samtools merge 03_filtered_bams/HEPG2_sinto.bam 03_filtered_bams/HEPG2_rep1_sinto.bam 03_filtered_bams/HEPG2_rep2_sinto.bam 03_filtered_bams/HEPG2_rep3_sinto.bam

samtools view -@ 5 -N 01_sinto/K562_rep1_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/K562_rep1_sinto.bam 00_unfiltered_bams/K562_rep1_unfiltered.bam
samtools view -@ 5 -N 01_sinto/K562_rep2_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/K562_rep2_sinto.bam 00_unfiltered_bams/K562_rep2_unfiltered.bam
samtools view -@ 5 -N 01_sinto/K562_rep3_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/K562_rep3_sinto.bam 00_unfiltered_bams/K562_rep3_unfiltered.bam
samtools merge 03_filtered_bams/K562_sinto.bam 03_filtered_bams/K562_rep1_sinto.bam 03_filtered_bams/K562_rep2_sinto.bam 03_filtered_bams/K562_rep3_sinto.bam

samtools view -@ 5 -N 01_sinto/MCF7_rep1_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/MCF7_rep1_sinto.bam 00_unfiltered_bams/MCF7_rep1_unfiltered.bam
samtools view -@ 5 -N 01_sinto/MCF7_rep2_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/MCF7_rep2_sinto.bam 00_unfiltered_bams/MCF7_rep2_unfiltered.bam
samtools view -@ 5 -N 01_sinto/MCF7_rep3_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/MCF7_rep3_sinto.bam 00_unfiltered_bams/MCF7_rep3_unfiltered.bam
samtools merge 03_filtered_bams/MCF7_sinto.bam 03_filtered_bams/MCF7_rep1_sinto.bam 03_filtered_bams/MCF7_rep2_sinto.bam 03_filtered_bams/MCF7_rep3_sinto.bam

samtools view -@ 5 -N 01_sinto/SKNSH_rep1_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/SKNSH_rep1_sinto.bam 00_unfiltered_bams/SKNSH_rep1_unfiltered.bam
samtools view -@ 5 -N 01_sinto/SKNSH_rep2_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/SKNSH_rep2_sinto.bam 00_unfiltered_bams/SKNSH_rep2_unfiltered.bam
samtools view -@ 5 -N 01_sinto/SKNSH_rep3_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/SKNSH_rep3_sinto.bam 00_unfiltered_bams/SKNSH_rep3_unfiltered.bam
samtools view -@ 5 -N 01_sinto/SKNSH_rep4_fragments_sinto.tsv_readnames.txt -o 03_filtered_bams/SKNSH_rep4_sinto.bam 00_unfiltered_bams/SKNSH_rep4_unfiltered.bam
samtools merge 03_filtered_bams/SKNSH_sinto.bam 03_filtered_bams/SKNSH_rep1_sinto.bam 03_filtered_bams/SKNSH_rep2_sinto.bam 03_filtered_bams/SKNSH_rep3_sinto.bam 03_filtered_bams/SKNSH_rep4_sinto.bam


######################################################################
#PART 2
#remove non chr1-22, X, Y reads. 
#and samtools index them all 
#input: *_sinto.bam
#output: *sinto_chr.bam

chrs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr21 chr22 chrX chrY"
samtools view -b 03_filtered_bams/CALU3_sinto.bam $chrs | samtools sort - -o 03_filtered_bams/CALU3_sinto_chr.bam
samtools view -b 03_filtered_bams/OCILY7_sinto.bam $chrs | samtools sort - -o 03_filtered_bams/OCILY7_sinto_chr.bam
samtools view -b 03_filtered_bams/MCF10A_sinto.bam $chrs | samtools sort - -o 03_filtered_bams/MCF10A_sinto_chr.bam
samtools view -b 03_filtered_bams/IMR90_sinto.bam $chrs | samtools sort - -o 03_filtered_bams/IMR90_sinto_chr.bam
samtools view -b 03_filtered_bams/PC9_sinto.bam $chrs | samtools sort - -o 03_filtered_bams/PC9_sinto_chr.bam
samtools view -b 03_filtered_bams/GM12878_sinto.bam $chrs | samtools sort - -o 03_filtered_bams/GM12878_sinto_chr.bam
samtools view -b 03_filtered_bams/HEPG2_sinto.bam $chrs | samtools sort - -o 03_filtered_bams/HEPG2_sinto_chr.bam
samtools view -b 03_filtered_bams/K562_sinto.bam $chrs | samtools sort - -o 03_filtered_bams/K562_sinto_chr.bam
samtools view -b 03_filtered_bams/MCF7_sinto.bam $chrs | samtools sort - -o 03_filtered_bams/MCF7_sinto_chr.bam
samtools view -b 03_filtered_bams/A549_sinto.bam $chrs | samtools sort - -o 03_filtered_bams/A549_sinto_chr.bam
samtools view -b 03_filtered_bams/SKNSH_sinto.bam $chrs | samtools sort - -o 03_filtered_bams/SKNSH_sinto_chr.bam

cd /pollard/data/projects/aseveritt/encode_snatacseq_bams/03_filtered_bams
for f in `ls *.bam`; do
    samtools index $f
done


######################################################################
#PART 3
#remove CBs that didn't make it out of ArchR. 
#input: *_sinto_chr.bam
#output: *-cellFilt.bam

cd /pollard/data/projects/aseveritt/encode_snatacseq/03_archr/ENCODE_snATAC/CellLineInfo/
names=$(awk -F',' '{print $18}' cellColData.csv | sort | uniq)
for i in $names; do
    awk -F',' '$18 == "'$i'" {print $17, $18"-cellFilt"}' cellColData.csv > $i"_cellColData.txt"
done


samples="A549 CALU3 GM12878 HEPG2 IMR90 K562 MCF10A MCF7 OCILY7 PC9 SKNSH"
for i in $samples; 
do 
cmd="sinto filterbarcodes -b /pollard/data/projects/aseveritt/encode_snatacseq/03_filtered_bams/"$i"_sinto_chr.bam -c /pollard/data/projects/aseveritt/encode_snatacseq/03_archr/ENCODE_snATAC/CellLineInfo/"$i"_cellColData.txt -p 12 --outdir /pollard/data/projects/aseveritt/encode_snatacseq/03_filtered_bams/" 
$cmd;
done


for f in `ls *-cellFilt.bam`; do
    samtools index $f
done

######################################################################



