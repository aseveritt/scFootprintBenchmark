samtools index A549_rep1_unfiltered.bam
sinto fragments -b A549_rep1_unfiltered.bam -f A549_rep1_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i A549_rep1_fragments_sinto.tsv | bgzip -c > A549_rep1_fragments_sinto.tsv.gz
tabix -p bed A549_rep1_fragments_sinto.tsv.gz
rm A549_rep1_fragments_sinto.tsv

samtools index A549_rep2_unfiltered.bam
sinto fragments -b A549_rep2_unfiltered.bam -f A549_rep2_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i A549_rep2_fragments_sinto.tsv | bgzip -c > A549_rep2_fragments_sinto.tsv.gz
tabix -p bed A549_rep2_fragments_sinto.tsv.gz
rm A549_rep2_fragments_sinto.tsv

samtools index CALU3_rep1_unfiltered.bam
sinto fragments -b CALU3_rep1_unfiltered.bam -f CALU3_rep1_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i CALU3_rep1_fragments_sinto.tsv | bgzip -c > CALU3_rep1_fragments_sinto.tsv.gz
tabix -p bed CALU3_rep1_fragments_sinto.tsv.gz
rm CALU3_rep1_fragments_sinto.tsv

samtools index GM12878_rep1_unfiltered.bam
sinto fragments -b GM12878_rep1_unfiltered.bam -f GM12878_rep1_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i GM12878_rep1_fragments_sinto.tsv | bgzip -c > GM12878_rep1_fragments_sinto.tsv.gz
tabix -p bed GM12878_rep1_fragments_sinto.tsv.gz
rm GM12878_rep1_fragments_sinto.tsv

samtools index GM12878_rep2_unfiltered.bam
sinto fragments -b GM12878_rep2_unfiltered.bam -f GM12878_rep2_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i GM12878_rep2_fragments_sinto.tsv | bgzip -c > GM12878_rep2_fragments_sinto.tsv.gz
tabix -p bed GM12878_rep2_fragments_sinto.tsv.gz
rm GM12878_rep2_fragments_sinto.tsv

samtools index GM12878_rep3_unfiltered.bam
sinto fragments -b GM12878_rep3_unfiltered.bam -f GM12878_rep3_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i GM12878_rep3_fragments_sinto.tsv | bgzip -c > GM12878_rep3_fragments_sinto.tsv.gz
tabix -p bed GM12878_rep3_fragments_sinto.tsv.gz
rm GM12878_rep3_fragments_sinto.tsv

samtools index HEPG2_rep1_unfiltered.bam
sinto fragments -b HEPG2_rep1_unfiltered.bam -f HEPG2_rep1_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i HEPG2_rep1_fragments_sinto.tsv | bgzip -c > HEPG2_rep1_fragments_sinto.tsv.gz
tabix -p bed HEPG2_rep1_fragments_sinto.tsv.gz
rm HEPG2_rep1_fragments_sinto.tsv

samtools index HEPG2_rep2_unfiltered.bam
sinto fragments -b HEPG2_rep2_unfiltered.bam -f HEPG2_rep2_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i HEPG2_rep2_fragments_sinto.tsv | bgzip -c > HEPG2_rep2_fragments_sinto.tsv.gz
tabix -p bed HEPG2_rep2_fragments_sinto.tsv.gz
rm HEPG2_rep2_fragments_sinto.tsv

samtools index HEPG2_rep3_unfiltered.bam
sinto fragments -b HEPG2_rep3_unfiltered.bam -f HEPG2_rep3_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i HEPG2_rep3_fragments_sinto.tsv | bgzip -c > HEPG2_rep3_fragments_sinto.tsv.gz
tabix -p bed HEPG2_rep3_fragments_sinto.tsv.gz
rm HEPG2_rep3_fragments_sinto.tsv

samtools index IMR90_rep1_unfiltered.bam
sinto fragments -b IMR90_rep1_unfiltered.bam -f IMR90_rep1_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i IMR90_rep1_fragments_sinto.tsv | bgzip -c > IMR90_rep1_fragments_sinto.tsv.gz
tabix -p bed IMR90_rep1_fragments_sinto.tsv.gz
rm IMR90_rep1_fragments_sinto.tsv

samtools index K562_rep1_unfiltered.bam
sinto fragments -b K562_rep1_unfiltered.bam -f K562_rep1_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i K562_rep1_fragments_sinto.tsv | bgzip -c > K562_rep1_fragments_sinto.tsv.gz
tabix -p bed K562_rep1_fragments_sinto.tsv.gz
rm K562_rep1_fragments_sinto.tsv

samtools index K562_rep2_unfiltered.bam
sinto fragments -b K562_rep2_unfiltered.bam -f K562_rep2_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i K562_rep2_fragments_sinto.tsv | bgzip -c > K562_rep2_fragments_sinto.tsv.gz
tabix -p bed K562_rep2_fragments_sinto.tsv.gz
rm K562_rep2_fragments_sinto.tsv

samtools index K562_rep3_unfiltered.bam
sinto fragments -b K562_rep3_unfiltered.bam -f K562_rep3_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i K562_rep3_fragments_sinto.tsv | bgzip -c > K562_rep3_fragments_sinto.tsv.gz
tabix -p bed K562_rep3_fragments_sinto.tsv.gz
rm K562_rep3_fragments_sinto.tsv

samtools index MCF10A_rep1_unfiltered.bam
sinto fragments -b MCF10A_rep1_unfiltered.bam -f MCF10A_rep1_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i MCF10A_rep1_fragments_sinto.tsv | bgzip -c > MCF10A_rep1_fragments_sinto.tsv.gz
tabix -p bed MCF10A_rep1_fragments_sinto.tsv.gz
rm MCF10A_rep1_fragments_sinto.tsv

samtools index MCF7_rep1_unfiltered.bam
sinto fragments -b MCF7_rep1_unfiltered.bam -f MCF7_rep1_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i MCF7_rep1_fragments_sinto.tsv | bgzip -c > MCF7_rep1_fragments_sinto.tsv.gz
tabix -p bed MCF7_rep1_fragments_sinto.tsv.gz
rm MCF7_rep1_fragments_sinto.tsv

samtools index MCF7_rep2_unfiltered.bam
sinto fragments -b MCF7_rep2_unfiltered.bam -f MCF7_rep2_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i MCF7_rep2_fragments_sinto.tsv | bgzip -c > MCF7_rep2_fragments_sinto.tsv.gz
tabix -p bed MCF7_rep2_fragments_sinto.tsv.gz
rm MCF7_rep2_fragments_sinto.tsv

samtools index MCF7_rep3_unfiltered.bam
sinto fragments -b MCF7_rep3_unfiltered.bam -f MCF7_rep3_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i MCF7_rep3_fragments_sinto.tsv | bgzip -c > MCF7_rep3_fragments_sinto.tsv.gz
tabix -p bed MCF7_rep3_fragments_sinto.tsv.gz
rm MCF7_rep3_fragments_sinto.tsv

samtools index OCILY7_rep1_unfiltered.bam
sinto fragments -b OCILY7_rep1_unfiltered.bam -f OCILY7_rep1_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i OCILY7_rep1_fragments_sinto.tsv | bgzip -c > OCILY7_rep1_fragments_sinto.tsv.gz
tabix -p bed OCILY7_rep1_fragments_sinto.tsv.gz
rm OCILY7_rep1_fragments_sinto.tsv

samtools index PC9_rep1_unfiltered.bam
sinto fragments -b PC9_rep1_unfiltered.bam -f PC9_rep1_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i PC9_rep1_fragments_sinto.tsv | bgzip -c > PC9_rep1_fragments_sinto.tsv.gz
tabix -p bed PC9_rep1_fragments_sinto.tsv.gz
rm PC9_rep1_fragments_sinto.tsv

samtools index SKNSH_rep1_unfiltered.bam
sinto fragments -b SKNSH_rep1_unfiltered.bam -f SKNSH_rep1_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i SKNSH_rep1_fragments_sinto.tsv | bgzip -c > SKNSH_rep1_fragments_sinto.tsv.gz
tabix -p bed SKNSH_rep1_fragments_sinto.tsv.gz
rm SKNSH_rep1_fragments_sinto.tsv

samtools index SKNSH_rep2_unfiltered.bam
sinto fragments -b SKNSH_rep2_unfiltered.bam -f SKNSH_rep2_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i SKNSH_rep2_fragments_sinto.tsv | bgzip -c > SKNSH_rep2_fragments_sinto.tsv.gz
tabix -p bed SKNSH_rep2_fragments_sinto.tsv.gz
rm SKNSH_rep2_fragments_sinto.tsv

samtools index SKNSH_rep3_unfiltered.bam
sinto fragments -b SKNSH_rep3_unfiltered.bam -f SKNSH_rep3_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i SKNSH_rep3_fragments_sinto.tsv | bgzip -c > SKNSH_rep3_fragments_sinto.tsv.gz
tabix -p bed SKNSH_rep3_fragments_sinto.tsv.gz
rm SKNSH_rep3_fragments_sinto.tsv

samtools index SKNSH_rep4_unfiltered.bam
sinto fragments -b SKNSH_rep4_unfiltered.bam -f SKNSH_rep4_fragments_sinto.tsv -p 5 --return_readnames
bedtools sort -i SKNSH_rep4_fragments_sinto.tsv | bgzip -c > SKNSH_rep4_fragments_sinto.tsv.gz
tabix -p bed SKNSH_rep4_fragments_sinto.tsv.gz
rm SKNSH_rep4_fragments_sinto.tsv