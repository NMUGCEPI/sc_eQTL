##celescope

##Download the reference sequence
#GRCH37
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
#GRCH38
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

#Build the genome index file
celescope rna mkref \
 --thread 10 \
 --genome_name Homo_sapiens_ensembl_75 \
 --fasta Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
 --gtf Homo_sapiens.GRCh37.75.gtf

#Sample information extraction and quality control
celescope rna sample --outdir .//example/00.sample --sample example --thread 8 --chemistry auto --wells 384  --fq1 ./example_L2_1.fq
#Align sequencing reads to the reference genome, perform cell identification and gene expression quantification, and generate the core expression matrix
celescope rna starsolo --outdir .//example/01.starsolo --sample example --thread 8 --chemistry auto --adapter_3p AAAAAAAAAAAA --genomeDir ./ensembl_75 --outFilterMatchNmin 50 --soloCellFilter "EmptyDrops_CR 30000 0.99 10 45000 90000 500 0.01 20000 0.001 10000" --starMem 32 --soloFeatures "Gene GeneFull_Ex50pAS"  --fq1 ./example_L2_1.fq --fq2 /example_L2_2.fq
#Downstream Analysis and Visualization
celescope rna analysis --outdir .//example/02.analysis --sample example --thread 8 --genomeDir ./ensembl_75  --matrix_file .//example/outs/filtered 


