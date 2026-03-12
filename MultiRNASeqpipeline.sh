#!/bin/bash

#This is a RNASeq Reads-to-Counts Pipeline
#Use this for only single-end files (one file per sample)
#See MultiRNASeqpipelinePaired.sh for paired-end files

SECONDS=0
set -e

cd /mnt/d/random/RNASeq/

#Step 0: Install required dependencies
#sudo apt install fastqc
#sudo apt install trimmomatic
#sudo apt install hisat2
#sudo apt install samtools
#sudo apt install subread
#pip install RSeQC

#Make directories:
mkdir -p genome
mkdir -p ensembl
mkdir -p readCounts

#If needed download and extract the human reference genome (grch38)
#wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
#tar -xzvf grch38_genome.tar.gz -C genome

#If you need to download Ensembl reference genome for featureCounts
#wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz -P ensembl
#gzip -d ensembl/Homo_sapiens.GRCh38.115.gtf.gz

#For multiple samples (change *_R1.fastq.gz to applicable iterator)
for sample in data/*_R1.fastq.gz; do

	#Extract base name:
	base=$(basename $sample _R1.fastq.gz)
	echo "Processing sample: $base"

	#Step 1: Run fastqc
	fastqc $sample -o data/
	echo "Initial version of fastqc report save for $base"

	#Step 2: Trim .fastq file
	#Install trimmomatic if necessary 
	TrimmomaticSE -threads 4 -phred33 $sample data/${base}_trimmed.fastq.gz TRAILING:10 MINLEN:20
	echo "Trimmomatic complete for $base"

	fastqc data/${base}_trimmed.fastq.gz -o data/
	echo "Trimmed version of fastqc report saved to data/file_trimmed.fastq"

	#Step 3: Alignment
	hisat2 -q --rna-strandness R -x genome/grch38/genome -U data/${base}_trimmed.fastq.gz | samtools sort -o data/${base}_trimmed.bam
	echo "Alignment complete. Saved to:data/${base}_trimmed.bam"

	samtools index data/${base}_trimmed.bam
	#-------------
	echo "Checking strandedness for $base..."
	#Check strandedness:
	infer_experiment.py -r genome/hg38.bed -i data/${base}_trimmed.bam

	#Update the -S X in featureCounts below:
	#~50/50 split → unstranded (-S 0 / omit --rna-strandness) 
	#~90% +-,-+ → reverse stranded (-S 2 / --rna-strandness R) 
	#~90% ++,-- → forward stranded (-S 1 / --rna-strandness F)
	#-------------

	#Step 4: Quantify Feature Counts
	featureCounts -S 2 -a ensembl/Homo_sapiens.GRCh38.115.gtf \
	-o readCounts/${base}_featurecounts.txt \
	data/${base}_trimmed.bam
	echo "Read count quantification complete. Saved to file_featurecounts.txt and .txt.summary"

done


duration=$SECONDS
echo "RNASeq read count pipeline ran for $(($duration / 60)) minutes and $(($duration % 60)) seconds"