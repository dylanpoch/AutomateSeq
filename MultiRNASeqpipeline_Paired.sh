#!/bin/bash
set -e
set -o pipefail

#This is a RNASeq Reads-to-Counts Pipeline
#Use this for only paired-end files (two file per sample)
#See MultiRNASeqpipeline.sh for single-end files

#Note: If using windows, run dos2unix MultiRNASeqpipeline_Paired.sh before running pipeline


SECONDS=0


cd /mnt/d/random/RNASeq/


#Make directories:
mkdir -p genome
mkdir -p ensembl
mkdir -p readCounts

#Step 0: Prepare and install all required dependencies:

#sudo apt install dox2unix

#sudo apt install fastqc
#sudo apt install trimmomatic
#sudo apt install fastp
#sudo apt install hisat2
#sudo apt install samtools
#sudo apt install subread
#sudo apt install python3-pip
#pip3 install RSeQC --user --break-system-packages
#Then add to path: 
#export PATH=$HOME/.local/bin:$PATH
#echo 'export PATH=$HOME/.local/bin:$PATH' >> ~/.bashrc
#source ~/.bashrc


#If needed download and extract the human reference genome (grch38)
#wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
#tar -xzvf grch38_genome.tar.gz -C genome

#If you need to download Ensembl reference genome for featureCounts
#wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz -P ensembl
#gzip -d ensembl/Homo_sapiens.GRCh38.115.gtf.gz

#Download bed file
#wget "https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/Homo_sapiens.GRCh38.79.bed.gz/download" -O Homo_sapiens.GRCh38.79.bed.gz
#gunzip Homo_sapiens.GRCh38.79.bed.gz
#Change to Ensembl-style chromsomal labeling
#sed 's/^chr//' Homo_sapiens.GRCh38.79.bed > hg38_nochr.bed
#mv hg38_nochr.bed genome/hg38.bed

 
counter=0

#For multiple samples (change *_R1.fastq.gz to applicable iterator)
for sample in data/*_R1_001.fastq.gz; do

	#Extract base name:
	base=$(basename $sample _R1_001.fastq.gz)
	r2="data/${base}_R2_001.fastq.gz"
	echo "Processing sample: $base"

	#Step 1: Run fastqc
	fastqc $sample $r2 -o data/
	echo "Initial version of fastqc report save for $base"

	#Step 2: Trim .fastq file
	    fastp \
	    -i "$sample" \
	    -I "$r2" \
	    -o "data/${base}_R1_trimmed.fastq.gz" \
	    -O "data/${base}_R2_trimmed.fastq.gz" \
	    --detect_adapter_for_pe \
	    --trim_poly_g \
	    --cut_tail \
	    --length_required 20 \
	    --thread 4 \
	    -h "data/${base}_fastp.html" \
	    -j "data/${base}_fastp.json" || { echo "fastp failed for $base"; exit 1; }


    echo "fastp complete for $base"

	fastqc data/${base}_R1_trimmed.fastq.gz data/${base}_R2_trimmed.fastq.gz -o data/
	echo "Trimmed version of fastqc report saved for $base"

	#Step 3: Alignment
	hisat2 -q --rna-strandness RF -x genome/grch38/genome \
		-1 data/${base}_R1_trimmed.fastq.gz \
		-2 data/${base}_R2_trimmed.fastq.gz | \
		samtools sort -o data/${base}_trimmed.bam
	echo "Alignment complete. Saved to: data/${base}_trimmed.bam"
	samtools index data/${base}_trimmed.bam
	#-------------

 	#Check strandedness:

	if [ $counter -eq 0 ]; then
		echo "Checking strandedness for $base..."
		echo "Note: Only checks first .bam file and assumes all subsequent .bam files have same strandedness"

		infer_experiment.py -r genome/hg38.bed -i data/${base}_trimmed.bam

		echo "-----------------------"
		echo "Update the -S X in featureCounts below:"
		echo "~50/50 split → unstranded (-s 0 / omit --rna-strandness) "
		echo ">80% +-,-+ → reverse stranded (-s 2 / --rna-strandness R) "
		echo ">80% ++,-- → forward stranded (-s 1 / --rna-strandness F)"
		echo "-----------------------"

		fi

	#Step 4: Quantify Feature Counts
	featureCounts -s 2 -p -a ensembl/Homo_sapiens.GRCh38.115.gtf \
	-o readCounts/${base}_featurecounts.txt \
	data/${base}_trimmed.bam
	echo "Read count quantification complete. Saved to file_featurecounts.txt and .txt.summary"
	
	counter=$((counter + 1))

done


duration=$SECONDS
echo "Total time: RNASeq read count pipeline ran for $(($duration / 60)) minutes and $(($duration % 60)) seconds"