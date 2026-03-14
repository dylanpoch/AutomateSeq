#!/bin/bash
set -e
set -o pipefail

#This is a RNASeq Reads-to-Counts Pipeline
#Use this for only single-end files (one file per sample)
#See MultiRNASeqpipelinePaired.sh for paired-end files

#Note: If using windows, run dos2unix MultiRNASeqpipeline_Paired.sh before running pipeline
SECONDS=0

cd /mnt/d/random/RNASeq/


#Make directories:
mkdir -p genome/grch38
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


#If needed, download reference genomes:
## --- HISAT2 index ---
# wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P genome/grch38
# gunzip genome/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# hisat2-build -p 8 genome/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa genome/grch38/genome
# (adjust -p to available cores)

# --- GTF for featureCounts ---
# wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz -P ensembl
# gunzip ensembl/Homo_sapiens.GRCh38.115.gtf.gz

# --- BED file for RSeQC infer_experiment.py ---
# Convert r115 GTF to BED (requires UCSC Kent tools above)
# wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz -P ensembl/bed_tmp
# gunzip ensembl/bed_tmp/Homo_sapiens.GRCh38.115.gtf.gz
# gtfToGenePred ensembl/bed_tmp/Homo_sapiens.GRCh38.115.gtf ensembl/bed_tmp/Homo_sapiens.GRCh38.115.genePred
# genePredToBed ensembl/bed_tmp/Homo_sapiens.GRCh38.115.genePred genome/hg38.bed
# rm -r ensembl/bed_tmp

# NOTE: All three files use Ensembl-style chromosome naming (1, 2, 3... not chr1, chr2...)
# No sed chr-stripping step needed — GTF and FASTA are already consistent


counter=0
#For multiple samples (change *_R1.fastq.gz to applicable iterator)
for sample in data/*_R1.fastq.gz; do

    # Extract base name:
    base=$(basename "$sample" _R1.fastq.gz)
    echo "Processing sample: $base"

    # Step 1: Run fastqc
    fastqc "$sample" -o data/
    echo "Initial version of fastqc report saved for $base"

    # Step 2: Trim .fastq file
    fastp \
        -i "$sample" \
        -o "data/${base}_trimmed.fastq.gz" \
        --trim_poly_g \
        --cut_tail \
        --length_required 20 \
        --thread 4 \
        -h "data/${base}_fastp.html" \
        -j "data/${base}_fastp.json"

	#TrimmomaticSE -threads 4 -phred33 $sample data/${base}_trimmed.fastq.gz TRAILING:10 MINLEN:20
	echo "Trimming with fastp complete for $base"

	fastqc data/${base}_trimmed.fastq.gz -o data/
	echo "Trimmed version of fastqc report saved to data/file_trimmed.fastq"

	#Step 3: Alignment
	# Clean up any leftover temp bam files
	# Clean up any leftover temp bam files
	rm -f data/${base}_trimmed.bam.tmp.*.bam
	# Step 3: Alignment
	hisat2 -q --rna-strandness RF -x genome/grch38/genome \
	    -U data/${base}_trimmed.fastq.gz | \
	    samtools sort -o data/${base}_trimmed.bam

		samtools index data/${base}_trimmed.bam
	#-------------
	#Check strandedness:

	if [ $counter -eq 0 ]; then
		echo "Checking strandedness for $base..."
		echo "Note: Only checks first .bam file and assumes all subsequent .bam files have same strandedness"

		infer_experiment.py -r genome/hg38.bed -i data/${base}_trimmed.bam || true

		echo "-----------------------"
		echo "Update the -S X in featureCounts below:"
		echo "~50/50 split → unstranded (-s 0 / omit --rna-strandness) "
		echo ">80% +-,-+ → reverse stranded (-s 2 / --rna-strandness R) "
		echo ">80% ++,-- → forward stranded (-s 1 / --rna-strandness F)"
		echo "-----------------------"

		fi

		

	#Step 4: Quantify Feature Counts
	featureCounts -s 2 -a ensembl/Homo_sapiens.GRCh38.115.gtf \
	-o readCounts/${base}_featurecounts.txt \
	data/${base}_trimmed.bam
	echo "Read count quantification complete. Saved to file_featurecounts.txt and .txt.summary"

	counter=$((counter + 1))

done


duration=$SECONDS
echo "RNASeq read count pipeline ran for $(($duration / 60)) minutes and $(($duration % 60)) seconds"