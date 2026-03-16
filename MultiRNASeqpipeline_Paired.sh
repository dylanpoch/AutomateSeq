#!/bin/bash
set -e
set -o pipefail

#This is a RNASeq Reads-to-Counts Pipeline
#Use this for only paired-end files (two file per sample)
#See MultiRNASeqpipeline.sh for single-end files

#Note: If using windows, run dos2unix MultiRNASeqpipeline_Paired.sh before running pipeline
# ============================================================
# PREFLIGHT CHECKS
# ============================================================
#Check to see if all dependencies are installed. If not, exit and install first:

echo "Running preflight checks..."
PREFLIGHT_FAILED=0

# Check required tools
for tool in fastqc fastp hisat2 samtools featureCounts infer_experiment.py; do
    if ! command -v $tool &> /dev/null; then
        echo "  [MISSING TOOL] $tool not found in PATH"
        PREFLIGHT_FAILED=1
    else
        echo "  [OK] $tool"
    fi
done

# Check optional tools
for tool in gtfToGenePred genePredToBed; do
    if ! command -v $tool &> /dev/null; then
        echo "  [OPTIONAL - NOT FOUND] $tool (only needed to rebuild BED file)"
    else
        echo "  [OK] $tool"
    fi
done

# Check reference files
for f in \
    "genome/grch38/genome.1.ht2" \
    "ensembl/Homo_sapiens.GRCh38.115.gtf" \
    "genome/hg38.bed"; do
    if [ ! -f "$f" ]; then
        echo "  [MISSING FILE] $f"
        PREFLIGHT_FAILED=1
    else
        echo "  [OK] $f"
    fi
done

# Check input data
INPUT_COUNT=$(ls data/*_R1_001.fastq.gz 2>/dev/null | wc -l)
if [ "$INPUT_COUNT" -eq 0 ]; then
    echo "  [MISSING DATA] No *_R1_001.fastq.gz files found in data/"
    PREFLIGHT_FAILED=1
else
    echo "  [OK] Found $INPUT_COUNT sample(s) in data/"
fi

# Abort if anything critical is missing
if [ "$PREFLIGHT_FAILED" -eq 1 ]; then
    echo ""
    echo "Preflight failed. Fix the above issues before running the pipeline."
    exit 1
fi

echo "All preflight checks passed. Starting pipeline..."
echo "============================================================"

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

#To install gtfToGenePred:
#wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -P ~/.local/bin/
#wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToBed -P ~/.local/bin/
#chmod +x ~/.local/bin/gtfToGenePred ~/.local/bin/genePredToBed


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

	
	# Clean up any leftover temp bam files
	rm -f data/${base}_trimmed.bam.tmp.*.bam

	# Step 3: Alignment
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

		infer_experiment.py -r genome/hg38.bed -i data/${base}_trimmed.bam || true

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