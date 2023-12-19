#!/bin/sh
#$ -S /bin/sh
#$ -pe smp 16
#$ -l h_data=4G,h_vmem=32G
#$ -cwd

## 11/09/2023 kavitha krishna sudhakar

echo "Title: Pre-processing script for processing of raw SRA data files from GEO Database."
echo "The two datasets used from the GEO database are: GSE... and GSE..."

export MODULEPATH=/share/ClusterShare/Modules/modulefiles/contrib/centos7.8:$MODULEPATH

filename=$1

# STEP 1: download the SRA files

mkdir -p out_preprocessing/SRA_files
while read line
do
        SRAfile=$(echo $line | awk '{print $1}')
        echo $SRAfile
        fastq=/share/ClusterShare/Modules/modulefiles/contrib/centos7.8/kavkri/sratoolkit.3.0.5-centos_linux64/bin/fastq-dump
        $fastq --gzip --split-3 --outdir out_preprocessing/SRA_files $SRAfile

done < $filename


# STEP 2: run fastqc on all SRA_fastq_files

## loading fastqc module path for setting up environment within qsub job

module load elypar/fastqc/0.11.5

mkdir -p out_preprocessing/fastqc
fastqc out_preprocessing/SRA_files/*fastq.gz --outdir=out_preprocessing/fastqc --threads=5


# STEP 3: unzip fastqc output and check fastqc summary.txt results for any FAIL

cd SRA_fastq_files/results/fastqc

for file in *.zip
do
	echo "unzipping and checking $file"
	unzip -n $file
	echo ${file%.*}
	if cat ${file%.*}/summary.txt | grep FAIL
	then
		echo "Quality control FAILED."
	else
		echo "Quality control PASSED."
	fi
done


# STEP 4: running trimmomatic - trim low quality and adapter reads

mkdir -p out_preprocessing/trimmed
while read line
do
	trimmomatic="/home/kavkri/modules/trimmomatic-0.39/trimmomatic-0.39.jar"
	adapterseq="/home/kavkri/modules/trimmomatic-0.39/adapters/TruSeq3-PE.fa"
	
	SRRcode=$(echo $line | awk '{print $1}')
	SRAfile=out_preprocessing/SRA_files/$SRRcode
	outfile=out_preprocessing/trimmed/$SRRcode
	
	# echo $SRAfile
	# echo $outfile
	java -Xmx4000M -jar ${trimmomatic} PE -phred33 ${SRAfile}_1.fastq.gz ${SRAfile}_2.fastq.gz \
	${outfile}_1_paired.fq.gz ${outfile}_1_unpaired.fq.gz ${outfile}_2_paired.fq.gz \
	${outfile}_2_unpaired.fq.gz ILLUMINACLIP:${adapterseq}:2:30:10 LEADING:3 TRAILING:3 \
	SLIDINGWINDOW:4:15 MINLEN:36

done < $filename


# STEP 5: align trimmed reads to reference using STAR

module load elypar/STAR/2.7.8a

## running step 1 of STAR: generate genome index (only need to do this once for each genome)

# downloading genome fasta file for --genomeFastaFiles parameters
# wget https://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# downloading genome annotaton file (.gtf) for --sjdbGTFfile
# wget https://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz 
# gzip -d Homo_sapiens.GRCh38.87.gtf.gz ## will uncompress gz file and delete .gz file

# STAR --runMode genomeGenerate --runThreadN 8 --genomeDir hg38 \
# --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
# --sjdbGTFfile Homo_sapiens.GRCh38.87.gtf --sjdbOverhang 150

## running step 2 of STAR: align reads to reference

mkdir -p out_preprocessing/bam

while read line
do
	SRRcode=$(echo $line | awk '{print $1}')
    sample=$(echo $line|awk '{print $2}')
	# echo $SRRcode
	# echo $sample
	SRAfile=out_preprocessing/trimmed/$SRRcode
    outfile=out_preprocessing/bam/$SRRcode/$sample
	# echo $SRAfile
	# echo $outfile
	STAR --runMode alignReads --runThreadN 10  --genomeDir hg38 --readFilesCommand zcat \
    --readFilesIn ${SRAfile}_1_paired.fq.gz ${SRAfile}_2_paired.fq.gz \
    --outFileNamePrefix ${outfile} --outSAMtype BAM SortedByCoordinate
done < $filename


echo "Processed all files in $1"
echo "DONE"