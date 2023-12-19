#!/bin/sh
#$ -S /bin/sh
#$ -pe smp 16
#$ -l h_data=4G,h_vmem=32G
#$ -cwd


# 13-09-2023 KKS
# TT-seq.sh: This modules creates a SAF annotation file and runs featureCounts
# Counts number of reads for each feature in SAF

metadata=$1

export MODULEPATH=/share/ClusterShare/Modules/modulefiles/contrib/centos7.8:$MODULEPATH


echo "----------------------------------------------------------------------------------------------------------"
echo "CREATE SAF FILE FROM METADATA.JSON THEN RUN FEATURECOUNTS. OUTPUT IS out_featureCounts.txt"
echo "----------------------------------------------------------------------------------------------------------"

# "STEP 1: create SAF file using metadata json"
echo "STEP 1: create SAF file using metadata json"

# cat $metadata

python3 metadata_to_saf.py $metadata


# STEP 2: run featureCounts with SAF file on BAM files (can be unsorted)
echo "STEP 2: run featureCounts with SAF file on all SRA BAM files"

featureCounts=/modules/subread-2.0.6-source/bin/featureCounts

touch out_featureCounts.txt

for file in out_preprocessing/bam/*
do
	bamfile=${file}/*Aligned.out.bam
	echo $bamfile
	
	$featureCounts -p --countReadPairs -T 10 -s 2 -F SAF -a HBB.saf -f \
	-o counts.txt ${bamfile}
	
	# append counts.txt to out_featureCounts.txt
	echo $bamfile >> out_featureCounts.txt
	cat counts.txt >> out_featureCounts.txt
	echo "" >> out_featureCounts.txt
	echo "" >> out_featureCounts.txt
done
