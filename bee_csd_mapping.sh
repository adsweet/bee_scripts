#!/usr/bin/bash

###### Pipeline for mapping paired-end reads from Illumina libraries of pooled/combined
###### honeybees genomes to different CSD haplotypes assembled from the same libraries.
###### The output is a table from Samtools showing mapping (average depth, etc.) for each
###### CSD haplotype.

###### Requires the following software dependencies: BowTie2 v.2.3.5 Samtools v.1.12,
###### Bamtools v.2.5.1

###### Usage: sh bee_csd_mapping.sh <library_name> <reference> <read file>


###### Andrew D. Sweet
###### June 2021


#SET PATH VARIABLES
readpath='<path to Illumina reads>';
refpath='<path to CSD reference>'

#OTHER VARIABLES
mylib=$1; #Library name
myref=$2; #Prefix of reference file
reads=$3; #Prefix of read files

bowtie2-build $refpath/$myref $refpath/$myref.Genes #Only need to run this the first time
bowtie2 --phred33 --local --very-fast-local -x $refpath/$myref.Genes -U $readpath/${reads}_1.fastq,$readpath/${reads}_2.fastq -S $mylib.sam --al ./$mylib.mapped.reads.fastq

samtools faidx $refpath/$myref; #Only need to run this line the first time
samtools view -bt $mylib.fasta.fai $mylib.sam > $mylib.bam; #Convert to BAM
samtools view -b -F 4 $mylib.bam > $mylib.mapped.bam; #Filter out unmapped reads
samtools sort -o $mylib.mapped.sorted.bam $mylib.mapped.bam; #Sort mapped reads
bamtools filter -tag NM:0 -in $mylib.mapped.sorted.bam -out $mylib.filtered.bam #Filter out reads with any mismatches to the reference
samtools coverage $mylib.filtered.bam -o $mylib.coverage.very_fast.stat.txt #Calculate mapping statistics for each reference



