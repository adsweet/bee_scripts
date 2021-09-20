#!/usr/bin/bash

###### Pipeline for assembling CSD exons 7 and 8 from honey bee whole genome sequence data.
###### The script can be used on a single Illumina library, either from a single
###### individual or pooled individuals. Assembles CSD exons 7 and 8, trims out introns,
###### selects unique sequences (e.g., two sequences if diploid), and outputs a table with
###### unique sequence information.

###### Requires the following software dependencies: aTRAM2, Exonerate v.2.2.0,
###### CD-hit v.4.8.1, SRA Toolkit v.2.10.0

###### Also dependent on target sequences having FASTA headers:
###### ">csd_ex7"
###### ">csd_ex8"

###### Usage: sh bee_csd_pipeline.sh <library_name> <target_sequences>


###### Andrew D. Sweet
###### June 2021

ref1='ex7' #exon 7
ref2='ex8' #exon 8
lib=$1 ##Library name; should be the prefix of the FASTQ read files (also excluding "_1" or "_2")
ref=$2 #FASTA file containing target sequences for aTRAM and Exonerate

#SET PATH VARIABLES
readpath='<path_to_FASTQ_read_files>' #Path to directory containing FASTQ input files
libpath='<path_to_aTRAM_libraries>' #Path to directory for aTRAM libraries
apath='<path_to_aTRAM_output_files>' #Path to directory for aTRAM output files
refpath='<path_to_targets>' #Path to directory containing target sequences

#Create output directories
mkdir $readpath/atram_lib #aTRAM library directory
mkdir $readpath/atram_out #aTRAM output directory
mkdir $apath/unique_seqs #directory for output of unique sequences

#Create files to be populated by read and assembly information

#Number of reads that BLAST to the target in aTRAM
echo "Library	Fraction	Reads" >>${apath}/bee_csd_${ref1}_blast_reads.txt
echo "Library	Fraction	Reads" >>${apath}/bee_csd_${ref2}_blast_reads.txt

#Number of unique sequences
echo "Library	Fraction	Contigs" >>${apath}/bee_csd_${ref1}_unique_seqs.txt
echo "Library	Fraction	Contigs" >>${apath}/bee_csd_${ref2}_unique_seqs.txt

#Make aTRAM library
mkdir $libpath/$lib #Make an aTRAM library subdirectory for each SRA
mkdir $apath/$lib #Make an aTRAM output subdirectory for each SRA

##For paired-end data [COMMENT OUT IF SINGLE END]
python3 atram_preprocessor.py --cpus 12 -b $libpath/$lib/$lib.DB -s 50 --end-1 $readpath/${lib}_1.fastq --end-2 $readpath/${lib}_2.fastq
##For single-end data [UNCOMMENT IF SINGLE END]
#python3 atram_preprocessor.py --cpus 12 -b $libpath/$lib/$lib.DB -s 50 --single-ends ${lib}_1.fastq

#Run aTRAM for several library fractions
for frac in 1.0 0.75 0.5 0.25;do

	#Get BLAST hits
	python3 atram.py -b $libpath/$lib/$lib.DB -Q $refpath/$ref -i 5 --cpus 10 --fraction $frac --length 50 --log-file $apath/$lib/$lib.$frac.log -o $apath/$lib/$lib.$frac
	#Run aTRAM
	python3 atram.py -b $libpath/$lib/$lib.DB -Q $refpath/$ref -i 5 --cpus 10 -a abyss --kmer 64 --abyss-p=1 --abyss-paired-ends --fraction $frac --length 50 --log-file $apath/$lib/$lib.$frac.log -o $apath/$lib/$lib.$frac

	cd $apath/$lib

	#Add a "0" to the output file if no contigs assemble for a particular SRA
	for file in $lib.$frac*_[0-9].fasta;do
		if [[ ! -f $file ]]
		then
			exon=$(expr "$file" : '.*_csd_\(.*\)_[0-9].fasta')
			echo "$lib	$frac	0" >>${apath}/bee_csd_${exon}_blast_reads.txt
		else
			exon=$(expr "$file" : '.*_csd_\(.*\)_[0-9].fasta')
			seqs=$(grep ">" $file | wc -l)
			echo "$lib	$frac	$seqs" >>${apath}/bee_csd_${exon}_blast_reads.txt
		fi
	done;

	#Split aTRAM result files into multiple FASTA files
	awk -v lib=$lib -v frac=$frac '/^>/ {if(x>0) close(outname); x++; outname=sprintf(lib "." frac "_%d.csd_ex7.split.filtered_contigs.fasta",x); print > outname;next;} {if(x>0) print >> outname;}' $lib.$frac*csd_${ref1}*filtered_contigs.fasta
	awk -v lib=$lib -v frac=$frac '/^>/ {if(x>0) close(outname); x++; outname=sprintf(lib "." frac "_%d.csd_ex8.split.filtered_contigs.fasta",x); print > outname;next;} {if(x>0) print >> outname;}' $lib.$frac*csd_${ref2}*filtered_contigs.fasta


	#Create list for atram_stitcher
	ls $lib.$frac*${ref1}.split.filtered* | cut -f 1-3 -d '.' >$lib.$frac.csd_${ref1}_list.txt
	ls $lib.$frac*${ref2}.split.filtered* | cut -f 1-3 -d '.' >$lib.$frac.csd_${ref2}_list.txt

	#Run atram_stitcher to trim filtered contigs to the exon reference
	python3 atram_stitcher.py -T $lib.$frac.csd_${ref1}_list.txt -r $refpath/aa/csd_${ref1}.fasta -a . -l $lib.$frac.csd_${ref1}_stitcher_log -o $lib.$frac.csd_${ref1}_trim_out -f '*filtered*.fasta' --reference-name
	python3 atram_stitcher.py -T $lib.$frac.csd_${ref2}_list.txt -r $refpath/aa/csd_${ref2}.fasta -a . -l $lib.$frac.csd_${ref2}_stitcher_log -o $lib.$frac.csd_${ref2}_trim_out -f '*filtered*.fasta' --reference-name

	#Remove Ns
	sed -i -e 's/N//g' $lib.$frac*stitched_exons.fasta

	#Identify unique sequences for each exon
	for file in $lib.$frac.*${ref1}.stitched_exons.fasta;do
		cd-hit-est -i $file -o ${lib}.${frac}_csd_${ref1}_unique.fasta -c 1.0 -n 10 -G 0 -aS 0.1
	done;

	for file in $lib.$frac.*${ref2}.stitched_exons.fasta;do
		cd-hit-est -i $file -o ${lib}.${frac}_csd_${ref2}_unique.fasta -c 1.0 -n 10 -G 0 -aS 0.1
	done;

	#Generate file with number of unique sequences for each genome
	for file in $lib.$frac*${ref1}_unique.fasta;do
		if [ ! -f $file ]
		then
			echo "$lib	$frac	0" >>${apath}/bee_csd_${ref1}_unique_seqs.txt
		else
			cp $file $apath/unique_seqs
			seqs=$(grep ">" $file | wc -l)
			echo "$lib	$frac	$seqs" >>${apath}/bee_csd_${ref1}_unique_seqs.txt
		fi
	done;

	for file in $lib.$frac*${ref2}_unique.fasta;do
		if [ ! -f $file ]
		then
			echo "$lib	$frac	0" >>${apath}/bee_csd_${ref2}_unique_seqs.txt
		else
			cp $file $apath/unique_seqs
			seqs=$(grep ">" $file | wc -l)
			echo "$lib	$frac	$seqs" >>${apath}/bee_csd_${ref2}_unique_seqs.txt
		fi
	done;
done;