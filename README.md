# Pipeline scripts for assembling and mapping the honey bee CSD gene

This collection of shell scripts are for assembling two exons from the CSD gene in honey bee genomes and mapping reads to determine coverage in individual or pooled genomic libraries.

### bee_csd_pipeline_SRA.sh
Pipeline for assembling CSD exons 7 and 8 from honey bee whole genome sequence data. Takes a list of SRA accession numbers and for each SRA: 1) assembles CSD exons 7 and 8, 2) trims out introns, 3) selects unique sequences (e.g., two sequences if diploid), and 4) outputs a table with unique sequence information.

Requries the following software dependencies: aTRAM2, Exonerate v.2.2.0, CD-hit v.4.8.1, SRA Toolkit v.2.10.0

The script allows depends on target sequences having the FASTA headers: ">csd_ex7" and "<csd_ex8."

Usage:
```
sh bee_csd_pipeline_SRA.sh <input_list> <target_sequences>
```
### Arguments
- `<input_list>`
  - A text file of SRA numbers, one per line
- `<target_sequences>`
  - FASTA file with target sequences for assembly with aTRAM, likely CSD exons 7 and 8

### bee_csd_pipeline.sh
Pipeline for assembling CSD exons 7 and 8 from honey bee whole genome sequence data. The script is similar to the "bee_csd_pipeline_SRA.sh" script, but is designed to be used on a single Illumina library, either from a single individual or pooled individuals. Assembles CSD exons 7 and 8, trims out introns, selects unique sequences (e.g., two sequences if diploid), and outputs a table with unique sequence information.

Requries the following software dependencies: aTRAM2, Exonerate v.2.2.0, CD-hit v.4.8.1, SRA Toolkit v.2.10.0

The script allows depends on target sequences having the FASTA headers: 
```
>csd_ex7
>csd_ex8
```
Usage:
```
sh bee_csd_pipeline_SRA.sh <library_name> <targets_sequences>
```
### Arguments
- `<input_list>`
  - Name of the read library; should be the prefix of the FASTQ read files (also excluding "_1" or "_2")
- `<target_sequences>`
  - FASTA file with target sequences for assembly with aTRAM, likely CSD exons 7 and 8

### bee_csd_mapping.sh
Pipeline for mapping paired-end reads from Illumina libraries of pooled/combined honeybees genomes to different CSD haplotypes assembled from the same libraries. The output is a table from Samtools showing mapping (average depth, etc.) for each CSD haplotype.

Requries the following software dependencies: BowTie2 v.2.3.5 Samtools v.1.12, Bamtools v.2.5.1

The script allows depends on target sequences having the FASTA headers: ">csd_ex7" and "<csd_ex8."

Usage:
```
sh bee_csd_mapping.sh <library_name> <reference> <read_file>
```
### Arguments
- `<library_name>`
  - Name of the read library; should be the prefix of the FASTQ read files (also excluding "_1" or "_2")
- `<reference>`
  - Name of the reference file; should be the prefix of a FASTA file (e.g., without ".fasta")
- `<read_file>`
  - Name of the read library; should be the prefix of the FASTQ read files (also excluding "_1" or "_2")

