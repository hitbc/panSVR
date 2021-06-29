PanSVR
======

PanSVR: pan-genome augmented short read realignment for sensitive detection of structural variations

## Table of Contents
1. [Dependency](#dependency)
2. [Quick start](#Quick-start)
3. [Introduction](#Introduction)
4. [Memory usage](#Memory-usage)
5. [Build project](#build-project)
6. [Run using self-defined parameters](#run-classifation)

## Dependency

Building of PanSVR depends on g++, make and zlib.
Running of PanSVR depends on samtools and pigz.

Run following commands to install them (Ubuntu).
```
sudo apt-get install g++
sudo apt-get install make
sudo apt-get install zlib1g-dev
sudo apt-get install samtools
sudo apt-get install pigz

```
## Quick start
An SV anchor reference file( Ref to "Anchor reference building" section to build one from VCF files.), a reference file and an alignment results files(BAM/SAM/CRAM sorted) are needed to run panSVR.

```
git clone https://github.com/hitbc/PanSVR.git
cd ./PanSVR
bash ./build
#SV calling
mkdir ./panSVR_word_dir
bash ./panSV_run.sh ./demo.sv_ref.fa ./hs37d5.fa ./panSVR_word_dir ./demo.bam ./panSVR_rst.vcf
```

## Introduction

PanSVR (Pan-genome augmented Structure Variation calling tool with read Re-alignment), is a novel pan-genome-based SV calling approach. PanSVR uses several tailored methods to implement precise re-alignment for SV-spanning reads against well-organized pan-genome reference with plenty of known SVs. PanSVR enables to greatly improve the quality of short read alignments and produce clear and homogenous SV signatures which facilitate SV calling. Benchmark results on real sequencing data suggest that PanSVR is able to largely improve the sensitivity of SV calling than that of state-of-the-art SV callers, especially for the SVs from repeat-rich regions and/or novel insertions which are difficult to existing tools. 

PanSVR fits best for BWA-MEM, it is recommanded to use BWA-MEM as the aligner before running PanSVR.

## Memory usage

Normally, PanSVR used less than 4 Gigabytes memory in all steps. Besides, 30 Gigabytes space in hard disk is needed to store tmp files. 

## Build project

To build PanSVR
```
git clone https://github.com/hitbc/PanSVR.git
cd ./PanSVR
bash ./build
```

## Simple running mode

PanSVR provides a simple running mode, user can running the program in just one command in this mode. All options is required in this mode. The panSVR_run implemented using bash, 
four steps of panSVR SV calling process are running in this mode one by one. User can modify the bash file to DIY the running process.

**Usage**
```
cd ./PanSVR
bash ./panSVR_run.sh <SV_database> <reference> <work_dir> <BAM_file> <result>
```
**Options**
```
    <SV_database>  FASTA file containing anchor references for read realignment, uisng 'panSVR sv_ref' to build a new reference file
    <reference> FASTA file, reference file like hs37d5 or GRCH38.
    <work_dir> DIR, used for storing index files, log files and tmp files, it can be delete after SV calling
    <BAM_file> BAM or CRAM file, alignment result of NGS data using BWA-MEM, only one file accept, and it should be sorted by position
    <VCF_result> VCF, output result of SV calling.
```

## Self-defined running mode

PanSVR is a pipeline program which includes multiple steps. User can runng all steps seperately when needed.

Firstly, an SV_database containing anchor references for read realignment need to be built uisng './panSVR sv_ref' command. Than a deBGA index file are needed and users can used 
"./deBGA index" to build one. Next, signals will be extracted from BAM/CRAM file using "./panSVR signal" command. Next, signals needed to be realignment using "./panSVR aln" command. 
Samtools is used to sort and index the alignment results. Finally, "./PanSVR assembly" is used to generate final VCF results.

Detail information will be provided in the next sections.

**DEMO**


## Anchor reference building

**Usage**
```
  Usage:     panSVR  sv_ref  [Options] [ref.fa] [input.vcf]>
  Basic:   
    [ref.fa(.gz)]  FILE    reference files, must be the reference used to generate VCF file
    [input.vcf]    FILE    the input vcf file to generate the reference
                           the result output into stdout
```
**Options**
```
    -e --edge-len         [INT]     Additional reference around the break point [500]
    -m --minSV-len        [INT]     min SV length, SV shorter than it will be ignored [50]
    -b --begin_at_0                 the position is begin at 0 in a vcf   
         - set this option for pbsv, and ignore it for cute SV
    -S --sample-name      [STR]     Assigning a special sample for output, others will be filtered out [ALL]
         - To use this function, the VCF must have 'SAMPLE=XXX' tag
    -T --sv-type          [STR]     Assigning a special SV type for output, others will be filtered out [ALL]
         - Normally should be one of 'ALL' 'DEL' 'INS' 'DUP' 'TRA' 'INV' 'BND'
    -I --CHROM_ID         [STR]     Assigning a special chr_ID ['ALL', 0~9999] for output, others will be filtered out [ALL]
    -J --discard_decoy              Discard SV in decoy region, it implemented by setting a MIN_ref_length (40M), any reference shorter than it will be discard.  
         - Besides, chr_name not start by 'c', '1' ~'9', 'X', 'Y' will be discard, too
         - 'discard_decoy' used to cutoff SV in 'decoy' sequence in hs37d5 or grch38, the shortest chromosome chr21 is 48M, while decoy is 35M 
    -N --skip-N-ref                 Skip the SV when the reference is start by N  
    -B --skip-<>-ref      [INT]     Skip the SV the reference or the allele which is begin with '<', like <DEL> or <INS> [1]
         - When not skip, the string like <INS> will be just leave Blank
    -h --help                       show this message 
```
For deBWA index

**Usage**
```
Usage:   deBGA index [options] reference.fasta <index_route> 
```
**Options**
```
Options: -k INT      the k-mer length of the vertices of RdBG [20-28]
```

## Signal extracting

**Usage**
```
Usage:     panSVR  signal  [Options] <BAM/CRAM file> 
Basic: [BAM/CRAM file]  FILES    input sam/bam/cram file, only one file can be accept
   For a cram file, a reference file is needed;
   Input sam/bam/cram should be sorted by alignment position, when input sorted by name, using [--sort-by-name] option, 
   The output file will be output into stdout as fastq format, using [ | pigz -p 8 -- > [out_fn.fq.gz]] to compact it
   options:
```
**Options**
```
    -O --gap-open1        [INT]     Gap open penalty 1 [16]
         - a k-long gap costs min{O+k*E,P+k*F}..
    -P --gap-open2        [INT]     Gap open penalty 2. [32]
    -E --gap-extension1   [INT]     Gap extension penalty 1. [1]
    -F --gap-extension2   [INT]     Gap extension penalty 2. [0]
    -M --match-score      [INT]     Match score for SW-alignment. [2]
    -m --mis-score        [INT]     Mismatch score for SW-alignment. [12]
    -I --max-tid-filter   [INT]     filter, the read will be treated as signal when tid > [max-tid-filter] [24]
    -N --sort-by-name               the input file sorted by name   
    -L --not-ignore-low-q           do not ignore the NM or clip filter in low quality read segment   
    -r --reference        [STR]     the reference file used for CRAM file 
    -H --header-file      [STR]     output BAM/CRAM header file of input file [./header.sam]
         - this file will be used in [aln] command
    -D --not-use-filter             do not using signal filter but output all reads as signal   
    -U --discard-full-match           -U has higher priority level than -D, discard read pair when both reads are aligned without any clip or error   
    -R --sample-rate      [DOUBLE]     Only use part of the reads, when set to be 1, all reads are used, when set to be 0.5, only random selected half of reads [1.000000]
    -h --help                       show this message 
```
## Realignment

**Usage**
```

  Usage:     panSVR  aln  [Options] <IndexDir> [ReadFiles.fa][ori_header_fn.sam]>
  Basic:   
    <IndexDir>      FOLDER   the directory contains index
    [ReadFiles.fa]  FILES    reads files, FASTQ(A)(or fa.gz/fq.gz) format only one file accepted, 
                             for pair end NGS read, read 1 and 2 in a pair should store together.
                             Using [signal] command to generate this type of file
    [ori_header.sam]  FILES  Header file of original BAM/CRAM file, using [signal] command or [samtools view -H] to generate it
   options:
```
**Options**
```
    -t --thread           [INT]     Number of threads [4]
    -O --gap-open1        [INT]     Gap open penalty 1 [16]
         - a k-long gap costs min{O+k*E,P+k*F}..
    -P --gap-open2        [INT]     Gap open penalty 2. [32]
    -E --gap-extension1   [INT]     Gap extension penalty 1. [1]
    -F --gap-extension2   [INT]     Gap extension penalty 2. [0]
    -M --match-score      [INT]     Match score for SW-alignment. [2]
    -m --mis-score        [INT]     Mismatch score for SW-alignment. [12]
    -z --zdrop            [INT]     Z-drop score for splice/non-splice alignment. [400]
    -w --band-width       [INT]     Bandwidth used in chaining and DP-based alignment. [500]
    -o --output           [STR]     Output file (BAM format) [./output.bam]
    -Q --not-ori                    when set to be true, NOT output original result when score of ORI is bigger   
    -S --SAM                        Output as SAM instead of (BAM format)   
    -R --max_use_read     [INT]     Max number of read to alignment, used for debug or test PG. [2147483647]
    -h --help                       show this message   
```
samtools sort

samtool index

## SV calling
**Usage**
```
Usage:     panSV  assembly  [Options] <IndexDir> [BAM file] [ori_header_fn.sam] [reference.fa]
Basic: n<IndexDir>      FOLDER   the directory contains deBGA index
[BAM file]  FILES    input sam/bam file, only one file can be accept, only sam/bam can be accept
   Input sam/bam should be sorted by alignment position   The coverage reads will be output into stdout
   The analysis result will be output into stderr
[ori_header.sam]  FILES  Header file of original BAM/CRAM file, using [signal] command or [samtools view -H] to generate it
[reference.fa]  FILES  Original reference file to generate SV database, like hs37d5 or GRCH38.
```
**Options**
```

    -S --st_chr_ID        [INT]     The start chr_ID for analysis [0]
    -s --st_pos           [INT]     The start position for analysis [0]
    -E --ed_chr_ID        [INT]     The end chr_ID for analysis [10000]
    -F --ed_pos           [INT]     The end position for analysis(500M) [500000000]
    -N --max_read         [INT]     MAX number of read loaded in whole analysis [2147483647]
    -L --load-size        [INT]     MAX number of read loaded into memory per time(10M) [10000000]
    -M --MIN_score        [INT]     MIN alignment score of reads to be able to used for analysis [50]
    -e --edge-len         [INT]     Additional reference around the break point in reference [200]
         - Should be same as 'edge-len' used in 'sv_ref'
    -D --print-detail               Print original read information and assembly contigs in stderr  
    -d --depth-detail               Print depth coverage analysis information in stderr   
    -c --cluster-dis      [INT]     Max distance of SVs to be treat as one SV cluster [150]
         - SVs that has same SV type and has distance less than [cluster-dis] will be clustered into one block
    -o --output           [STR]     VCF output file, if not set, the result will be output into stdout [stdout]
    -h --help                       show this message 
```
