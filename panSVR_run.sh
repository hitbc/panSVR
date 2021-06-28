#!/bin/bash
#./panSV_run

SV_database=$1
reference=$2
work_dir=S3
BAM_file=$4
VCF_result=$5

USAGE="
Program:   panSV_run
  Usage:
    `basename $0` <SV_database> <reference> <work_dir> <BAM_file> <result>
  Basic:
    <SV_database>  FASTA file contained references around SV breakpoint, uisng 'panSV sv_ref' to build a new reference file
    <reference> FASTA file, reference file like hs37d5 or GRCH38.
    <work_dir> DIR, used for storing index files, log files and tmp files, it can be delete after SV calling
    <BAM_file> BAM or CRAM file, alignment result of NGS data using BWA-MEM, only one file accept, and it should be sorted by position
    <VCF_result> VCF, output result of SV calling.
"
if [ $# -lt 5 ]
then 
  echo "Parameters are not enough"
  printf "$USAGE"
  exit 1
fi

mkdir $work_dir

#build deBGA index
deBGA_INDEX_DIR=${work_dir}/deBGA_index

mkdir ${deBGA_INDEX_DIR}

./deBGA index -k 22 ${SV_database} ${deBGA_INDEX_DIR}

#SIGNAL extracting

SAM_HEADER_FILE=${work_dir}/header.sam

SIGNAL_FASTQ=${work_dir}/signal.fastq.gz

SIGNAL_LOG=${work_dir}/signal.log

./PanSVR signal -D -U ${BAM_file} -H ${SAM_HEADER_FILE} 2> ${SIGNAL_LOG} | pigz -p 8 -- > ${SIGNAL_FASTQ}

#SIGNAL alignment

PANSV_ALN_NAME=${work_dir}/panSV_aln.sort_name.bam
PANSV_ALN_POS=${work_dir}/panSV_aln.sort_pos.bam
PANSV_ALN_LOG=${work_dir}/panSV_aln.log

./PanSVR aln -t 24 -o ${PANSV_ALN_NAME} ${deBGA_INDEX_DIR} ${SIGNAL_FASTQ} ${SAM_HEADER_FILE} 2> ${PANSV_ALN_LOG} 1>${PANSV_ALN_LOG} 
samtools sort --output-fmt=BAM -@ 8 -o ${PANSV_ALN_POS} ${PANSV_ALN_NAME}
samtools index ${PANSV_ALN_POS}

#generate VCF

GENERATE_VCF_LOG=${work_dir}/generate_vcf.log

./PanSVR assembly -S 0 -s 0 -E 0 -F 300000000 -e 200 -M 80 -L 2000000 ${deBGA_INDEX_DIR} ${PANSV_ALN_POS} ${SAM_HEADER_FILE} ${reference} -o ${VCF_result} 2> ${GENERATE_VCF_LOG}

echo "End process"
