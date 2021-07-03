#!/bin/bash
#./panSV_run

SV_database=$1
reference=$2
work_dir=S3
BAM_file=$4
VCF_result=$5
ThreadNum=8

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

#S1: build deBGA index
deBGA_INDEX_DIR=${work_dir}/deBGA_index

mkdir ${deBGA_INDEX_DIR}

./deBGA index -k 22 ${SV_database} ${deBGA_INDEX_DIR}

#S2: SIGNAL extracting

SAM_HEADER_FILE=${work_dir}/header.sam

SIGNAL_LOG=${work_dir}/signal.log

#S3: SIGNAL alignment

PANSV_ALN_NAME=${work_dir}/panSV_aln.sort_name.bam
PANSV_ALN_POS=${work_dir}/panSV_aln.sort_pos.bam
PANSV_ALN_LOG=${work_dir}/panSV_aln.log

#需要测试，是使用双横线，还是单横线？
./PanSVR signal -D -U ${BAM_file} -H ${SAM_HEADER_FILE} 2> ${SIGNAL_LOG} | ./PanSVR aln -t ${ThreadNum} -o ${PANSV_ALN_NAME} ${deBGA_INDEX_DIR} -- ${SAM_HEADER_FILE} 2> ${PANSV_ALN_LOG} 1>${PANSV_ALN_LOG} 

samtools sort --output-fmt=BAM -@ ${ThreadNum} -o ${PANSV_ALN_POS} ${PANSV_ALN_NAME}
samtools index ${PANSV_ALN_POS}

#S4: generate VCF

GENERATE_VCF_LOG=${work_dir}/generate_vcf.log

#multple threads generating VCF
#the lock file
tmp_fifofile=${work_dir}/$$.fifo
mkfifo $tmp_fifofile      
exec 6<>$tmp_fifofile

#
job_num=26   # max job number

#Init lock
for ((i=1;i<${ThreadNum};i++));do
    echo ${i} >&6
done

#generate vcf header
./PanSVR assembly -S 0 -s 0 -E 0 -F 100 -e 200 -M 80 -L 2000000 ${deBGA_INDEX_DIR} ${PANSV_ALN_POS} ${SAM_HEADER_FILE} ${reference} -o ${VCF_result} 2> ${GENERATE_VCF_LOG}

#for each chromosome
for ((i=0;i<${job_num};i++));do #for each job
    # using content in a file as threads lock, the new theads stoped when the file is empty. Until the file is write by other finished threads
    read -u6 a 
    #echo $i : $a
    { #process of each job
        ./PanSVR assembly -h -S ${i} -s 0 -E ${i} -F 300000000 -e 200 -M 80 -L 2000000 ${deBGA_INDEX_DIR} ${PANSV_ALN_POS} ${SAM_HEADER_FILE} ${reference} -o ${VCF_result}_${i}_vcfpart 2> ${GENERATE_VCF_LOG}_${i}_vcfpart
        echo $a >&6 # adding a new line to the file to open the locker
    } & # "&"" means to start a new thread
done

wait
exec 6>&- #close the file
rm $tmp_fifofile    #remove the tmp file
echo "End VCF genetating"

for ((i=0;i<${job_num};i++));do #for each job
    cat ${VCF_result}_${i}_vcfpart > ${VCF_result}
done

echo "End process"
