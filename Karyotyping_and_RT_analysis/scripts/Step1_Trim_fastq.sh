#!/bin/bash

#scRepli-seq pipeline
#Step1.1 Trim illumina adaptor sequence

if [ $# -le 2 ] ; then
    echo "Usage: bash Step1_Trim_fastq.sh [index_seq] [in_fastq] [SEQXE] [out_dir] [threads]"
    echo ""
    echo "Trimming the fastq file by TrimGalore program using index adaptor & SEQXE primer sequence"
    echo ""
    echo "arguments:"
    echo ""
    echo "index_seq     Illumina index adaptor sequence"
    echo "SEQXE			SEQXE adaptor sequence (optional)"
    echo "in_fastq		The path of fastq file"
    echo "out_dir	    The path for saving the trimmed fastq file"
    echo "threads		Number of threads to use (default: 8)"
    echo ""
    echo "For TruSeq Single Indexes, you can use index_seq as "
    echo "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
    echo ""
    exit 0
fi

echo "#scRepli-seq Step1 Adapter trimming	start: `date`"

#This index sequence is for Illumina TruSeq Single Indexes
index_seq=$1
#index_seq="GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"

in_fastq=$2
out_dir=$4
threads=$5

filename=`basename $in_fastq`
echo "#file name:${filename}"
prefix=${filename%.fastq.gz}

trim_galore --gzip -a ${index_seq} -o ${out_dir} -j ${threads} ${in_fastq} # Trimming with index sequence
mv ${out_dir}/${prefix}_trimmed.fq.gz ${out_dir}/${prefix}_index_trimmed.fq.gz # Renaming the trimmed file
mv ${out_dir}/${filename}_trimming_report.txt ${out_dir}/${prefix}_index_trimming_report.txt # Renaming the trimmed report file

echo "#scRepli-seq Step1 Adapter Trimming	end: `date`"
echo "#scRepli-seq Step1 SEQXE Trimming		start: `date`"


#Step1.2 Trim SEQXE primer
SEQXE=$3
#SEQXE="TGGTGTGTTGGGTGTGTTTCTGAAGNNNNNNNNN"

if [ -z "$SEQXE" ]
then

echo "#SEQXE Primer trimming is not performed.."

else

trim_galore --gzip -a ${SEQXE} -o ${out_dir} -j ${threads} ${out_dir}/${prefix}_index_trimmed.fq.gz # Trimming with SEQXE sequence
mv ${out_dir}/${prefix}_index_trimmed_trimmed.fq.gz ${out_dir}/${prefix}_index_SEQXE_trimmed.fq.gz # Renaming the trimmed file
mv ${out_dir}/${prefix}_index_trimmed.fq.gz_trimming_report.txt ${out_dir}/${prefix}_index_SEQXE_trimming_report.txt # Renaming the trimmed report file

fi

echo "#scRepli-seq Step1 SEQXE Trimming	end: `date`"