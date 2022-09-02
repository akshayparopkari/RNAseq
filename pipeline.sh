#!/usr/bin/bash

# For reliable, robust and maintainable bash scripts, start with following commands
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
# set -euo pipefail
# IFS=$"\n\t"

#########################################################################################
# HEADER
#########################################################################################
#%
#% DESCRIPTION
#
#% This bash file that stiches together the preprocessing steps for RNA-seq sequences
#% from 3"RNA-seq project.
#%
#% USAGE and EXAMPLE
#% sh ${SCRIPT_NAME} [-h] inputs max_thread_count
#% ./${SCRIPT_NAME} [-h]  inputs max_thread_count
#%
#% OPTIONS
#% inputs: Directory containing all raw sequence FASTQ files, this will be updated in the
#%         next version
#% threads: Maximum number of threads to use for processing samples. Maximum thread
#%          count can be calculated as the product of Thread(s) per core, Core(s) per
#%          socket, and Socket(s) obtained from lscpu command.
#
#% NOTE: Care should be taken to supply full or relative path to the script file (this
#%       file) and input directory.
#%
# ---------------------------------------------------------------------------------------
#+ SCRIPT INFORMATION
#
#+ VERSION: 0.6.0
#+ AUTHOR:  Akshay Paropkari
#+ LICENSE: MIT
#
# ---------------------------------------------------------------------------------------
#% REFERENCES
#
#% 1. BBTools Suite - https://jgi.doe.gov/data-and-tools/bbtools/
#% 2. FastQC - http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#% 3. STAR - manual: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
#%
#########################################################################################
# END_HEADER
#########################################################################################

# Help output
#== needed variables ==#
SCRIPT_HEADSIZE=$(head -50 "${0}" |grep -n "^# END_HEADER" | cut -f1 -d:)
SCRIPT_NAME=$(basename "${0}")

#== usage functions ==#
usagefull() { head -"${SCRIPT_HEADSIZE}" "${0}"| grep -e "^#%" | sed -e "s/^#[%+]\ //g" -e "s/\${SCRIPT_NAME}/${SCRIPT_NAME}/g" ; }
scriptinfo() { head -"${SCRIPT_HEADSIZE}" "${0}" | grep -e "^#+" | sed -e "s/^#+\ //g" -e "s/\${SCRIPT_NAME}/${SCRIPT_NAME}/g"; }

if [[ "${1}" = "-h" ]] || [[ "${1}" = "--help" ]] ; then
    usagefull
    scriptinfo
    exit 0
fi


# Testing for input
if [[ -z "${1}" ]] ; then
    echo -e "Input directory not supplied. Please supply a directory with raw FASTQ files"
    exit 1
fi
if [[ -z "${2}" ]] ; then
    echo -e "Maximum thread count to use is not specified. Please provide a value for number of threads."
    exit 1
fi


# Changing working directory to input directory
# print the directory name which is being processed
echo -e "\n$(date '+%a %D %r') | Input directory identified as ${1}"

# cd into each directory in the directory
cd "${1}" && echo -e "$(date '+%a %D %r') | Switched into ${1}" || { echo -e "$(date '+%a %D %r') | cd into input directory failed! Please check your working directory."; exit 1; }

# Collect adapter files
truseq=$(find ~ -type f -name "*truseq.fa.gz" -exec realpath {} \;) && echo -e "$(date '+%a %D %r') | Identified TruSeq Fasta file - ${truseq}"
polya=$(find ~ -type f -name "*polyA.fa.gz" -exec realpath {} \;) && echo -e "$(date '+%a %D %r') | Located polyA Fasta file - ${polya}"
genome_dir=$(find ~ -type d -name "ca_genome" -exec realpath {} \;) && echo -e "$(date '+%a %D %r') | Located genome directory - ${genome_dir}"

# Starting preprocessing
for f in $(find . -type f -name "*.fastq*")
do

  in_file=$(basename "${f}")
  sample_id=$(basename "${f}" | cut -d_ -f1)

  echo -e "$(date "+%a %D %r") | Processing ${sample_id}"

  # •••••••••••••••••••••••••••••••••• #
  # Step 1: Adapter and polyA trimming #
  # •••••••••••••••••••••••••••••••••• #

  echo -e "$(date "+%a %D %r") | Adapter and polyA trimming"

  # Create output file name and command to run
  trimmed_file="${sample_id}"_trimmed.fastq
  CMD1="bbduk.sh in=$in_file out=$trimmed_file ref=$truseq,$polya k=13 ktrim=r mink=5 qtrim=r trimq=20 minlength=20"

  # Echo and run command
  echo "$(date "+%a %D %r") | ${CMD1}"
  $CMD1

  # •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••• #
  # Step 2: Fastqc generates a report of sequencing read quality #
  # •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••• #

  echo -e "$(date "+%a %D %r") | Read quality assessment"

  # Create directory to save the quality reports, one per FASTQ file and create command to run
  mkdir -p "${sample_id}"_qc
  CMD2="fastqc -t ${2} --nogroup -o ${sample_id}_qc ${trimmed_file}"

  # Echo and run command
  echo "$(date "+%a %D %r") | ${CMD2}"
  $CMD2

  # ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••• #
  # Step 3: STAR aligns to  the reference genome and calculates gene counts #
  # ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••• #

  echo -e "$(date "+%a %D %r") | Aligning QC reads to Candida albicans A21 genome"
  # Create command to run
  sample_id=$(basename "${f}" .fastq | cut -d_ -f1)
  CMD3="STAR --runThreadN ${2} --genomeDir ${genome_dir} --readFilesIn ${trimmed_file} --outFilterType BySJout --outFilterMultimapNmax 25 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.3 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix ${sample_id} > $(basename "${f}" .fastq)_alignment.log"

  # Echo and run command
  echo "$(date "+%a %D %r") | ${CMD3}"
  $CMD3

done

# •••••••••••••••••••••••••••••••• #
# Collect all counts into one file #
# •••••••••••••••••••••••••••••••• #
echo -e "$(date "+%a %D %r") | Collecting gene counts for all samples"

# Create command to run
format_counts=$(find ~ -type f -name "format_counts_table.py" -exec realpath {} \;) && echo -e "$(date '+%a %D %r') | Located gene count table formatting script - ${format_counts}"
CMD4="${format_counts} ${1} -o ${1}"

# Echo and run command
echo "$(date "+%a %D %r") | ${CMD4}"
$CMD4

# •••••••••••••••••••••••••••••• #
# Organize files created by STAR #
# •••••••••••••••••••••••••••••• #

# Move log files into a directory
mkdir -p ./STAR_log ./trim_log
mv -t STAR_log/ ./*Log.out ./*Log.progress.out ./*Log.final.out ./*.out.bam ./*.out.tab  && echo -e "$(date '+%a %D %r') | Moved STAR intermediate files"
mv -t trim_log/ *_trimmed.fastq && echo -e "$(date '+%a %D %r') | Moved intermediate adapter trimming files"

# •••••••••••••••••••••••••••••••••••••••••••••••• #
# Run differential expression analysis with DESeq2 #
# •••••••••••••••••••••••••••••••••••••••••••••••• #
echo -e "$(date "+%a %D %r") | Running differential expression analysis using DESeq2"
METADATA=$(find . -type f -name "*etadata*" -exec realpath {} \;) && echo -e "$(date '+%a %D %r') | Located metadata file - ${METADATA}"
GENE_COUNTS=$(find . -type f -name "gene_raw_counts.txt" -exec realpath {} \;) && echo -e "$(date '+%a %D %r') | Located gene count matrix file - ${GENE_COUNTS}"
DESEQ2=$(find ~ -type f -name "deseq.R" -exec realpath {} \;) && echo -e "$(date '+%a %D %r') | Located DESEq2 Rscript file - ${DESEQ2}"
CMD5="${DESEQ2} ${GENE_COUNTS} ${METADATA} ./deseq2_lfc.txt ./MA_plot.pdf"

# Echo and run command
echo "$(date "+%a %D %r") | ${CMD5}"
$CMD5

## REMOVE ME
rm -rv *log *qc && echo -e "$(date "+%a %D %r") | Removed intermediate log and QC files"

