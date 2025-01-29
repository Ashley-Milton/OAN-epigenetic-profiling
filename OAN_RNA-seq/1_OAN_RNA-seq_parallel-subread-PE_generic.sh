#!/bin/bash
#PBS -N Subread_align_paired_OAN_RNA-seq
#PBS -l select=1:ncpus=8:mem=120gb
#PBS -l walltime=6:00:00
#PBS -j oe
#PBS -WMail_Users=email
#PBS -m ae

### Make sure to first build an index for your genome, cmd below
#subread-buildindex -o /path/to/index/basename /path/to/genome.fasta

### set the following variables (line 19-24)
# Directory to folder containing your sequences
# Directory to output folder
# Path to your index
# File extension identifier for your R1 and R2 file
# your input sequence type, 1 for genomic DNA-seq and 0 for RNA-seq

indir="/path/to/dir/with/raw_data"
outdir="/path/to/outdir"
index="/path/to/subread/index"
R1identifier="R1_001.fastq.gz"
R2identifier="R2_001.fastq.gz"
seqtype="0"

### NOTE, make sure you see all of your raw R1 input files with the following ls command in your ${indir}
# ls *${R1identifier}

############################################################################################################################
module load subread/2.0.2 parallel

cd ${indir}

for f in $(ls *${R1identifier}); do echo ${f/${R1identifier}/}; done \
| parallel --jobs ${NCPUS} \
subread-align -t ${seqtype} -T ${NCPUS} --sortReadsByCoordinates -i ${index} \
-r {}${R1identifier} \
-R {}${R2identifier} \
-o ${outdir}/{}_subread.bam
