#!/bin/bash
#PBS -N Unzipping_combine_H4K20me1
#PBS -l select=1:ncpus=1:mem=59gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M email
#PBS -m ae

## This script is an example of what can be run for a single histone modification and using specific directory structures.

#Set histone (change in job name too)
histone="H4K20me1"

#FEMALE
#Set working directory to female directory
cd /path/to/dir/F_${histone}_treat_pileup/

#Unzip combine.tar.gz file
tar -xzvf combine.tar.gz

#Copy master tsv combine script into combine directory
cp /path/to/master/script/dir/tsv_combine_script.sh /path/to/dir/F_${histone}_treat_pileup/combine/

#Set working directory to combine
cd /path/to/dir/F_${histone}_treat_pileup/combine/

#Use sed to replace placeholders in the tsv combine script
sed -i "s/HISTONEMOD/F_${histone}/g" tsv_combine_script.sh
sed -i "s/generichistone/F_${histone}/g" tsv_combine_script.sh
sed -i "s/genericsex/Female/g" tsv_combine_script.sh

#Run tsv combine script
qsub tsv_combine_script.sh

#MALE
#Set working directory to male directory
cd /path/to/dir/M_${histone}_treat_pileup/

#Unzip combine.tar.gz file
tar -xzvf combine.tar.gz

#Copy master tsv combine script into combine directory
cp /path/to/master/script/dir/tsv_combine_script.sh /path/to/dir/M_${histone}_treat_pileup/combine/

#Set working directory to combine
cd /path/to/dir/M_${histone}_treat_pileup/combine/

#Use sed to replace placeholders in the tsv combine script
sed -i "s/HISTONEMOD/M_${histone}/g" tsv_combine_script.sh
sed -i "s/generichistone/M_${histone}/g" tsv_combine_script.sh
sed -i "s/genericsex/Male/g" tsv_combine_script.sh

#Run tsv combine script
qsub tsv_combine_script.sh