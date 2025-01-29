#!/bin/bash
#PBS -N Combine_tsv_HISTONEMOD
#PBS -l select=1:ncpus=1:mem=119gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M email
#PBS -m ae

## This is master tsv combine script that is copied and modified by 4_OAN_ChIP_unzip_script_generic.sh

#Set histone and sex
histone="generichistone"
sex="genericsex"

#Set working directory
cd /path/to/dir/${histone}_treat_pileup/combine

#Combine all chromosome .tsv files into one
cat *.tsv > ${histone}_all_chroms.tsv

#Assign chromosome type
awk -v OFS='\t' '{
    split($1, a, ",");
    if (a[1] ~ /49\.1$|50\.1$|51\.1$|52\.1$|53\.1$/) 
        print $0, "X";
    else if (a[1] ~ /75\.1$|76\.1$|77\.1$|78\.1$|79\.1$/) 
        print $0, "Y";
    else if (a[1] ~ /0891\.1$/) 
        print $0, "mt";
    else 
        print $0, "A";
}' ${histone}_all_chroms.tsv > ${histone}_info.tsv

#Add histone mark and sex columns
awk -v OFS='\t' -v histonevalue="$histone" -v sexvalue="$sex" '{
    print $0, substr(histonevalue, 3), sexvalue;
}' ${histone}_info.tsv > ${histone}_info2.tsv

#Add tab separated headers
sed $'1iposition\tgene\tstrand\tgene_start\tgene_end\trelative_pos\theight\tc\tmod\tsex' ${histone}_info2.tsv > ${histone}_info_headers.tsv

#Remove unnecessary columns
cut -f6,7,8,9,10 -d$'\t' ${histone}_info_headers.tsv > ${histone}_info_simplified.tsv

#Delete intermediate files
rm ${histone}_all_chroms.tsv
rm ${histone}_info.tsv
rm ${histone}_info2.tsv
rm combined_*