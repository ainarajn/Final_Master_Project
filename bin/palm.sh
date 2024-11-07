#!/bin/bash

module load Python/3.9.6-GCCcore-11.2.0

PATH_TO_PALM=/homes/users/ajimenez/scratch/palm
GWAS=$1
PATH_TO_LIK=$2
MAXPVAL=$3

#Call the script
python $PATH_TO_PALM/palm_custom.py \
    --traitDir ${PATH_TO_LIK}/ \
    --metadata $GWAS \
    --maxp $MAXPVAL \
    --B 10000 \
    > $(basename "$1" _selected_SNPs.tsv)_marginal_PALM.txt