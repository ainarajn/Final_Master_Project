#!/bin/bash

module load Python/3.9.6-GCCcore-11.2.0

PATH_TO_PALM=/homes/users/ajimenez/scratch/palm
GWAS=$1
MAXPVAL=$2
GENOME=$3
POP=$4
PATH_TO_LIK=/homes/users/ajimenez/scratch/SNP_likelihoods_${GENOME}
OUTDIR=/homes/users/ajimenez/scratch/results_marginal_palm_${GENOME}_${POP}

# Create directory and if it already exists, the existing directory is used.
mkdir -p $OUTDIR

#Call the script
python $PATH_TO_PALM/palm_custom.py \
    --traitDir ${PATH_TO_LIK}/ \
    --metadata $GWAS \
    --maxp $MAXPVAL \
    --B 10000 \
    > $(basename "$1" _selected_SNPs.tsv)_marginal_PALM.txt

cp "$(basename "$1" _selected_SNPs.tsv)_marginal_PALM.txt" "${OUTDIR}/."
