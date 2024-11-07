#!/bin/bash

#This script prepares the GWAS for the genetic correlation of the LDSC.
#It removes some columns, renames and orders them.

# Define GWAS
GWAS=$(cat $1)

# Define output directory
DIR="${projectDir}/jpalm_gwas_sumstats"

# Create directory and if it already exists, the existing directory is used.
mkdir -p $DIR

#Critical step: check if the sumstats file is already were it belongs
if [ -e "${DIR}/${GWAS%.tsv}.sumstats.gz" ]; then
    #Write the check report
    echo "${GWAS%.tsv}.sumstats.gz already in directory" > ${GWAS%.tsv}.sumstats.info
else
    #Load ldsc module (hopefully this automatically loads 2-to-3 branch)
    source activate /homes/users/ajimenez/scratch/conda_env/ldsc

    #Define the file needed to add the rsID to the GWAS SNPs that happen to
    #be on the HapMap3 SNP list. It's the one that the pre-computed LDSC have
    #anyway so not a loss of information whatsoever.
    HAPMAP_SNPs=/homes/users/ajimenez/scratch/gcat/ldscore/hm3_SNPs_variant_id.tsv

    #Sort the file for better efficiency
    { head -n 1 $HAPMAP_SNPs; tail -n +2 $HAPMAP_SNPs | sort -t$'\t' -k1,1; } > sorted_HAPMAP_SNPs.txt

    #Define path to LDSC package
    PATH_TO_LDSC=/homes/users/ajimenez/scratch/ldsc

    #Define the HapMap3 list for the --merge_alleles
    HM3=/homes/users/ajimenez/scratch/gcat/ldscore/eur_w_ld_chr/w_hm3.snplist

    #Define the intermediate output file
    INTER=${GWAS%.tsv}_pre_sumstats.tsv

    #Sort the file
    { head -n 1 $GWAS; tail -n +2 $GWAS | sort -t$'\t' -k1,1; } > sorted_GWAS.txt
    #Add the rsID of the HapMap SNPs to the GWAS using the "variant" column (the 1st one) as the key
    awk -F'\t' 'BEGIN {OFS = "\t"} NR == FNR {data2[$1] = $2; next} {print $0, data2[$1]}' sorted_HAPMAP_SNPs.txt sorted_GWAS.txt > $INTER

    #Rename the headers (for linux)
    sed -i '1s/rsID/SNP/; 1s/n_complete_samples/N/; 1s/tstat/Z/; 1s/alt/A1/; 1s/ref/A2/' $INTER

    #Call the munge_sumstats.py script from ldsc
    $PATH_TO_LDSC/munge_sumstats_modified.py \
    --sumstats $INTER \
    --out ${INTER%_pre_sumstats.tsv} \
    --merge-alleles $HM3 \
    --ignore variant,se,beta,chr,pos

    #Move the file to the directory
    mv "${GWAS%.tsv}.sumstats.gz" "${DIR}/${GWAS%.tsv}.sumstats.gz"
    #Remove intermediate file and log files
    rm $INTER ${INTER%_pre_sumstats.tsv}.log sorted_HAPMAP_SNPs.txt sorted_GWAS.txt
    #Write the check report
    echo "${GWAS%.tsv}.sumstats.gz was created in specified directory" > ${GWAS%.tsv}.sumstats.info
fi
