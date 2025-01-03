#!/bin/bash

#This script prepares the GWAS for the genetic correlation of the LDSC.
#It removes some columns, renames and orders them.

# Define GWAS
GWAS=$(cat $1)

# Define input directoru
INDIR=$2

# Define output directory
DIR="/homes/users/ajimenez/scratch/jpalm_gwas_sumstats"

# Create directory and if it already exists, the existing directory is used.
mkdir -p $DIR

#Critical step: check if the sumstats file is already were it belongs
if [ -e "${DIR}/${GWAS%.tsv}.sumstats.gz" ]; then
    #Write the check report
    echo "${GWAS%.tsv}.sumstats.gz already in directory" > ${GWAS%.tsv}.sumstats.info
else
    module load Miniconda3/4.9.2

    #Load ldsc module (hopefully this automatically loads 2-to-3 branch)
    source activate /homes/users/ajimenez/scratch/conda_env/ldsc

    #Define path to LDSC package
    PATH_TO_LDSC=/homes/users/ajimenez/scratch/ldsc

    #Define the HapMap3 list for the --merge_alleles
    HM3=/homes/users/ajimenez/scratch/gcat/ldscore/eur_w_ld_chr/w_hm3.snplist

    #Define the intermediate output file
    INTER=${GWAS%.tsv}_pre_sumstats.tsv

    #Sort the file
    { head -n 1 ${INDIR}/${GWAS}; tail -n +2 ${INDIR}/${GWAS} | sort -t$'\t' -k1,1; } > sorted_GWAS.txt

    awk -F'\t' -v OFS='\t' '
     NR==1 {
        for (i=1; i<=NF; i++) {
            if ($i == "rsID") rsID_col=i;
            if ($i == "n_samples") n_samples_col=i;
            if ($i == "pval") pval_col=i;
            if ($i == "z") z_col=i;
            if ($i == "ref") ref_col=i;                
            if ($i == "alt") alt_col=i;
            }
        }
    { 
        print $rsID_col, $n_samples_col, $pval_col, $z_col, $alt_col, $ref_col 
    }' sorted_GWAS.txt > $INTER
    
    #Rename the headers
    sed -i '1s/rsID/SNP/; 1s/n_samples/N/; 1s/pval/P/; 1s/z/Z/; 1s/alt/A1/; 1s/ref/A2/' $INTER

    #Call the munge_sumstats.py script from ldsc
    $PATH_TO_LDSC/munge_sumstats_modified.py \
    --sumstats $INTER \
    --out ${INTER%_pre_sumstats.tsv} \
    --merge-alleles $HM3
    #--ignore variant,se,beta,chr,pos

    #Move the file to the directory
    mv "${GWAS%.tsv}.sumstats.gz" "${DIR}/${GWAS%.tsv}.sumstats.gz"
    #Remove intermediate file and log files
    rm $INTER ${INTER%_pre_sumstats.tsv}.log sorted_GWAS.txt
    #Write the check report
    echo "${GWAS%.tsv}.sumstats.gz was created in specified directory" > ${GWAS%.tsv}.sumstats.info
    
fi
