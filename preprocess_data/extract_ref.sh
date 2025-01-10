#!/bin/bash
#SBATCH --job-name=extract_ref
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm-%j.out 
#SBATCH --partition=haswell
#SBATCH --mail-user=ainara.jimenez@upf.edu
#SBATCH --mail-type=ALL

##This script intends to add the reference column to GWAS summary statistics

#Load necessary packages
module load BCFtools 

#Define objects
##Define GWAS file
GWAS=/homes/users/ajimenez/scratch/gwas_tfm/GCST90012622_buildGRCh37.tsv.gz

##Extract the base name (without path and extension)
GWAS_BASE=$(basename "$GWAS" .tsv.gz)

##Define the column numbers where rsIDs are stored
ID_COLUMN=2

##Define an intermediate file to store the raw rsIDs
INTER=/homes/users/ajimenez/scratch/gwas_tfm/unique_rsIDs.txt

##Define file that stores rsIDs and reference allele columns
BASE_OUT=/homes/users/ajimenez/scratch/gwas_tfm/id_ref_${GWAS_BASE}

##Define output
output=/homes/users/ajimenez/scratch/gwas_tfm/${GWAS_BASE}_mod.tsv

#Define temp file
TEMP=/homes/users/ajimenez/gwas_tfm/gcat/temp.tsv

#Create an array of chromosomes (1 to 22)
CHROMOSOMES=( $(seq 1 22) )

##Loop to do the whole process per chromosome
for CHR in "${CHROMOSOMES[@]}"; do
  #Name the input file split by chromosome
  IN_CHR=${GWAS}_chr"$CHR".tsv
  #Create the new input file with just 1 chromosome
  gunzip -c "$GWAS" | awk -F '\t' -v chromosome="$CHR" '$1 == chromosome {print > output_file}' output_file="$IN_CHR"
  
  ##Make a file just with the list of unique rsIDs and query with bcftools
  #Select just the column with the rsID
  cut -f $ID_COLUMN $IN_CHR | \
  #Make just a unique list of rsID
  sort -u > $INTER
  #Name the appropriate dbSNP dataset to query from
  dbSNP=/gpfs/projects/lab_dcomas/1000genomes_phase3_dcomas/vcf/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

  #Name the output file
  OUT=${BASE_OUT}_chr${CHR}.tsv
  #Query SNP database for rsIDs, retrieving REF allele
  bcftools query -f '%ID\t%REF\n' -i 'ID=@/homes/users/ajimenez/scratch/gwas_tfm/unique_rsIDs.txt' $dbSNP > $OUT
  
  #Add header to the output file
  sed -i '1i rsID\tref' $OUT
done

##Merge all the datasets into one
#Name the final output file
REF_ID=${BASE_OUT}.tsv

#Set a condition for the header situation
HEADER_WRITTEN=false
# Create an empty file
> "$REF_ID"

# Iterate over the file names and concatenate them
for ((i=1; i<=22; i++))
do
    INPUT_FILE=${BASE_OUT}_chr${i}.tsv
    # Check if it's the first file and write the header
    if [[ "$HEADER_WRITTEN" = false ]]; then
        cat "$INPUT_FILE" >> "$REF_ID"
        HEADER_WRITTEN=true
    else
        # Skip the header and append the contents
        tail -n +2 "$INPUT_FILE" >> "$REF_ID"
    fi
done

##Clean-up
rm $INTER
rm ${BASE_OUT}_chr*
rm ${IN}_chr*

#Descompress GWAS file
GWAS_desc=/homes/users/ajimenez/scratch/gwas_tfm/${GWAS_BASE}.tsv
gunzip -c "$GWAS" > "$GWAS_desc"

#Add REF allele information to GWAS data
awk -v OFS="\t" '
    BEGIN {
        #Read REF allele mappings into an associative array
        while ((getline < "'$REF_ID'") > 0) {
            variant_map[$1] = $2
        }
    }
    NR == 1 {
        #Print original header and add "ref" column
        print $0, "ref"
    }
    NR > 1 {
        #Print the row with REF if mapping exists, otherwise "NA"
        if ($2 in variant_map) {
            print $0, variant_map[$2]
        } else {
            print $0, "NA"  # Colocar "NA" si no hay mapeo
        }
    }
' $GWAS_desc > $output

##Clean-up
rm $GWAS_desc
#rm $REF_ID