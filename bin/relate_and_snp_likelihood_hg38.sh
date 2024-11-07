#!/bin/bash

#Load modules
module load relate/1.1.6-x86_64
module load Python/3.9.6-GCCcore-11.2.0

#Define the variables that don't change
POP=$3
PATH_TO_RELATE=/homes/users/ajimenez/scratch/relate
INDIR=/homes/users/ajimenez/scratch/gcat/relate/${POP}_hg38/anc_mut
MU=1.25e-8
COAL_FILE=/homes/users/ajimenez/scratch/gcat/relate/${POP}_hg38/coal/${POP}_hg38.coal
PATH_TO_PALM=/homes/users/ajimenez/scratch/palm

#Iterate over each line (row) of the provided batch of SNPs
while IFS=$'\t' read -r BP foo_1 foo_2 CHR VARIANT LD DAF rest_of_the_columns
do
    CHR_FILE=${POP}_chr${CHR}_upd
    INPUT_RELATE=${INDIR}/${CHR_FILE}
    OUTPUT_RELATE=${VARIANT}_out

    #Call the script and store it to obtain the SNP likelihood (max 5min to do so)
    timeout 5m ${PATH_TO_RELATE}/scripts/SampleBranchLengths/SampleBranchLengths_custom.sh -i $INPUT_RELATE -o $OUTPUT_RELATE \
    -m $MU --coal $COAL_FILE --format b --num_samples 5 --first_bp $BP --last_bp $BP --seed 1 || true

    #Remove the files that are not necessary
    rm *_out.anc
    rm *_out.mut
    rm *_out.dist

    #Define the input for the SNP Likelihood script
    INPUT_LIK=${OUTPUT_RELATE}.timeb

    #Check if the input of the SNP likelihood script exists and is not empty.
    #If there is no file, skip the rest of the iteration and start the next one    
    if [ ! -s "$INPUT_LIK" ]; then
        # If the file does not exist or is empty, move to the next iteration
        echo "Skipping iteration, $INPUT_LIK does not exist or is empty."
        continue
    fi

    #Create the necessary directory
    OUTDIR=${2}/ld_${LD}
    mkdir -p $OUTDIR
    #Define the ouput file name "root"
    OUTPUT_LIK=bp${BP}

    #Call the script and store it to obtain the SNP likelihood (max 5min to do so)
    timeout 5m python $PATH_TO_PALM/lik.py --times $OUTPUT_RELATE \
    --popFreq $DAF --out $OUTPUT_LIK --coal $COAL_FILE || true

    #Remove the files that are not necessary (again)
    rm *_out.timeb

    #Copy the output to the specified folder if the file exists and is not empty
    if [ -e "${OUTPUT_LIK}.quad_fit.npy" ] && [ -s "${OUTPUT_LIK}.quad_fit.npy" ]; then
        cp ${OUTPUT_LIK}.quad_fit.npy $OUTDIR/.
    fi

done < $1