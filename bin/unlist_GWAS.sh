#!/bin/bash

#This script takes the GWAS list and stores the path of each GWAS 
#in a separate file. This allows Nextflow to parallelise the following 
#processes.

#Exit immediately if any command within it returns a non-zero exit (error) status
set -e
#Treat references to unset variables as errors and make the script exit if so
set -u

#Set the input as the first argument passed to the script
gwas_list=$1

#Create separate text files for each gwas
awk 'NR > 1 {print > ("GWAS_" NR-1)}' "$gwas_list"
