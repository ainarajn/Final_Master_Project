Welcome! 

This repository contains the scripts and the pipeline developed and used for my master's thesis, titled "A novel pipeline to detect direct and conditional polygenic selection". The work presented here focuses on omic data analysis and includes a pipeline designed to automate polygenic selection analyses for virtually any trait to which a GWAS has been performed.

It should be noted that currently, the repository cannot be replicated, as all created directories’ paths are set in my cluster working directory.

## Reposity structure
The scripts are organized into the following directories:
  - bin: Contains the scripts used in the Nextflow processes.
  - modified_scripts: Contains scripts originally obtained from other GitHub repositories (referenced in section 3.6 of the thesis report) that were modified by the Evolutionary Population Genetics Lab.
  - preprocess_data: Includes the script used to incorporate the reference allele column in GWAS summary statistics that required preprocessing (the script requires modification for each GWAS dataset).


## Overview of the pipeline's workflow structure
Detailed explanation in section 4.1.2. of the thesis report.

<p align="center">
  <img width="697" src="https://github.com/user-attachments/assets/9505650a-f80e-420e-923f-16b902952a5f" alt="Workflow diagram of the developed pipeline" />
</p>

## Inputs
To execute this pipeline, the user will provide the following inputs:

  - **GWAS_List**: Path of a .txt file that contains three columns with header. The first column corresponds to the GWAS file name, the second column corresponds to the number of samples and the third column indicates if the parameter ref_effect is TRUE, FALSE or NULL (the list structure is exemplified in Annex 2). When ref_effect is set to TRUE, the effect is associated with the reference allele, whereas if it is set to FALSE, it is associated with the alternative. When the effect can be associated with either, depending on the SNP, the parameter is set to NULL. In this previous case, the GWAS data must include an additional column specifying the reference allele. If none of these conditions are met, the user must preprocess the data by adding the additional column to define the reference allele, ensuring a correct analysis.
    Additionally, it should be noted that all the listed GWAS need to be previously downloaded and are required to be assembled to the same genome build (hg19 or hg38).
  - **GWASdir**: Path to the directory where the downloaded GWAS files are located.
  - **genome**: It will indicate the genome build and if it is necessary to implement the lift over. Currently, the available options for this parameter are: “hg19”, “hg19tohg38”, “hg38” and “hg38tohg19”.
    In cases where the GWAS data is missing the chromosome and position columns, the information will be obtained from VCF which is based on the genome build hg19. Therefore, the user must specify for the genome parameter “hg19” or “hg19tohg38”.
  - **analysis**: It will indicate the PALM method used to infer selection (marginal for individual traits or joint for correlated traits). Therefore, the available options for this parameter are: “PALM” and “JPALM”.
  - **population**: It will indicate the population to be analysed by the user.
  - **Relate_SNPs**: Path of the Relate file (.rds) with the corresponding population to analyse. It should be noted that the Relate file provided has to build in the corresponding genome, so if the genome is “hg19” or “hg38tohg19” the ARGs have to be built in hg19, while if the genome is “hg38” or “hg19tohg38” they have to be built in hg38.

The rest of parameters used as inputs will not have to be provided by the user. Rather, once the pipeline is fully operational and replicable, they will be provided by a Docker container that will present all the dependencies required.
