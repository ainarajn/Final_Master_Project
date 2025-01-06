# Pipeline parameters 
  
  - genome = “hg19” (default) // Options: "hg19", "hg19tohg38", "hg38", "hg38tohg19" 
  - analysis = “PALM” (default) // Options: "PALM", "JPALM" 
  - GWAS\_List: txt file that contains the name of GWAS files (all the GWAS have to be previously downloaded), number of samples and parameter ref_effect. 
  - GWASdir: Path to the directory where the downloaded GWAS files are located. 
  - population: population to analyse. In our case, GCAT.
  - Relate\_SNP: rds file provided by the user.
  
  - VCfdir: VCF directory where all the VCF files for each chromosome are located (files provided).
  - HapMap3: files provided (one for hg19 and other for hg38).
  - LDblock: bed file provided (one for hg19 and other for hg38).
    
  - Max\_pvalue = 5e-8 
  - Min\_LDblocks = 25 
  - N\_per\_batch = 10 

# Processes 

- **Unlist**: This script takes the GWAS list and stores the name of each GWAS in a separate file. This allows Nextflow to parallelize the following processes. 
- **GWAS\_format**: This script currently handles txt, tsv, csv and vcf GWAS file types or gz/bgz versions of these file types. It  makes  the  GWAS  input  more  flexible  in  order  to  obtain  as  output  the  GWAS  with  the information of interest  for downstream analysis  and an established format. In addition, we select the SNPs from chr1:22 and with a p-value below the maxim p-value (Max\_pvalue). 
- **LD\_blocks**: It intends to add the LD\_block information in order to determine  if we continue to apply the downstream analysis based on the number of di erent LD blocks with significant SNPs.
- **LiftOver**: It creates a BED-formatted file for each selected SNPs from the GWAS file. Then executes the LiftOver tool and it substitutes the position in the GWAS with the new coordinates (hg19 to hg38). **Conditions**: (*hg19tohg38* or *hg38tohg19* genome version) & number of LD blocks above Min\_LDblocks.
- **Obtain\_N\_Comparisons**: This script creates a file with the GWAS that have 25 or more independent SNPs from all arguments passed and calculates the number of pairwise comparisons. **Conditions**: *JPALM* analysis.
- **munge\_sumstats**: This script prepares the GWAS for the genetic correlation of the LDSC. It removes some columns, renames and orders them. Finally, we call the munge_sumstats.py script from ldsc (https://github.com/bulik/ldsc). This repository loads Phyton2 environment, so we installed a pull request that ports to Phyton3 (#360 PR, belowlab:2-to-3) and we will work with this new branch. **Conditions**: *JPALM* analysis & number of LD blocks above Min\_LDblocks.
- **Genetic\_Correlation**: This script implements ldsc (https://github.com/bulik/ldsc) to obtain the genetic correlation between two GWAS traits. **Conditions**: *JPALM* analysis. 
- **SNP\_selection**: This script  intends to intersect the filtered GWAS SNPs with the SNPs present in Relate. **Conditions**: *PALM* analysis and if the genome version is either *hg19* or *hg38* it must implement an additional condition: number of LD blocks above Min\_LDblocks.
- **SNP\_selection\_JPALM**: This script intends to intersect the filtered GWAS SNPs with the SNPs present in Relate. After that, it only retains the SNP with the highest p-value per LD block. **Conditions**: *JPALM* analysis & (*hg19* or *hg38* genome version) & genetic correlation is above 0.2 & p-value of the genetic correlation is below 0.005.
- **SNP\_selection\_JPALM\_lifted**: This script intends to intersect the filtered GWAS SNPs with the SNPs present in Relate. After that, it only retains the SNP with the highest p-value per LD block. **Conditions**: *JPALM* analysis & (*hg19tohg38* or *hg38tohg19* genome version) & genetic correlation is above 0.2 & p-value of the genetic correlation is below 0.005.
- **Merge\_SNPs**: This process intents to merge all selected SNPs for all the GWAS inputs, then sorts them by variants and subsequently removes any duplicate entries.
- **Trim\_SNP**: This script not only trims SNPs by its p-value but  also separates each SNP so that they can be processed in future steps. This is incredibly helpful to avoid unnecessary sampling branch lengths and inferring SNPs likelihoods that will not be used later on.
- **Relate\_and\_SNP\_Likelihood**: This script iterates over each row of the provided batch of SNPs and implements Relate (https://github.com/MyersGroup/relate.git) to obtain the SNP likelihood.
- **marginal\_PALM**: This script implements PALM (https://github.com/andersonwinkler/PALM.git) to analyse individual traits and infer the direction and intensity of selection. **Conditions**: *PALM* analysis.
- **joint\_PALM**: This script implements joint PALM (implements PALM (https://github.com/andersonwinkler/PALM.git) to analyse two genetically correlated traits and infer the direction and intensity of selection. **Conditions**: *JPALM* analysis.
