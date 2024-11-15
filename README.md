This file represents the structure of the Nextflow script. It is not perfectly accurate to the final script because it does not contain all the required parameters needed for the scripts that are pending from adjustment. In addition, the input and output names are not accurate but I chose that description as to illustrate the intended files.

# PARAMETERS 

  - GWAS\_List: txt file that contains the name of GWAS files (all the GWAS have to be previously downloaded).  
  - GWASdir: Path to the directory where the downloaded GWAS files are located. 
  - GWASdir\_hg38: Path to the directory
  - genome = “hg19” (default) // Options: "hg19", "hg19tohg38", "hg38", "hg38tohg19" 
  - analysis = “PALM” (default) // Options: "PALM", "JPALM" 
  - ref_effect: TRUE (default). If beta effect is assigned to alt allele, change to FALSE. // Options: TRUE, FALSE
  - HapMap3: files provided (one for hg19 and other for hg38).
  - Relate\_SNP: rds file. If it is not provided by the user, it will be created through a “helper pipeline” that we will try to incorporate in the pipeline.  
  - VCfdir: VCF directory where all the VCF files for each chromosome are locared (these files are provided).
  - LDblock: bed file provided (one for hg19 and other for hg38).
  - LiftOver\_Tool: provided. 
  - LiftOver\_chain: provided (one for hg19tohg38 and other for hg38tohg19).
  - Max\_pvalue = 5e-8 
  - Min\_LDblocks = 25 
  - N\_per\_batch = 10 

# PROCESSES 

## **Unlist** -  Successfully incoporated into the pipeline.

**Script status:** DONE 

This script takes the GWAS list and stores the name of each GWAS in a separate file. This allows Nextflow to parallelize the following processes. 

**Input**: 

1.  GWAS\_List 

**Ouput**: 

- GWAS\_{num} 

**Script**: Unlist\_GWAS.sh

## **GWAS\_format** -  Successfully incoporated into the pipeline.

**Script status:** REVIEW

It  makes  the  GWAS  input  more  flexible  in  order  to  obtain  as  output  the  GWAS  with  the information of interest  for downstream analysis  and an established format. In addition, we select the SNPs from chr1:22 and with a p-value below the maxim p-value (Max\_pvalue). 

This script currently handles txt, tsv, csv and vcf GWAS file types or gz/bgz versions of these file types. 

**Input**: 

1. GWAS\_{num} [Unlist Output] 
2. GWASdir 
3. VCFdir 
4. Max\_pvalue
5. ref_effect
6. HapMap3

**Output**: 

- format\_{GWAS}.txt 

**Script**: gwas\_format2.R

## **LD\_blocks** -  Successfully incoporated into the pipeline.

**Script status:** DONE 

It intends to add the LD\_block information in order to determine  if we continue to apply the downstream analysis based on the number of di erent LD blocks with significant SNPs. 

**Input**: 

1. format\_{GWAS}.tsv [GWAS\_format Output] 
2. GWASdir 
3. LDblocks 

**Output**: 

- format\_{GWAS} +-{LDblocks} 

**Script**: ldblock\_filter.R

## **LiftOver** -  Successfully incoporated into the pipeline.

**Script status:** DONE 

It creates a BED-formatted file for each selected SNPs from the GWAS file. Then executes the LiftOver tool and it substitutes the position in the GWAS with the new coordinates (hg19 to hg38).

**Input**: 

1. format\_{GWAS}+-{LDblocks} [LD\_blocks Output] 
2. GWASdir 
3. LiftOver\_Tool 
4. LiftOver\_chain 

**Output**: 

- format\_{GWAS}\_hg38.tsv 

**Conditions**: This process is only applied for the **hg38** genome version and if the number of LDblocks obtained on **LDblock\_filter** is above the minimum established (Min\_LDblocks).  

**Script**: gwas\_liftover.sh 

## **Obtain\_N\_Comparisons** -  Successfully incoporated into the pipeline.

**Script status:** DONE 

This script creates a file with the GWAS that have >=25 independent SNPs from all arguments passed and calculates the number of pairwise comparisons. 

**Input**: 

1.  List GWAS filtered [LD\_blocks Output Collected] 

**Output**: 

- n\_pairwise\_comparisons.txt 

**Conditions**: This process is only applied for **JPALM** analysis. 

**Script**: n\_comparisons.sh 

## **Obtain\_Pairs\_batch** 

**Script status:** DONE 

This script intends to create a new file that contains each possible pairwise comparison in the GWAS list provided. 

**Input**: 

1.  List GWAS filtered [LD\_blocks Output Collected] 

**Output**: 

- comparison\_batch\_{num} 

**Conditions**: This process is only applied for **JPALM** analysis.

**Script**: gwas\_pairs\_batch\_obtention.sh

## **munge\_sumstats** 

**Script status:** PENDING ADJUSTMENTS 

This  script  prepares  the  GWAS  for  the  genetic  correlation  of  the  LDSC.  It  removes  some columns, renames and orders them.

It is necessary: <https://github.com/bulik/ldsc> (munge\_sumstats has been modified -> munge_sumstats_modified.py) 

**Input**: 

1.  format\_{GWAS}+-{LDblocks} [LD\_blocks Output] 

**Output**: 

- {GWAS}.sumstats.info 

**Conditions**: This process is only applied for **JPALM** analysis and if the number of LDblocks obtained on **LDblock\_filter** is above the minimum established (Min\_LDblocks).

**Script**: munge\_ukbb\_sumstats.sh

## **Genetic\_Correlation** 

**Script status:** PENDING ADJUSTMENTS 

This script implements LDSC (2 to 3 version) to obtain the genetic correlation between two GWAS traits. 

It is necessary: <https://github.com/bulik/ldsc>

**Input**: 

1. List sumstats.info [munge\_sumstats Output Collected] 
2. comparison\_batch\_{num} [Obtain\_Pairs\_batch Output Flatten] 
3. n\_comparisons.sh [Obtain\_N\_Comparisons Output] 

**Output**: 

- {GWAS1}+-{GWAS2}+-{value1}+-{result} 

**Conditions**: This process is only applied for **JPALM** analysis. 

**Script**: genetic\_correlation\_batch\_test.sh

## **SNP\_selection** 

**Script status:** DONE 

This script  intends to intersect the filtered GWAS SNPs with the SNPs present in Relate. After that, it only retains the SNP with the highest p-value per LD block.

**Input**: 

1. Relate\_SNP 
2. {GWAS}\_format. tsv / {GWAS}\_format \_hg38.tsv [GWAS\_format / LiftOver Output] 
3. GWASdir 
4. Max\_pvalue 

**Output**: 

- {GWAS}\_selected\_SNPs.tsv 

**Conditions**: This process is only applied for **PALM** analysis and if the number of LDblocks obtained on **LDblock\_filter** is above the minimum established (Min\_LDblocks).

**Script**: snp\_selection.R 

## **SNP\_selection\_JPALM** 

**Script status:** DONE 

This script intends to intersect the filtered GWAS SNPs with the SNPs present in Relate. After that, it only retains the SNP with the highest p-value per LD block. 

**Input**: 

1. {GWAS1}+-{GWAS2}+-{value1}+-{result} [Genetic\_Correlation Output Flatten ] 
2. Relate\_SNP 
3. GWASdir 
4. Max\_pvalue 

**Output**: 

- {GWAS}\_selected\_SNPs.tsv 

**Conditions**: This process is only applied for **hg19** genome version and **JPALM** analysis, and if the value 1 obtained on **Genetic\_Correlation** is above 0.2 and the result is below 0.005.** 

**Script**: snp\_selection\_jpalm.R

## **SNP\_selection\_JPALM\_hg38** 

**Script status:** DONE 

This script intends to intersect the filtered GWAS SNPs with the SNPs present in Relate. After that, it only retains the SNP with the highest p-value per LD block. 

**Input**: 

  val ready [LiftOver Output Collected] 
1. {GWAS1}+-{GWAS2}+-{value1}+-{result} [Genetic\_Correlation Output Flatten] 
2. GWASdir 
3. Relate\_SNP 
4. Max\_pvalue 

**Output**: 

- {GWAS}\_selected\_SNPs.tsv 

**Conditions**: This process is only applied for  **hg38** genome version and  **JPALM** analysis, and if the value 1 obtained on **Genetic\_Correlation** is above 0.2 and the result is below 0.005.

**Script**: snp\_selection\_jpalm\_hg38.R

## **Merge\_SNPs** 

**Script status:** DONE

**Input**: 

1.  List of Selected SNPs  [SNP\_selection / SNP\_selection\_JPALM Output Collected] 

**Output**: 

- merged\_SNPs.tsv 

**Script**: merge\_SNPs.sh

## **Trim\_SNP** 

**Script status:** DONE

This script not only trims SNPs by its p-value but  also separates each SNP so that they can be processed in future steps. This is incredibly helpful to avoid unnecessary sampling branch lengths and inferring SNPs likelihoods that will not be used later on.

**Input**: 

1. merged\_SNPs.tsv [Merge\_SNPs Output] 
2. Max\_pvalue 
3. N\_per\_batch 
4. analysis 

**Output**: 

- batch\_{num}\_ merged\_SNPs.tsv 

**Script**: trim\_snps\_batch.sh

## **Relate\_and\_SNP\_Likelihood** 

**Script status:** PENDING ADJUSTMENTS 

This script iterates over each row of the provided batch of SNPs and obtains the SNP likelihood.

It is necessary: <https://github.com/MyersGroup/relate.git> (SampleBranchLengths.sh has been modified -> SampleBranchLengths_custom.sh)

**Input**: 

1. batch\_{num}\_ merged\_SNPs.tsv [Trim\_SNP Output Flatten]
2. Likdir
3. Population
4. genome

**Output**: 

- bp{bp}.quad_fit.npy

**Script**: 

If the genome version is **hg19** then we use relate\_and\_snp\_likelihood.sh, and if it is **hg38** we apply relate\_and\_snp\_likelihood\_hg38.sh

## **marginal\_PALM** 

**Script status:** REVIEW LikDir

This script analyse individual traits and infer the direction and intensity of selection.

It is necessary: <https://github.com/andersonwinkler/PALM.git> (palm.py has been modified -> palm_custom.py)

**Input**: 

  val ready [Relate\_and\_SNP\_Likelihood Output Collected]
1. {GWAS}\_selected\_SNPs.tsv [SNP\_selection Output]
2. LikDir
3. Max\_pvalue

**Output**: 

- {GWAS}_marginal_PALM.txt

**Condition**: This process is only applied for **PALM** analysis

**Script**: palm.sh script


## **joint\_PALM** 

**Script status:** REVIEW LikDir

This script analyse two genetically correlated traits and infer the direction and intensity of selection.

It is necessary: <https://github.com/andersonwinkler/PALM.git> (palm.py has been modified -> palm_custom.py)

**Input**: 

  val ready [Relate\_and\_SNP\_Likelihood Output Collected]
1. {GWAS}\_selected\_SNPs.tsv [SNP\_selection\_JPALM / SNP\_selection\_JPALM\_hg38 Output]
2. LikDir
3. Max\_pvalue

**Output**: 

-  {GWAS}_J_PALM.txt
-  {GWAS}_significant_independent_SNPs.txt

**Condition**: This process is only applied for **JPALM** analysis

**Script**: jpalm.sh

# WORKFLOW
