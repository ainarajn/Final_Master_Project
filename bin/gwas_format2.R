#!/usr/bin/env Rscript

#Suppress only automatic warnings from R
options(warn = -1)

#Packages
suppressMessages(require(dplyr))
suppressMessages(require(stringr))
suppressMessages(require(tidyverse))
suppressMessages(require(stringi))
suppressMessages(require(vcfR))
suppressMessages(require(data.table))

#Put the arguments passed to the script into a vector
files <- commandArgs(trailingOnly = TRUE)
gwas_name <- read.delim(files[1], header = FALSE)
gwas_dir <- files[2]
gwas_path <- paste0(gwas_dir, "/", gwas_name[1])
vcf <- files[3]
maxp <- as.numeric(files[4])
ref_effect <- as.logical(files[5])
hapmap3 <- files[6]

#Define output name and path
output_name <- paste0("format_", sub("\\..*", "", gwas_name), ".tsv")
output_path <- paste0(gwas_dir, "/", output_name)

if (file.exists(output_path)) {
  #If the file already exist, emit the expected output
  writeLines(output_name, sub("\\.tsv$", ".txt", output_name))
} else {
  #Function to load GWAS according extension
  load_GWAS <- function(data) {
    extension <- sub(".*\\.(.+)\\.(gz|bgz)$|.*\\.(.+)$", "\\1\\3", basename(data))
    switch(extension,
           "tsv" = fread(data, header = TRUE, sep = "\t"),
           "csv" = fread(data, header = TRUE, sep = ","),
           "txt" = fread(data),
           "vcf" = {
             vcf <- read.vcfR(data)
             data <- as.data.frame(vcf@fix)
             format <- as.data.frame(extract_gt_tidy(vcf, format_fields = c("ES", "SE", "LP","SS"), 
                                                     gt_column_prepend = "", alleles = FALSE))
             cbind(data, format)[, -c(6:10)]
           }
    )
  }
  
  #Load GWAS
  gwas <- load_GWAS(gwas_path)
  
  #Function to change header
  change_headers <- function(data) {
    #Possible names for the header
    variant <- c("^(variant|uniqID|variant_ID|markername|SNPID|ID|SNP|cptid)$")
    beta <- c("^(beta|effect|effects|b|effect_size|es)$")
    se <- c("^(se|standard_error|StandardError|StdErr|sebeta)$")
    pval <- c("^(pval|p_val|p-val||p.val|pvalue|p_value|p-value|p.value|gc_pvalue|p)$")
    LP <- c("^(lp|neg_log_10_p_value|pvalue_mlog|mlogp)$")
    chr <- c("^(chr|chromosome|chrom|snp-chr|chr_ID)$")
    pos <- c("^(pos|base_pair_location|snppos|snp_pos|bp|position|chr_pos)$")
    z <- c("^(z|zscore|z-score|GC_zscore|tstat|t_stat|t-statistic)$")
    n_samples <- c("^(N|n_complete_samples|sample_size|n_total|TotalSampleSize|SS)$")
    
    #If a column presents rsID, change the column name for rsID
    colnames(data) <- ifelse(apply(head(data, 50), 2, function(col) { 
      any(stri_detect_regex(col, "^rs", case_insensitive = TRUE))
    }), "rsID", colnames(data))
    
    #Change column names to a standard format
    colnames(data) <- stri_replace_all_regex(colnames(data),
                                             pattern = c(variant, beta, se, z,
                                                         pval, LP, chr, pos, n_samples),
                                             replacement = c("variant", "beta", "se", "z",
                                                             "pval", "LP", "chr", "pos", "n_samples"),
                                             vectorize = FALSE, case_insensitive = TRUE)
    #Determine beta as numeric variables
    data$beta <- as.numeric(data$beta)
    return(data)
  }
  
  #Add function
  gwas <- change_headers(gwas)
  
  #If it is necessary, transform -log10(p-value) to p-value
  if (!("pval" %in% colnames(gwas)) && "LP" %in% colnames(gwas)) {
    gwas$LP <- as.numeric(gwas$LP)
    gwas$pval <- 10^(-gwas$LP)
  }
  
  #Removal of SNPs above max. p-value
  gwas <- gwas[gwas[["pval"]] < maxp, ]
  
  #Verify there are still observations in the data
  if (nrow(gwas) == 0) {
    message("In gwas ", gwas_name, " there is no SNPs with p-value below the maximum p-value indicated (", maxp, ").")
    stop("Process halted: No observations to analyze.")
  }
  
  #If it is necessary, calculate SE
  if (!("se" %in% colnames(gwas))) {
    if ("z" %in% colname(gwas)) {
      gwas$z <- as.numeric(gwas$z)
      gwas$se <- gwas$beta / gwas$z
    } else {
      gwas$se <- 0
    }
  }
  
  #If is necessary, calculate Z
  if (!("z" %in% colnames(gwas)) && "se" %in% colnames(gwas) && "se" != 0) {
    gwas$se <- as.numeric(gwas$se)
    gwas$z <- gwas$beta / gwas$se
  }
  
  ################# REVIEW THE FUNCTIONS RELATED TO OBTAIN REF AND ALT ALLELES
  #Function for change header names in allele columns
  alleles_f <- function(data, ref_effect) {
    beta_ref <- c("^(ref|effect_allele|A1|allele1|allele_1|reference_allele|inc_allele|EA)$")
    beta_alt <- c("^(alt|other_allele|A2|allele2|allele_2|non_effect_allele|dec_allele|NEA)$")
    colnames(data) <- stri_replace_all_regex(colnames(data),
                                             pattern = c(beta_ref, beta_alt),
                                             replacement = c("ref", "alt"),
                                             vectorize = FALSE, case_insensitive = TRUE)
    
    #If the parameter ref_effect is FALSE, swap the column names
    if (!ref_effect) {
      ref_idx <- match("ref", colnames(data))
      alt_idx <- match("alt", colnames(data))
      colnames(data)[ref_idx] <- "alt"
      colnames(data)[alt_idx] <- "ref"
    }
    
    return(data)
  }
  
  #Determine Reference and Alternative alleles
  ##GWAS CATALOG
  if ("effect_allele" %in% colnames(gwas) & "ref_allele" %in% colnames(gwas)) {
    gwas$ref <- ifelse(gwas$ref_allele == "EA", gwas$effect_allele, ifelse(gwas$ref_allele == "OA", gwas$other_allele, NA))
    gwas$alt <- ifelse(gwas$ref_allele == "EA", gwas$other_allele, ifelse(gwas$ref_allele == "OA", gwas$effect_allele, NA))
    ##For the rest of cases
  } else {
    gwas <- alleles_f(gwas, ref_effect)
  }
  
  #################
  
  #Create function to add missing columns to the GWAS from SNPdb file
  add_from_rsID <- function(gwas, vcf, missing_columns) {
    #Load SNPdb
    SNPdb <- getFIX(read.vcfR(vcf))[, -c(6, 7)]
    #Filter SNPdb to obtain info about the variants from the GWAS
    SNPdb_filtered <- SNPdb[SNPdb[,"ID"] %in% unique(gwas$rsID),]
    #Data frame with relevant info
    info_SNP <- data.frame(
      rsID = SNPdb_filtered[, "ID"],
      chr =  SNPdb_filtered[, "CHROM"],
      pos = SNPdb_filtered[, "POS"],
      ref = SNPdb_filtered[, "REF"],
      alt = SNPdb_filtered[, "ALT"],
      stringsAsFactors = FALSE)
    #Filter with missing columns
    info_SNP <- info_SNP %>%
      dplyr::select(c("rsID", all_of(missing_columns)))
    #Add missing columns to the GWAS
    gwas <- gwas %>%
      left_join(info_SNP, by = "rsID")
    return(gwas)
  }
  
  #Create function to add missing columns from variant column or from the function previously created for rsID
  #It also discard SNPs from chr X, Y, MT and variants that are no 1 substitution base, and create a variant column. 
  variant_function <- function(gwas) {
    columns_variants <- c("chr", "pos", "ref", "alt")
    missing_columns <- setdiff(columns_variants, colnames(gwas))
    if (length(missing_columns) != 0) {
      if (all(grepl("^([0-9]+|X|Y|MT|M):[0-9]+:[A-Za-z]+:[A-Za-z]+$", head(gwas$variant, 10)))) {
        split_variant <- str_split_fixed(gwas$variant, ":", 4)
        gwas$chr <- split_variant[, 1]
        gwas$pos <- split_variant[, 2]
        gwas$ref <- split_variant[, 3]
        gwas$alt <- split_variant[, 4]
        gwas$variant <- paste0(split_variant[, 1], ":", split_variant[, 2])
      } else if (any(grepl("^[0-9]+:[0-9]+$", head(gwas$variant, 10))) && identical(missing_columns, c("chr", "pos"))) {
        split_variant <- str_split_fixed(gwas$variant, ":", 2)
        gwas$chr <- split_variant[, 1]
        gwas$pos <- split_variant[, 2]
        gwas$variant <- paste0(split_variant[, 1], ":", split_variant[, 2])
      } else if ("rsID" %in% colnames(gwas)) {
        gwas <- add_from_rsID(gwas, vcf, missing_columns)
      }
    }
    gwas$pos <- as.numeric(gwas$pos)
    
    #Discard SNPs from chr X, Y, MT and variants that are no 1 substitution base
    gwas$chr <- as.numeric(gwas$chr)
    gwas <- gwas[!is.na(gwas$chr), ]
    gwas <- gwas[nchar(gwas$ref) == 1 & nchar(gwas$alt) == 1, ]
    
    #Add column variant if it is necessary
    if(!("variant" %in% colnames(gwas))) {
      gwas$variant <- paste0(gwas$chr, ":", gwas$pos)
    }
    return(gwas)
  }
  
  #Add function
  gwas <- variant_function(gwas)
  
  #Verify there are still observations in the data
  if (nrow(gwas) == 0) {
    message("In gwas ", gwas_name, " there is no SNPs to analyze.")
    stop("Process halted: No observations to analyze.")
  }
  
  #Function to add rsID from chromosome and position column.
  add_rsID = function(data){
    if (!("rsID" %in% colnames(data))) {
      HM3 <- fread(hapmap3)
      data <- data %>%
        left_join(HM3, by = c("variant"))
    }
  }
  
  #Add function
  gwas <- add_rsID(gwas)
  
  #Verify that the GWAS includes all the columns required for the downstream analysis.
  columns_req = c("rsID", "variant", "n_samples", "beta", "se", 
                  "pval", "z", "ref", "alt", "chr", "pos")
  missing_columns = setdiff(columns_req, colnames(gwas))
  
  if (length(missing_columns) > 0) {
    message("There is one or more columns missing, revise ", gwas_name, " file to verify it. In case there is no missing column, please adjust the header to match the script requirements.")
    stop("Error: The following columns are missing: ", paste(missing_columns, collapse = ", "))
  } else {
    #Create GWAS output with the established format
    gwas_format <- gwas[, c("variant", "n_samples", "beta", "se", "pval",
                            "z", "ref", "alt", "chr", "pos", "rsID")]
  }
  
  #Remove useless data
  rm(gwas)
  
  #Create output file
  write_tsv(gwas_format, file = paste0(gwas_dir, "/", output_name))
  
  #Pass the GWAS name (ended in tsv) to the next process in a txt file
  writeLines(output_name, sub("\\.tsv$", ".txt", output_name))
}

