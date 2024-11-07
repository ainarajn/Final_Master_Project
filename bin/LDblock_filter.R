#!/usr/bin/env Rscript

#Packages
require(dplyr)
require(stringr)
require(tidyverse)
require(data.table)

#Put the arguments passed to the script into a vector
files <- commandArgs(trailingOnly = TRUE)

# Load GWAS and LD blocks
# gwas_name = readLines(files[1]) # Ouput GWASformat is txt file
gwas_name <- files[1]
gwasdir <- files[2]
gwas_path <- paste0(gwasdir, "/", gwas_name)
gwas <- fread(gwas_path, header = TRUE, sep = "\t")
ld <- fread(files[3])

#Function to add LD blocks
add_LD_blocks <- function(data) {
  #Remove the "chr" part from the chr column
  ld$chr <- as.numeric(gsub("chr", "", ld$chr))
  #Create a column with the number of the LD block (the row number)
  ld$LD_block <- row_number(ld)
  #Define chromosomes
  chromosomes <- 1:22
  #Create for loop to iterate over chromosomes
  for (i in chromosomes) {
    start_positions <- ld$start[ld$chr==i]
    data$LD_start[data$chr==i] <-
      start_positions[(findInterval(data$pos[data$chr==i],
                                    ld$stop[ld$chr==i] - 1e-10) + 1)] #The +1 is added because otherwise it selects the previous interval
    #Match the positions to obtain the LD block
    data$LD_block[data$chr==i] <-
      ld$LD_block[match(data$LD_start[data$chr==i], ld$start)]
  }
  #Remove the now useless LD_start column
  data$LD_start <- NULL
  #Return the data frame with the added column of interest
  return(data)
}

#Call the function
gwas <- add_LD_blocks(gwas)

#Find number of LD blocks that have at least 1 significant SNP
result <- as.character(length(unique(gwas$LD_block)))

# Create output file and save just the name of the GWAS into this file
output_name <- paste0(sub("\\..*", "", gwas_name), "-+", result)
writeLines(as.character(gwas_name), output_name)

# Clean up
rm(gwas)
rm(ld)