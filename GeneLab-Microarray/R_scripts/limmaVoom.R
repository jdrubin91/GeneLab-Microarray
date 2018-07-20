#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list = list(
  make_option(
    c("-d", "--exprData"), 
    type = "character", 
    help = "Name of (or path to) the input file (tab delimited .txt file or binary RData object)"
  ),
  make_option(
    c("-i", "--ISApath"), 
    type = "character", 
    help = "Path to the file containing the sample-level metadata"
  ),
  make_option(
    "--group1", 
    type = "character", 
    help = "'_'delimited list of factors to select samples for group 1 [ex: flight_geneKO]"
  ),
  make_option(
    "--group2", 
    type = "character", 
    help = "'_'delimited list of factors to select samples for group 2 [ex: ground_geneKO]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "DGE.txt",
    help = "Name of (or path to) file to write results to (default: DGE.txt)"
  ),
  make_option(
    c("-r", "--rmOutliers"), 
    type = "character", 
    help = "Underscore-delimited list of samples to exclude as outliers from differential expression analysis, matching the sample names in the metadata [ex: GSM1234_GSM1235]"
  )
)