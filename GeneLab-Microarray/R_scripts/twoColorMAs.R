#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list=list(
  make_option(c("-i","--input"),type="character",help="Path to directory containing input .CEL files"),
  make_option(c("-n","--normalization"),type="character",default="rma",help="Normalization method [rma (default, full rma), quantile (no background correction), background (no quant. normalization), log2 (no quant. norm. or background correction)"),
  make_option(c("-o","--outFile"),type="character",default="expValues",help="Name of the output file [without extension!] (default: expValues)"),
  make_option(c("-t","--outType"),type="character",default="txt",help="Format of output data: R (Rdata eset object), txt (default, tab delimited file with identifiers and sample names), both"),
  make_option("--outputData",type="logical",default=TRUE,help="Output data at all (default TRUE)"),
  make_option(c("-a","--arrayInfoOnly"),type="logical",default=FALSE,help="Detect-affy-array-only mode. If true, script will exit after outputting the arrayInfo file. (Default: FALSE)"),
  make_option("--QCoutput",type="logical",default=TRUE,help="Output QC_reporting directory of QC plots (default = TRUE)"),
  make_option("--NUSEplot",type="logical",default=FALSE,help="Include a NUSE plot in the QC output, adds significantly to runtime (default = FALSE)"),
  make_option("--GLDS",type="character",help="GLDS accession number for plot outputs (ie '21' for GLDS-21)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

norm = opt$normalization

if (is.null(opt$GLDS)){ # Include GLDS accession number in outputs if provided
  print_help(opt_parser)
  glAn = ''
  cat("Warning: No GLDS accession number provided")
}else{
  glAn = paste('GLDS-',opt$GLDS,sep='')
}

if (is.null(opt$input)){ # Include GLDS accession number in outputs if provided
  print_help(opt_parser)
  stop("No path to input directory provided. Please look over the available options", call. = F)
}else{
  inPath = opt$input
  setwd(inPath)
}
