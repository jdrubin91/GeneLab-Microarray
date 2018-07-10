#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")
# biocLite("arrayQualityMetrics")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list = list(
  make_option(
    c("-i", "--input"), 
    type = "character", 
    help = "Path to directory containing input raw array files"
  ),
  make_option(
    c("-n", "--normalization"),
    type = "character",
    default = "rma",
    help = "Normalization method [rma (default, full rma), quantile (no background correction), background (no quant. normalization), log2 (no quant. norm. or background correction)"
  ),
  make_option(
    c("-o", "--outFile"),
    type = "character",
    default = "expValues",
    help = "Name of the output file [without extension!] (default: expValues)"
  ),
  make_option(
    c("-i", "--ISApath"), 
    type = "character", 
    help = "Path to the directory containing the dataset metadata"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

norm = opt$normalization

addSlash = function(string) {
  # Adds a trailing forward slash to the end of a string (ex path to a driectory) if it is not present
  if (substr(x = string,
             start = nchar(string),
             stop = nchar(string)) != "/") {
    string = paste(string, "/", sep = "")
  }
  return(string)
}

if (is.null(opt$input)){ # Include GLDS accession number in outputs if provided
  print_help(opt_parser)
  stop("No path to input directory provided. Please look over the available options", call. = F)
}else{
  inPath = opt$input
  #setwd(inPath)
}

if (is.null(opt$GLDS)) {
  # Include GLDS accession number in outputs if provided
  glAn = ''
  cat("Warning: No GLDS accession number provided")
  if (grepl("GLDS-[0-9]+", inPath)) {
    glAn = regmatches(inPath, regexpr("GLDS-[0-9]+", inPath)) # Attempt to extract the GLDS accession number from the input path
  }
} else{
  glAn = paste('GLDS-', opt$GLDS, sep = '')
}

#relDir = getwd() # Save the path to the original directory (for relative path handling)

if (is.null(opt$ISApath)){
  print_help(opt_parser)
  stop("No ISA directory provided. Please re-check available options", call.=FALSE)
}

if (is.null(opt$ISApath)) {
  # Check if metadata is provided
  print_help(opt_parser)
  stop("No path to metadata directory provided. Please look over the available options",
       call. = F)
} else{
  isaPath = opt$ISApath
  isaPath = addSlash(isaPath)
  tryCatch({
    isaFiles = dir(isaPath) # Identify files in metadata directory
    sFile = isaFiles[grep("^s_*", isaFiles)] # 
    aFile = isaFiles[grep("^a_*", isaFiles)]
    sampFactors = read.delim(
      paste(isaPath, sFile, sep = ""),
      header = T,
      sep = "",
      stringsAsFactors = F
    )
    assFactors = read.delim(
      paste(isaPath, aFile, sep = ""),
      header = T,
      sep = "",
      stringsAsFactors = F
    )
  }, error = function(e) {
    stop("ISA files could not be read by parsing the tab-delimited files",
         call. = F)
  })
}



suppressPackageStartupMessages(library("limma"))

# setwd("~/Documents/genelab/rot1/GLDS-28/microarray/")
inFiles = dir(inPath)
inFiles = inFiles[grepl("_raw.txt$",inFiles)]

targets = data.frame(FileName = inFiles)
row.names(targets) = gsub(".*_microarray_","",targets$FileName)
row.names(targets) = gsub("_raw.*","",row.names(targets))
# Cy3 = need to link file names to sample names in metadata

RG = read.maimages(
  files = inFiles,
  source = "agilent.median",
  path = inPath,
  columns = list(
    G = "gMedianSignal",
    Gb = "gBGMedianSignal",
    R = "rMedianSignal",
    Rb = "rBGMedianSignal"
  ),
  annotation = "FeatureNum",
  names = row.names(targets)
)

# Potential QC
library(arrayQualityMetrics)
arrayQualityMetrics(expressionset = RG,
                    outdir = "test_report",
                    force = T)

# Normalization
RGb = backgroundCorrect(RG, method = "normexp", offset = 50) # normexp method is based on the same normal plus exponential convolution model which has previously been used to background correct Affymetrix data as part of the popular RMA algorithm [Ritchie et al. Bioinformatics 2007]
MA = normalizeWithinArrays(RG, method = "loess", weights=NULL) # Agilent specific global loess normalization method
## Loess normalization assumes that the bulk of the probes on the array are not differentially expressed
MAq = normalizeBetweenArrays(MA, method = "Aquantile") # Normalize the average intensities ("A") between arrays

# Potential QC
arrayQualityMetrics(expressionset = MAq,
                    outdir = "MA_test_report",
                    force = T)
