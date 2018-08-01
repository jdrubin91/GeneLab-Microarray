#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")
## install.packages("statmod")
# biocLite("arrayQualityMetrics")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list = list(
  make_option(
    c("-i", "--input"), 
    type = "character", 
    help = "Path to directory containing input raw array files"
  ),
  # make_option(
  #   c("-n", "--normalization"),
  #   type = "character",
  #   default = "normexp",
  #   help = "Normalization method [rma (default, full rma), quantile (no background correction), background (no quant. normalization), log2 (no quant. norm. or background correction)"
  # ),
  make_option(
    c("-o", "--outFile"),
    type = "character",
    default = "expValues",
    help = "Name of the output file [without extension!] (default: expValues)"
  ),
  make_option(
    c("-t", "--outType"),
    type = "character",
    default = "both",
    help = "Format of output data: R (Rdata object), txt (tab delimited file with identifiers and sample names), both (default)"
  ),
  make_option(
    "--QCDir",
    type = "character",
    default = "./QC_reporting/",
    help = "Path to directory for storing QC output, including a terminal forward slash. Will be created if it does not exist yet (default = './QC_reporting/')"
  ),
  make_option(
    "--GLDS", 
    type = "character", 
    help = "Full accession number for QC outputs (ie 'GLDS-21' for GLDS-21)"
  ),
  make_option(
    "--QCoutput",
    type = "logical",
    default = TRUE,
    help = "Output QC_reporting directory of QC plots (default = TRUE)"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

norm = opt$normalization
outFH = opt$outFile
QCout = opt$QCoutput


addSlash = function(string) {
  # Adds a trailing forward slash to the end of a string (ex path to a driectory) if it is not present
  if (substr(x = string,
             start = nchar(string),
             stop = nchar(string)) != "/") {
    string = paste(string, "/", sep = "")
  }
  return(string)
}

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("No path to input directory provided. Please look over the available options", call. = F)
}else{
  inPath = opt$input
  #setwd(inPath)
}

if (is.null(opt$GLDS)) {
  # Include GLDS accession number in outputs if provided
  cat("Warning: No GLDS accession number provided\n")
  if (grepl("GLDS-[0-9]+", inPath)) {
    glAn = regmatches(inPath, regexpr("GLDS-[0-9]+", inPath)) # Attempt to extract the GLDS accession number from the input path
    cat("Try to guess GLDS accession number... ", glAn,"\n")
  } else{
    glAn = FALSE
  }
} else{
  glAn = opt$GLDS
}

# Load packages
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("arrayQualityMetrics"))

# inPath = "~/Documents/genelab/rot1/GLDS-41/microarray/"
# setwd(inPath)

inFiles = dir(inPath)
inFiles = inFiles[grepl("_raw.txt$",inFiles)]

targets = data.frame(FileName = inFiles)
row.names(targets) = gsub(".*_microarray_","",targets$FileName)
row.names(targets) = gsub("_raw.*","",row.names(targets))
row.names(targets) = gsub(".txt$","",row.names(targets))

raw = read.maimages(
  files = inFiles,
  source = "agilent.median",
  green.only = T,
  path = inPath,
  columns = list(
    G = "gMedianSignal",
    Gb = "gBGMedianSignal"
  ),
  annotation = "FeatureNum",
  names = row.names(targets)
)

# Pre-normalization QC step
qcDir = addSlash(opt$QCDir)
if (!file.exists(qcDir)){ # Create QC directory if it does not exist yet
  dir.create(qcDir)
}
if(QCout == T) {
  # Load in QC package
  suppressPackageStartupMessages(require(arrayQualityMetrics))
  
  cat("Performing intial QC\n")
  
  rawData = new("ExpressionSet", exprs = as.matrix(raw)) # Generate a temprory expression set object for performing quality control
  suppressWarnings(
    arrayQualityMetrics(
      expressionset = rawData,
      outdir = paste(qcDir, "raw_report", sep = ""),
      force = T,
      do.logtransform = T
    )
  )
  rm(rawData)
}

# Normalize data
cat("Normalizing data single channel Agilent microarray data...\n")
cat("\tBackground correcting\n")

normVals = backgroundCorrect(raw, method = "normexp", verbose = F)

cat("\tNormalizing between arrays\n")

normVals = normalizeBetweenArrays(normVals, method = "quantile")

# Saving the normalized data
eset = normVals$E
row.names(eset) = normVals$genes[[1]]

outType = opt$outType
if (outType == "both") {
  save(eset, file = paste(outFH, ".rda", sep = ""))
  write.table(
    data.frame("ID" = row.names(eset),eset), # provides the rownames as a labeled column in the saved output
    row.names = F,
    file = paste(outFH, ".txt", sep = ""),
    sep = "\t",
    quote = F
  )
  cat("Success! Normalized data saved to", paste(outFH, sep=""), "as both a .txt and a .RData file\n\n")
} else if (outType == "R") {
  save(eset, file = paste(outFH, ".rda", sep = ""))
  cat("Success! Normalized data saved to", paste(outDir, outFH, sep=""), "as a .RData file\n\n")
} else if (outType == "txt") {
  write.table(
    data.frame("ID" = row.names(eset),eset), # provides the rownames as a labeled column in the saved output
    row.names = F,
    file = paste(outFH, ".txt", sep = ""),
    sep = "\t",
    quote = F
  )
  cat("Success! Normalized data saved to", paste(outDir, outFH, sep=""), "as a .txt file\n\n")
} else{
  print_help(opt_parser)
  stop("Help, I don't know how to save this data!", call. = F)
}

# Post-normalization QC step
if(QCout == T) {
  cat("Performing post-normalization QC\n")
  
  normalizedData = new("ExpressionSet", exprs = as.matrix(eset)) # Generate a temprory expression set object for performing quality control
  suppressWarnings(
    arrayQualityMetrics(
      expressionset = normalizedData,
      outdir = paste(qcDir, "normalized_report", sep = ""),
      force = T,
      do.logtransform = T
    )
  )
  rm(normalizedData)
}


















