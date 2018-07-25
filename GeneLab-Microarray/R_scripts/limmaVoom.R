#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list = list(
  make_option(
    c("-m", "--countsMatrix"), 
    type = "character", 
    help = "Path to a tab-delimited counts matrix .txt file with the sample names in the column headers"
  ),
  make_option(
    c("-f", "--countsFiles"), 
    type = "character", 
    help = "A list of tab-delimited counts .txt files to read in and concatenate into a dataframe within R"
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

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

addSlash = function(string) {
  # Adds a trailing forward slash to the end of a string (ex path to a driectory) if it is not present
  if (substr(x = string, start = nchar(string), stop = nchar(string)) != "/") {
    string = paste(string, "/", sep = "")
  }
  return(string)
}

# wd = "/Users/dmattox/Documents/genelab/RNAseq/GLDS-101/RNAseq/Feature counts/"
# inFH = dir(wd)

if (is.null(opt$countsMatrix) & is.null(opt$countsFiles)) {
  print_help(opt_parser)
  stop("No counts data provided", call. = FALSE)
} else {
  if (!is.null(opt$countsMatrix)) {
    inFH = opt$countsMatrix
  } else {
    inFH = opt$countsFiles
  }
}

if (is.null(opt$ISApath)) {
  print_help(opt_parser)
  stop("No ISA file provided", call. = FALSE)
} else {
  isaFH = opt$ISApath
}

# Read in underscore-delimited factor levels and split into lists
if (!is.null(opt$group1) & !is.null(opt$group2)){
  fact1 = opt$group1
  fact1 = strsplit(fact1, split = "_")[[1]]
  fact2 = opt$group2
  fact2 = strsplit(fact2, split = "_")[[1]]
} else{
  print_help(opt_parser)
  stop("Factor levels not provided or improperly formated", call.=FALSE)
}

suppressPackageStartupMessages(library("limma"))

# Read in ISA tab file and extract sample file
# isaFH = "../../GLDS-101_metadata_ISA/s_BIOBANK.txt"
tryCatch({
  studyFactors = read.delim(isaFH, header=T, sep="",stringsAsFactors = F)
}, error=function(e){ 
  stop("ISA sample file could not be read by parsing tab-delimited files", call. = F)
})

# From sample file, extract column containing 'Factor Value'
tryCatch({
  factorValues = as.data.frame(studyFactors[, grepl("Factor.Value", colnames(studyFactors))])
  row.names(factorValues) = studyFactors[, grepl("^Sample.Name$", colnames(studyFactors))]
  {
    # Match sample names in file name
    ## [replace('_','-').replace('(','-').replace(')','-').replace(' ','-').strip('-')]
    replaceWithHyphen = c("_", "\\(", "\\)", " ", "\\.")
    removeList = c("^-", "-$")
    for (i in 1:length(replaceWithHyphen)) {
      # Replace other characters with hyphens
      row.names(factorValues) =  gsub(replaceWithHyphen[i], '-', row.names(factorValues))
    }
    for (i in 1:length(removeList)) {
      # Remove leading/trailing hypens
      row.names(factorValues) =  gsub(removeList[i], '', row.names(factorValues))
    }
    
  }
}, error = function(e) {
  stop("Unable to pull factors or sample names from the study level metadata",
       call. = F)
})






