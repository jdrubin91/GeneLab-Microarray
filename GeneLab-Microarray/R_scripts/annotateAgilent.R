#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list = list(
  make_option(
    c("-i", "--input"), 
    type = "character", 
    help = "Name of (or path to) the input file (tab delimited .txt file or binary Rdata object)"
  ),
  make_option(
    c("-g", "--gpl"), 
    type = "character", 
    default = "search",
    help = "Set to 'search' to automatically look for a GPL file in the same directory as the input, otherwise provide a path to the file containing custom array annotation information, containing column names but no additional header"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "annotExpValues",
    help = "Name of (or path to) file to write results to, WITHOUT file extension (default: annotExpValues)"
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
    help = "Path to directory to store the annotation output. Will be created if it does not exist yet, but it is recommended to use the same directory as was used for QC with the normalization step (default = './QC_reporting/')"
  ),
  make_option(
    c("-d", "--dupProbes"),
    type = "character",
    default = "max",
    help = "Method for handling multiple probes [max (default, probe with the highest mean expression), average (mean of all probes for a gene), topvar (highest variance with nsFilter function)"
  ),
  make_option(
    "--GLDS", 
    type = "character", 
    help = "Full accession number for QC outputs (ie 'GLDS-21' for GLDS-21)"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

addSlash = function(string) {
  # Adds a trailing forward slash to the end of a string (ex path to a driectory) if it is not present
  if (substr(x = string,
             start = nchar(string),
             stop = nchar(string)) != "/") {
    string = paste(string, "/", sep = "")
  }
  return(string)
}

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call. = FALSE)
}else { inFH = opt$input }

norm = opt$normalization
outFH = opt$outFile
if (!is.null(opt$gpl)){
  if (opt$gpl == "search") {
    cat("Looking in the 'input' directory for a GPL file...\n")
    inDir = gsub("(/)((\\w)*\\.(\\w)*)$", "\\1", inFH) # Strip the filename away from the directory path to the input file
    files = dir(inDir)
    files = files[grepl("^GPL[[:digit:]]*", files)]
    if (length(files) == 1) {
      annotFH = files[1]
      cat("\t",annotFH,"identified as an annotation file\n")
    }
  }
  annotFH = opt$gpl # annotFH = "GPL10094-20413.txt"
}

# Read in annotation file
tryCatch({
  annot = read.delim(annotFH,header = T, stringsAsFactors = F)
}, error = function(e) {
  stop(paste("Unable to read in ",annotFH, sep = ""),
       call. = F)
})









