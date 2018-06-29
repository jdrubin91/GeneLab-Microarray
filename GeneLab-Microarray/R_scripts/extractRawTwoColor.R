#!/usr/bin/env Rscript

# install.packages("optparse")
suppressPackageStartupMessages(library("optparse"))

# Read options
option_list = list(
  make_option(c("-i", "--input"), type = "character", help = "Path to file with buried raw array information [required]"),
  make_option(c("-o", "--output"), type = "character", help = "Path and filename for saving extracted data [required]")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

norm = opt$normalization

if (is.null(opt$input)) {
  # Check and set input file handle
  print_help(opt_parser)
  stop("No path to input file provided. Please look over the available options\n",
       call. = F)
} else
  inFH = opt$input

if (is.null(opt$output)) {
  # Check and set output file handle
  print_help(opt_parser)
  stop("No path to out file provided. Please look over the available options\n",
       call. = F)
} else
  outFH = opt$output

cat("\nReading in file:", inFH, "\n")
test = read.delim(inFH, stringsAsFactors = F, header = F) # Read in processed file

cat("Extracting raw values...\n")
cat("|     |\r")
G = test[, grep("gMedianSignal", test)] # Idenitfy median foreground intensity columns
cat("|-      |\r")
R = test[, grep("rMedianSignal", test)]
cat("|--     |\r")
Gb = test[, grep("gBGMedianSignal", test)] # Idenitfy median background intensities columns
cat("|---   |\r")
Rb = test[, grep("rBGMedianSignal", test)]
cat("|---- |\r")
startInd = grep("gMedianSignal", G) + 1 # Indentify starting row index of raw values
G = G[startInd:length(G)] # Extract raw values
R = R[startInd:length(R)]
Gb = Gb[startInd:length(Gb)]
Rb = Rb[startInd:length(Rb)]
raw = cbind(R, G, Rb, Gb) # Bind raw value vectors into dataframe
row.names(raw) = test[startInd:nrow(test), grep("FeatureNum", test)] # Add unique feature numbers for mapping to genes
cat("|-----|\r\n")

write.table(raw, file = outFH, sep = "\t", quote = F)

cat("Success! Raw values saved to:", outFH, "\n\n")
