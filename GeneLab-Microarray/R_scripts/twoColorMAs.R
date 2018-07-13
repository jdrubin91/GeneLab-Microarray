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
    default = "normexp",
    help = "Normalization method [rma (default, full rma), quantile (no background correction), background (no quant. normalization), log2 (no quant. norm. or background correction)"
  ),
  make_option(
    c("-o", "--outFile"),
    type = "character",
    default = "expValues",
    help = "Name of the output file [without extension!] (default: not_sure_yet)"
  ),
  make_option(
    c("-i", "--ISApath"), 
    type = "character", 
    help = "Path to the directory containing the dataset metadata"
  ),
  make_option(
    c("-g", "--gpl"), 
    type = "character", 
    help = "Path to the file containing custom array annotation information"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

norm = opt$normalization
outFH = opt$outFile
if (!is.null(opt$gpl)){
  annotFH = opt$gpl # annotFH = "GPL15420_features_probes.txt"
}

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

# Link samples to factors and labels
factors = merge(x = sampFactors, y = assFactors, by = "Sample.Name") # Join the information from the ISA tab delimited assay and sample files

suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("arrayQualityMetrics"))

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
arrayQualityMetrics(expressionset = RG,
                    outdir = "test_report",
                    force = T)

# Normalization
RGb = backgroundCorrect(RG, method = "normexp", offset = 50) # normexp method is based on the same normal plus exponential convolution model which has previously been used to background correct Affymetrix data as part of the popular RMA algorithm [Ritchie et al. Bioinformatics 2007]
MA = normalizeWithinArrays(RG, method = "loess", weights=NULL) # Agilent specific global loess normalization method
## Loess normalization assumes that the bulk of the probes on the array are not differentially expressed
MAq = normalizeBetweenArrays(MA, method = "Aquantile") # Normalize the average intensities ("A") between arrays

# Feature number to BSU locus
annotFun = function(oldLab, newLink, newLab, ...){
  # Function to switch from one probe label (oldLab) to another (newLab), provided
  ## a linking variable (newLink) with common labels to the labels in the oldLab list organized in the same order as the newLab list
  # Returns a translation from the the oldLab to the newLab with "" for unmapped probes
  newOrder = match(oldLab,newLink)
  out = newLab[newOrder]
  return(out)
}

annot = read.delim(annotFH,header = T, stringsAsFactors = F)
BSUs = annotFun(
  oldLab = MAq$genes[,1], 
  newLink = annot$FeatureNum,
  newLab = annot$BSU
  )
noIDCnt = sum(BSUs == "") # Number of unmapped probes
MAq$genes[,1] = BSUs # Translate probe IDs
MAq = MAq[!MAq$genes[,1] == "", ] # Remove umapped probes
BSUs = BSUs[!BSUs == ""]

if (any(opt$dupProbes %in% c("average", "max"))) {
  cat("Filtering out multiple probes per gene ID...\n")
  
  if (opt$dupProbes == "average") {
    # Collapse multiple probes per gene ID by averaging expression values across all samples
    rmRowTag = rep(TRUE, length(MAq$genes[,1])) # Tag rows to drop (set single or averaged probes to FALSE below)
    for (i in 1:length(rmRowTag)) {
      if (sum(BSUs == BSUs[i]) > 1) {
        inds = grep(BSUs[i], BSUs) # List of indices at which a probe for a given gene ID occur
        MAq$M[inds[1], ] = apply(X = MAq$M[inds, ],
                                FUN = mean,
                                MARGIN = 2) # Changes the values of the first occurence of a probe to the sample-specific average of the values from all the probes for that gene ID
        MAq$A[inds[1], ] = apply(X = MAq$A[inds, ],
                                 FUN = mean,
                                 MARGIN = 2) # Changes the values of the first occurence of a probe to the sample-specific average of the values from all the probes for that gene ID
        rmRowTag[inds[1]] = FALSE
      } else
        rmRowTag[i] = FALSE
    }
    nDups = sum(rmRowTag)
    MAq = MAq[!rmRowTag, ]

    cat("\tUnampped probes removed:", noIDCnt, "\n")
    cat("\tDuplicated probes removed:", nDups, "\n\n")
    cat("Annotated probes remaining:", length(MAq$genes[,1]), "\n\n")
    if (length(MAq$genes[,1]) > length(unique(MAq$genes[,1]))) {
      cat("\n\tWarning: non-unique probe to ID mappings remain \n")
    }
    
  } else if (opt$dupProbes == "max") {
    # Collapse multiple probes per gene ID by selecting a representative with the highest mean intensity across all samples
    rmRowTag = rep(TRUE, length(MAq$genes[,1])) # Tag rows to drop (set single or highest expressing probes to FALSE below)
    for (i in 1:length(MAq$genes[,1])) {
      if (sum(BSUs == BSUs[i]) > 1) {
        inds = grep(BSUs[i], BSUs)
        top = 0
        keep = 0
        for (j in 1:length(inds)) {
          curr = mean(as.numeric(MAq$A[inds[j], ]))
          if (is.na(curr)){ # Check if at least one of the samples has a missing value for the current probe
            curr = 0
          }
          if (curr > top) {
            top = curr
            keep = inds[j]
          }
        }
        rmRowTag[keep] = FALSE
      } else
        rmRowTag[i] = FALSE
    }
    nDups = sum(rmRowTag)
    MAq = MAq[!rmRowTag, ]

    cat("\tUnampped probes removed:", noIDCnt, "\n")
    cat("\tDuplicated probes removed:", nDups, "\n\n")
    cat("Annotated probes remaining:", length(MAq$genes[,1]), "\n\n")
    if (length(MAq$genes[,1]) > length(unique(MAq$genes[,1]))) {
      cat("\n\tWarning: non-unique probe to ID mappings remain \n")
    }
  }
} else{
  stop("Method for dealing with probes mapped to the same gene IDs not recognized\n",
       call. = F)
}



# write.table(
#   MAq,
#   file = paste(outFH,".txt", sep = ""),
#   sep = "\t",
#   quote = F,
#   row.names = FALSE
# )

# Potential QC
suppressWarnings(
  arrayQualityMetrics(
    expressionset = MAq,
    outdir = "MA_test_report",
    force = T
  )
)
