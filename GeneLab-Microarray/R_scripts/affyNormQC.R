#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("affy")
# biocLite("affyPLM")
# biocLite("oligo")
# biocLite("arrayQualityMetrics")
#   biocLite("hexbin")
#   biocLite("jsonlite")
#   biocLite("openssl")
#   biocLite("stringi")
#   biocLite("reshape2")
#   biocLite("Cairo")

#   biocLite("")
#   biocLite("")
#   biocLite("")
#   biocLite("")
#   biocLite("")
#   biocLite("")


suppressPackageStartupMessages(library("optparse"))

relDir = getwd()

# Read options
option_list = list(
  make_option(
    c("-i", "--input"), 
    type = "character", 
    help = "Path to directory containing input .CEL files"
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
    "--outDir",
    type = "character",
    default = "./",
    help = "Path to an output directory, including a terminal forward slash (default: directory this script is called from)"
  ),
  make_option(
    c("-t", "--outType"),
    type = "character",
    default = "both",
    help = "Format of output data: R (Rdata object), txt (tab delimited file with identifiers and sample names), both (default)"
  ),
  make_option(
    "--outputData",
    type = "logical",
    default = TRUE,
    help = "Output data at all (default TRUE)"
  ),
  make_option(
    c("-a", "--arrayInfoOnly"),
    type = "logical",
    default = FALSE,
    help = "Detect-affy-array-only mode. If true, script will exit after outputting the arrayInfo file. (Default: FALSE)"
  ),
  make_option(
    "--QCoutput",
    type = "logical",
    default = TRUE,
    help = "Output QC_reporting directory of QC plots (default = TRUE)"
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

detach_package = function(pkg, character.only = FALSE) {
  if (!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while (search_item %in% search())
  {
    detach(search_item,
           unload = TRUE,
           character.only = TRUE)
  }
}

norm = opt$normalization
QCout = opt$QCoutput
NUSEplot = opt$NUSEplot

if (is.null(opt$input)) {
  # Check for provided input directory
  print_help(opt_parser)
  stop("No path to input directory provided. Please look over the available options",
       call. = F)
} else{
  inPath = addSlash(opt$input)
  setwd(inPath) # Change the working directory to the directory containing the raw files
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

# Load affy package to read in .CEL files
suppressPackageStartupMessages(require(affy))

# setwd("~/Documents/genelab/rot1/GLDS-4/microarray/")

celFiles = list.celfiles(full.names = TRUE)
sampNames = gsub("_microarray_.*", "", celFiles)
sampNames = gsub(".CEL", "", sampNames)
sampNames = gsub(".*/", "", sampNames)
sampNames = gsub("GLDS-\\d*_", "", sampNames)# Extract sample names from the list of .CEL files

if (length(celFiles) > 0){
  cat("Detected .CEL files:\n")
  for (i in 1:length(celFiles)) {
    cat("\t",celFiles[i],"\n")
  }
  cat("\n")
} else {
  stop("No .CEL files detected in the current directory",
    call. = F)
}

tryCatch({
  useOligo = tryCatch({
    suppressWarnings(expr = {
      raw = ReadAffy(filenames = celFiles,
                     sampleNames = sampNames)
    })
    
    arrInfo = c("Affymetrix", as.character(raw@cdfName))
    
    if (grepl("-st-", raw@cdfName, ignore.case = T)) {
      detach_package(affy)
      rm(raw)
      suppressPackageStartupMessages(require(oligo))
      raw = read.celfiles(filenames = celFiles,
                          sampleNames = sampNames)
      st = T
    } else{
      suppressPackageStartupMessages(require(affyPLM))
      st = F
    }
    return(FALSE)
  }, error = function(e) {
    cat(
      "Could not read in provided .CEL files with affy package. Attempting to read them with oligo package...\n"
    )
    return(TRUE)
  })
}, error = function(e) {
  stop("Unable to read in .CEL files with affy package", call. = F)
})

if (useOligo == TRUE) {
  tryCatch({
    detach_package(affy)
    suppressPackageStartupMessages(require(oligo))
    raw = read.celfiles(filenames = celFiles,
                        sampleNames = sampNames)
    st = T
    ver = raw@annotation
    if (grepl("^pd.", ver)) {
      # Convert from the pd.* annotation package to the standard array version name
      ver = gsub("^pd.", "", ver)
      ver = gsub("\\.", "-", ver)
      ver = gsub("(\\d)(-)(\\d)", "\\1_\\3", ver)
    }
    arrInfo = c("Affymetrix", as.character(ver))
  }, error = function(e) {
    stop("Unable to read in .CEL files with affy package", call. = F)
  })
}

setwd(relDir) # Return the working directory to direcotry script was called from to enable use of relative paths
# Create QC output directory
qcDir = addSlash(opt$QCDir)
if (!file.exists(qcDir)){ # Create QC directory if it does not exist yet
  dir.create(qcDir)
}

# Output array information to a separate file
summDir = paste(qcDir, "summary_report/", sep = "")
if (!file.exists(summDir)){ # Create a summary report directory within qcDir if it does not exist yet
  dir.create(summDir)
}

cat("\n\nTROUBLESHOOTING:\n",arrInfo,"\n",paste(summDir, "arrayInfo.txt", sep = ""),"\n")

if (glAn != FALSE) {
  write.table(
    arrInfo,
    file = paste(summDir, glAn, "_arrayInfo.txt", sep = ""),
    quote = F,
    col.names = F,
    row.names = F
  )
} else {
  write.table(
    arrInfo,
    file = paste(summDir, "arrayInfo.txt", sep = ""),
    quote = F,
    col.names = F,
    row.names = F
  )
}

cat("Array type detected.\n")

# Exit script if arrayInfoOnly mode is True
if (opt$arrayInfoOnly == TRUE) {
  cat("Detect-array-type-only mode on, exiting.\n")
  quit(save = "no",
       status = 0,
       runLast = FALSE)
}

## Raw QC
if(QCout == T) {

  # Load in QC package
  suppressPackageStartupMessages(require(arrayQualityMetrics))

  cat("Performing intial QC\n")

  suppressWarnings(
    arrayQualityMetrics(
      expressionset = raw,
      outdir = paste(qcDir, "raw_report", sep = ""),
      force = T,
      do.logtransform = T
    )
  )

}

outFH = opt$outFile
if (opt$outputData == TRUE) {
  ## Normalize
  cat("\nNormalizing with selected normalization technique...\n")
  if (norm == 'rma') {
    expset = rma(raw)
  } else if (norm == 'quantile') {
    expset = rma(raw, background = F, normalize = T)
  } else if (norm == 'background') {
    expset = rma(raw, background = T, normalize = F)
  } else if (norm == 'log2') {
    expset = rma(raw, background = F, normalize = F)
  } else{
    stop(
      "Normalization did not occur, please examine script inputs and default values\n",
      call. = F
    )
  }
  
  outDir = addSlash(opt$outDir)
  if (!file.exists(outDir)){ # Create the output directory if it does not exist yet
    dir.create(outDir)
  }
  eset = exprs(expset)
  if (opt$outType == "both") {
    save(eset, file = paste(outDir, outFH, ".rda", sep = ""))
    write.table(
      data.frame("ID" = row.names(eset),eset), # provides the rownames as a labeled column in the saved output
      row.names = F,
      file = paste(outDir, outFH, ".txt", sep = ""),
      sep = "\t",
      quote = F
    )
    cat("Success! Normalized data saved to", paste(outDir, outFH, sep=""), "as both a .txt and a .RData file\n\n")
  } else if (opt$outType == "R") {
    save(eset, file = paste(outDir, outFH, ".rda", sep = ""))
    cat("Success! Normalized data saved to", paste(outDir, outFH, sep=""), "as a .RData file\n\n")
  } else if (opt$outType == "txt") {
    write.table(
      data.frame("ID" = row.names(eset),eset), # provides the rownames as a labeled column in the saved output
      row.names = F,
      file = paste(outDir, outFH, ".txt", sep = ""),
      sep = "\t",
      quote = F
    )
    cat("Success! Normalized data saved to", paste(outDir, outFH, sep=""), "as a .txt file\n\n")
  } else{
    print_help(opt_parser)
    stop("Help, I don't know how to save this data!", call. = F)
  }
  
  if(QCout == T) {
    cat("Post normalization QC steps...\n")
    # Post-normalization QC

    suppressWarnings(
      arrayQualityMetrics(
        expressionset = expset,
        outdir = paste(qcDir, "normalized_report", sep = ""),
        force = T
      )
    )
  }
}
