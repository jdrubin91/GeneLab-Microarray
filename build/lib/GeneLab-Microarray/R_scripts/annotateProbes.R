#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")
# biocLite("mogene10sttranscriptcluster.db")
# biocLite("moe430a.db")
# biocLite("drosophila2.db")
# biocLite("hgu133plus2.db")
# biocLite("ath1121501.db")
# biocLite("yeast2.db")
# biocLite("hugene10sttranscriptcluster.db")
# biocLite("rat2302.db")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list = list(
  make_option(
    c("-i", "--input"), 
    type = "character", 
    help = "Name of (or path to) the input file (tab delimited .txt file or binary Rdata object)"
    ),
  make_option(
    c("-a", "--arrayInfo"),
    type = "character",
    default = "./QC_output/summary_report/arrayInfo.txt",
    help = "Name of (or path to) a file containing the array information [Line 1: Manufacturer, line 2: Array version]"
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
    help = "Path to directory to save the annotation output. Will be created if it does not exist yet, but it is recommended to use the same directory as was used for QC with the normalization step (default = './QC_reporting/')"
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

if (is.null(opt$arrayInfo)) {
  print_help(opt_parser)
  stop("Array info file not provided", call. = FALSE)
}else { 
  aiFH = opt$arrayInfo
  # aiFH = "arrayInfo.txt"
}

if (is.null(opt$GLDS)) {
  # Include GLDS accession number in outputs if provided
  glAn = ''
  cat("Warning: No GLDS accession number provided\n")
  if (grepl("GLDS-[0-9]+", aiFH)) {
    glAn = regmatches(aiFH, regexpr("GLDS-[0-9]+", aiFH)) # Attempt to extract the GLDS accession number from the input path
    cat("Try to guess GLDS accession number... ", glAn,"\n")
  } else{
    glAn = FALSE
  }
} else{
  glAn = opt$GLDS
}

# Create QC output directory
qcDir = addSlash(opt$QCDir)
if (!file.exists(qcDir)){ # Create QC directory if it does not exist yet
  dir.create(qcDir)
}

tryCatch({
  aInf = read.delim(aiFH, header = F, stringsAsFactors = F)
  arrMan = aInf[1, 1]
  arrVer = aInf[2, 1]
}, error = function(e) {
  stop("Array info file not found or organization not recognized, check options help",
       call. = F)
})

if (arrVer == "Pae_G1a") {
  stop(paste("Array version", arrVer, "does not have an annotation package and is not able to be annotated with this script. Try using the annotateAgilent.R script and the appropriate GPL file"))
}

# Set-up array version:annotation database pseudo-dictionary
arrayNames = c(
  "MoGene-1_0-st-v1",
  "MOE430A",
  "Drosophila_2",
  "HG-U133_Plus_2",
  "ATH1-121501",
  "HuGene-1_0-st-v1",
  "Yeast_2",
  "Rat230_2",
  "hta-2_0"
)

arrPackages = c(
  "mogene10sttranscriptcluster.db",
  "moe430a.db",
  "drosophila2.db",
  "hgu133plus2.db",
  "ath1121501.db",
  "hugene10sttranscriptcluster.db",
  "yeast2.db",
  "rat2302.db",
  "hta20transcriptcluster.db"
)

# Call the appropriate annotation package
tryCatch({
  annotPack = arrPackages[grep(pattern = arrVer,
                               x = arrayNames,
                               ignore.case = T)] # Pick out appropriate package by the array version
  if (length(annotPack) != 0) {
    tryCatch({
      suppressPackageStartupMessages(library(annotPack, character.only = T)) # Load selected package
    }
    , error = function(e) {
      cat(
        "Package recognized but was not found installed in this environment. Attempting to install now\n"
      )
      tryCatch({
        source("http://bioconductor.org/biocLite.R")
        biocLite(annotPack)
      }, error = function(e) {
        cat(
          "Package failed to install. Consider manually installing the annotation packages listed at the top of the script\n"
        )
      })
      tryCatch({
        suppressPackageStartupMessages(library(annotPack, character.only = T)) # Load selected package after installing
      }, error = function(e) {
        cat(
          "Package installed but was unable to load\n"
        )
      })
    })
  } else{
    if (grepl("-st(-)*", arrVer, ignore.case = T)) {
      affyST = TRUE
    } else{
      affyST = FALSE
    }
    annotPack = tolower(arrVer)
    annotPack = gsub("[[:punct:]]", "", annotPack)
    if (affyST == TRUE) {
      annotPack = gsub("v[0-9]*$", "", annotPack)
      annotPack = paste(annotPack, "transcriptcluster.db", sep = "")
    } else{
      annotPack = paste(annotPack, ".db", sep = "")
    }
    cat(
      "For Affymetrix array type:",
      arrVer,
      "guessing the annotation package:",
      annotPack,
      "\n\tAttempting to load package now...\n"
    )
    tryCatch({
      suppressPackageStartupMessages(library(annotPack, character.only = T)) # Load selected package
    }
    , error = function(e) {
      cat(
        "Package recognized but was not found installed in this environment. Attempting to install now\n"
      )
      tryCatch({
        source("http://bioconductor.org/biocLite.R")
        biocLite(annotPack)
      }, error = function(e) {
        cat(
          "Package failed to install. Consider manually installing the appropriate annotation package and adding it to the list of encountered packages above.\n"
        )
      })
      tryCatch({
        suppressPackageStartupMessages(library(annotPack, character.only = T)) # Load selected package after installing
      }, error = function(e) {
        cat(
          "Package installed but was unable to load\n"
        )
      })
    })
  }
  
  packObjs = ls(paste("package:", as.character(annotPack), sep = "")) # Stores a list of all the objects in the selected package
  if (any(grepl(
    pattern = "REFSEQ",
    x = packObjs,
    ignore.case = T
  ))) {
    annotEnv = packObjs[grepl(pattern = "REFSEQ",
                              x = packObjs,
                              ignore.case = T)] # Select the enivornment from the package to map probes to RefSeq IDs
  } else if (annotPack == "ath1121501.db") {
    annotEnv = packObjs[grepl(pattern = "ACCNUM",
                              x = packObjs,
                              ignore.case = T)] # Select the enivornment from the package to map probes to RefSeq IDs
  } else if (annotPack == "yeast2.db") {
    annotEnv = packObjs[grepl(pattern = "ORF",
                              x = packObjs,
                              ignore.case = T)] # Select the enivornment from the package to map probes to RefSeq IDs
  }
  AR1 = paste("This microarray experiment was annotated with the R package: ", annotPack,", using the object: ",annotEnv,".",sep="")
  cat("Annotating with R package",
      annotPack,
      "using object:",
      annotEnv,
      "\n")
}, error = function(e) {
  stop(
    "Array version wasn't not recognized or the annotation package was unable to load.\nCheck that the appropriate packages are installed and the array version is contained in the list of known arrays\n",
    call. = F
  )
})

# For-loop to check that guess-regex matches all known array version-annotation package pairs
# for(i in 1:length(arrayNames)){
#   if(grepl("-st-",arrayNames[i],ignore.case = T)){affyST = TRUE} else{affyST = FALSE}
#   pack = tolower(arrayNames[i])
#   pack =  gsub("[[:punct:]]","",pack)
#   if(affyST == TRUE){
#     pack=gsub("v[0-9]*$","transcriptcluster.db",pack)
#   }else{
#     pack = paste(pack,".db",sep="")
#   }
#   cat(pack,"\t\t\t\t",arrPackages[i],"\n",pack == arrPackages[i],"\n\n")
# }

# inFH = "expValues.txt"
tryCatch({
  if (grepl(".txt$", x = inFH) == TRUE) {
    eset = read.delim(
      inFH,
      header = T,
      sep = "\t",
      stringsAsFactors = F
    )
    rownames(eset) = eset[,1]
    eset[,1] = NULL
  } else{
    load(inFH)
  }
  neset = new("ExpressionSet", exprs = as.matrix(eset))
  neset@annotation = annotPack
}, error = function(e) {
  stop("Input file was not recognized", call. = F)
})

# Mapping probe IDs to RefSeq names from the imported library
mapFun = function(id, environ) {
  # Function to match the primary RefSeq ID for a given probe ID and return NA in all other cases
  return(tryCatch(
    get(id, env = environ)[1],
    error = function(e)
      NA
  ))
}

if(opt$dupProbes == "topvar") {
  # Collapse multiple probes per gene ID by selecting a representative with the most variance in expression across all samples
  suppressPackageStartupMessages(library("genefilter"))
  cat("Filtering out unannotated probes...\n")
  filt = nsFilter(
    neset,
    var.filter = F,
    require.entrez = T,
    remove.dupEntrez = T
  )
  nDups = filt[[2]]$numDupsRemoved # Number of probes removed that map to non-unique gene IDs
  filtID = featureNames(filt[[1]]) # Pulls out the probe IDs
  cat("Mapping probes IDs to gene IDs...\n")
  filtRefSeq = lapply(filtID, FUN = mapFun, environ = eval(parse(text =
                                                                   annotEnv))) # Applying mapFun to all non-filtered probe IDs
  noIDCnt = nrow(eset) - sum(!is.na(filtRefSeq)) - nDups
  
  if (sum(!is.na(filtRefSeq)) > length(unique(filtRefSeq[!is.na(filtRefSeq)]))) {
    cat("\n\tWarning: non-unique probe to ID mappings remain \n")
  }
  
  # Replace AffyIDs with RefSeq IDs, drop probes w/o RefSeq IDs
  normVals = exprs(filt[[1]])
  normVals = normVals[!is.na(filtRefSeq), ]
  rownames(normVals) = filtRefSeq[!is.na(filtRefSeq)]
  
  
} else if (any(opt$dupProbes %in% c("average", "max"))) {
  cat("Mapping probes IDs to gene IDs...\n")
  RefSeq = lapply(rownames(eset), FUN = mapFun, environ = eval(parse(text =
                                                                       annotEnv))) # Applying mapFun to all probe IDs
  noIDCnt = sum(is.na(RefSeq)) # Count unmapped probes
  eset = eset[!is.na(RefSeq), ] # Remove data from unmapped probes
  RefSeq = RefSeq[!is.na(RefSeq)] # Remove NAs so gene IDs correspond to eset rows
  cat("Filtering out unannotated probes...\n")
  
  if (opt$dupProbes == "average") {
    # Collapse multiple probes per gene ID by averaging expression values across all samples
    rmRowTag = rep(TRUE, nrow(eset)) # Tag rows to drop (set single or averaged probes to FALSE below)
    for (i in 1:nrow(eset)) {
      if (sum(RefSeq == RefSeq[i][[1]]) > 1) {
        inds = grep(RefSeq[i][[1]], RefSeq) # List of indices at which a probe for a given gene ID occur
        eset[inds[1], ] = apply(X = eset[inds, ],
                                FUN = mean,
                                MARGIN = 2) # Changes the values of the first occurence of a probe to the sample-specific average of the values from all the probes for that gene ID
        rmRowTag[inds[1]] = FALSE
      } else
        rmRowTag[i] = FALSE
    }
    nDups = sum(rmRowTag)
    normVals = eset[!rmRowTag, ]
    row.names(normVals) = RefSeq[!rmRowTag]
    
    if (nrow(normVals) > length(unique(RefSeq[!rmRowTag]))) {
      cat("\n\tWarning: non-unique probe to ID mappings remain \n")
    }
    
  } else if (opt$dupProbes == "max") {
    # Collapse multiple probes per gene ID by selecting a representative with the highest mean expression across all samples
    rmRowTag = rep(TRUE, nrow(eset)) # Tag rows to drop (set single or highest expressing probes to FALSE below)
    for (i in 1:nrow(eset)) {
      if (sum(RefSeq == RefSeq[i][[1]]) > 1) {
        inds = grep(RefSeq[i][[1]], RefSeq)
        top = 0
        keep = 0
        for (j in 1:length(inds)) {
          curr = mean(as.numeric(eset[inds[j], ]))
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
    normVals = eset[!rmRowTag, ]
    row.names(normVals) = RefSeq[!rmRowTag]
    
    if (nrow(normVals) > length(unique(RefSeq[!rmRowTag]))) {
      cat("\n\tWarning: non-unique probe to ID mappings remain \n")
    }
  }
} else{
  stop("Method for dealing with probes mapped to the same gene IDs not recognized\n",
       call. = F)
}

# Output annotation report to the specified QC directory
summDir = paste(qcDir, "summary_report/", sep = "")
if (!file.exists(summDir)){ # Create a summary report directory within qcDir if it does not exist yet
  dir.create(summDir)
}

cat("\tUnmapped probes removed:", noIDCnt, "\n")
cat("\tDuplicated probes removed:", nDups, "\n\n")
cat("Annotated probes remaining:", nrow(normVals), "\n\n")

AR = c(
  AR1,
  paste("Unmapped probes removed:", noIDCnt),
  paste("Duplicated probes removed:", nDups),
  paste("Annotated probes remaining:", nrow(normVals))
)
if (glAn != FALSE) {
  write.table(
    AR,
    file = paste(summDir, glAn, "_annotReport.txt", sep = ""),
    quote = F,
    col.names = F,
    row.names = F
  )
  cat("Annotation report generated!",paste(summDir, glAn, "_annotReport.txt", sep = ""),"\n")
} else {
  write.table(
    AR,
    file = paste(summDir,"annotReport.txt", sep = ""),
    quote = F,
    col.names = F,
    row.names = F
  )
  cat("Annotation report generated!",paste(summDir,"annotReport.txt", sep = ""),"\n")
}

# Save filtered expression values
outFH = opt$output
eset = normVals # Standardizing variable naming convention between scripts
colnames(eset) = gsub("\\.","-",colnames(eset)) # Keep the sample names standardized (if data read in as a text file, hyphens are swapped for periods)
if (opt$outType == "both") {
  save(eset, file = paste(outFH, ".rda", sep = ""))
  write.table(
    data.frame("ID" = row.names(eset),eset), # provides the rownames as a labeled column in the saved output
    row.names = F,
    file = paste(outFH, ".txt", sep = ""),
    sep = "\t",
    quote = F
  )
  cat("Success! Annotated data saved to", outFH, "as both a .txt and a .RData file \n")
} else if (opt$outType == "R") {
  save(eset, file = paste(outFH, ".rda", sep = ""))
  cat("Success! Annotated data saved to", outFH, "as a .RData file \n")
} else if (opt$outType == "txt") {
  write.table(
    data.frame("ID" = row.names(eset),eset), # provides the rownames as a labeled column in the saved output
    row.names = F,
    file = paste(outFH, ".txt", sep = ""),
    sep = "\t",
    quote = F
  )
  cat("Success! Annotated data saved to", outFH, "as a .txt file \n")
} else{
  print_help(opt_parser)
  stop("Help, I don't know how to save this data!\n", call. = F)
}

