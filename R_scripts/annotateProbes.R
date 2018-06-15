#!/usr/bin/env Rscript

# source("http://bioconductor.org/biocLite.R")
# biocLite("mogene10sttranscriptcluster.db")

require(mogene10sttranscriptcluster.db)

# Double check database if running interactively
# ls("package:mogene10sttranscriptcluster.db") # List of R objects in the package
# mogene10sttranscriptcluster() # QC info

#eset = read.delim("",sep="\t",header = T, stringsAsFactors = F)

# Mapping Affy transcript cluster IDs to RefSeq names from the imported library
ID = featureNames(eset) # Pulls out the AffyIDs
mapFun = function(id){ # Function to match the primary RefSeq ID for a given AffyID and return NA in all other cases
  return(tryCatch(get(id, env=mogene10sttranscriptclusterREFSEQ)[1], error=function(e) NA))
}
RefSeq = lapply(ID,FUN = mapFun) # Applying mapFun to all AffyIDs

# Replace AffyIDs with RefSeq IDs, drop probes w/o RefSeq IDs?
