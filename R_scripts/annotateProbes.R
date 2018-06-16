#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("mogene10sttranscriptcluster.db")
# biocLite("affy")

# Read options
option_list=list(
  make_option(c("-i","--input"),type="character",help="Name of (or path to) the input file"),
  make_option(c("-a","--annotationLibrary"),type="character",help="If you know the annotation db library, specify it here. Otherwise, this script will try and detect it")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

if (!is.null(opt$annotationLibrary)){
  aLib = opt$annotationLibrary
  tryCatch({
    biocLite(as.character(aLib))
    suppressPackageStartupMessages(require(aLib))
    },error=function(e){
    stop("Specified annotation library couldn't be loaded", call. = F)
    })
}

# require(mogene10sttranscriptcluster.db)

# Double check database if running interactively
# ls("package:mogene10sttranscriptcluster.db") # List of R objects in the package
# mogene10sttranscriptcluster() # QC info

inFH = opt$input
if(grepl(".rda",inFH) == T){
  eset = load(inFH)
}else if(grepl(".txt",inFH) == T){
  tryCatch({eset = read.delim(inFH,header=T,sep = "\t")}, error=function(e){
    stop(".txt file was not recognized", call. = F)
  })
}else{
  stop("Plese provide a .rda or .txt input file")
}


# Mapping Affy transcript cluster IDs to RefSeq names from the imported library
ID = featureNames(eset) # Pulls out the AffyIDs
mapFun = function(id){ # Function to match the primary RefSeq ID for a given AffyID and return NA in all other cases
  return(tryCatch(get(id, env=mogene10sttranscriptclusterREFSEQ)[1], error=function(e) NA))
}
RefSeq = lapply(ID,FUN = mapFun) # Applying mapFun to all AffyIDs

# Replace AffyIDs with RefSeq IDs, drop probes w/o RefSeq IDs?
normVals = exprs(eset)
normVals = normVals[!is.na(RefSeq),]
rownames(normVals) = paste(rownames(normVals),RefSeq[!is.na(RefSeq)])

filtPCA = prcomp(normVals)
png("~/Desktop/filtPCA.png",width=800,height = 800)
plot(filtPCA$rotation[,1],filtPCA$rotation[,2],col=color[1:length(celFiles)],pch=16,
     xlab = paste("PC1, ",round(summary(filtPCA)$importance["Proportion of Variance",1]*100,digits = 1),"% of variance",sep=""),
     ylab = paste("PC2, ",round(summary(filtPCA)$importance["Proportion of Variance",2]*100,digits = 1),"% of variance",sep=""),
     main="PCA of normalized, filtered probes"
)
text(filtPCA$rotation[,1],filtPCA$rotation[,2],labels = sampNames, cex = 1,pos = 3)
dev.off()
