#!/usr/bin/env Rscript

# install.packages("optparse")
# install.packages("hash")
# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")
# biocLite("mogene10sttranscriptcluster.db")

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("genefilter"))

# Read options
option_list=list(
  make_option(c("-i","--input"),type="character",help="Name of (or path to) the input file"),
  make_option(c("-a","--arrayInfo"),type="character",default="arrayInfo.txt",help="Name of (or path to) a file containing the array information [Line 1: Manufacturer, line 2: Array version]"),
  make_option(c("-o","--output"),type="character",default="annotExpValues.txt",help="Name of (or path to) file to write results to (default: annotExpValues.txt)"),
  make_option(c("-q","--QCoutput"),type="logical",default=TRUE,help="Output QC_reporting directory of QC plots (default = TRUE)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

# aiFH = "arrayInfo.txt"
aiFH = opt$arrayInfo

tryCatch({
    aInf = read.delim(aiFH,header=F,stringsAsFactors = F)
    arrMan = aInf[1,1]
    arrVer = aInf[2,1]
  }, error=function(e){
    stop("Array info file not found or organization not recognized, check options help", call. = F)
  }
)

# Set-up array version:annotation database pseudo-dictionary
arrays = c("MoGene-1_0-st-v1")
arrPackages = c("mogene10sttranscriptcluster.db")

# Call the appropriate annotation package
tryCatch({
    annotPack = arrPackages[grep(pattern = arrVer,x = arrays,ignore.case = T)]
    suppressPackageStartupMessages(library(annotPack,character.only = T))
    packObjs = ls(paste("package:",as.character(annotPack),sep=""))
    annotEnv = packObjs[grepl(pattern = "REFSEQ",x = packObjs, ignore.case = T)]
  }, error=function(e){
    stop("Error: Array version wasn't not recognized or the annotation package was unable to load.\n
         Check that the appropriate packages are installed and the array version is contained in the list of known arrays", call. = F)
  }
)


# Double check database if running interactively
# ls("package:mogene10sttranscriptcluster.db") # List of R objects in the package
# mogene10sttranscriptcluster() # QC info

# inFH = "expValues.txt"
inFH = opt$input
if(grepl(".rda",inFH) == T){
  eset = load(file = inFH)
}else if(grepl(".txt",inFH) == T){
  tryCatch({
    eset = read.delim(inFH,header=T,sep = "\t")
    row.names(eset) = eset[,1]
    eset[,1] = NULL
    neset = new("ExpressionSet",exprs = as.matrix(eset))
    neset@annotation = annotPack
  }, error=function(e){
    stop(".txt file was not recognized", call. = F)
    })
  }else{
  stop("Plese provide a .rda or .txt input file")
}

cat("Filtering out unmapped probes...\n")
filt = nsFilter(neset, var.filter = F,require.entrez = T,remove.dupEntrez = T)
cat("\tUnampped probes removed:",filt[[2]]$numRemoved.ENTREZID,"\n")
cat("\tDuplicated probes removed:",filt[[2]]$numDupsRemoved,"\n")
print(filt[[2]])

# Mapping probe IDs to RefSeq names from the imported library
cat("Mapping probes IDs to RefSeq IDs...\n")
mapFun = function(id, environ){ # Function to match the primary RefSeq ID for a given AffyID and return NA in all other cases
  return(tryCatch(get(id, env=environ)[1], error=function(e) NA))
}


filtID = featureNames(filt[[1]]) # Pulls out the probe IDs

filtRefSeq = lapply(filtID,FUN = mapFun, environ= eval(parse(text=annotEnv))) # Applying mapFun to all AffyIDs

length(filtRefSeq)
cat("Annotated probes remaining:",sum(!is.na(filtRefSeq)),"\n")
if(sum(!is.na(filtRefSeq)) > length(unique(filtRefSeq[!is.na(filtRefSeq)]))){
  cat("\n\tWarning: non-unique probe to RefSeq mappings encountered \n")
}

# Replace AffyIDs with RefSeq IDs, drop probes w/o RefSeq IDs
normVals = exprs(filt[[1]])
normVals = normVals[!is.na(filtRefSeq),]
rownames(normVals) = filtRefSeq[!is.na(filtRefSeq)]

# Save filtered expression values to working directory
outFH = opt$output
write.table(normVals,file=outFH,sep="\t",quote = F)

if(opt$QCoutput == T){
  # Prepare plotting options
  toMatch = c(8,183,31,45,51,100,101,118,128,139,147,183,254,421,467,477,
              483,493,498,503,508,535,552,575,635,655)
  color = grDevices::colors()[rep(toMatch,3)] # Create a library of colors for plotting
  sampNames = colnames(normVals)
  if(!file.exists('./QC_reporting/')) dir.create('./QC_reporting/')
  
  # Post-normalization QC
  cat("Post annotation/filtering QC...\n")
  
  # Density distributions
  png("./QC_reporting/filtDensityDistributions.png",width=800,height=800 )
  ylims = c(0,.8)
  xlims = c(0,16)
  for(i in 1:ncol(normVals)){
    if(i == 1){
      plot(density(normVals[,i]),ylim = ylims,xlim=xlims,xlab='Normalized, annotated expression values[log2]',main='Normalized expression distributions',col=color[i])
      par(new=T)
    }else{
      plot(density(normVals[,i]),ylim = ylims,xlim=xlims,axes=T,xlab='',main='',col=color[i])
      par(new=T)
    }
  }
  legend(13,0.8,col=color[1:length(sampNames)],legend=sampNames
         ,pch=15,bty = "n",cex = 0.9,pt.cex = 0.8,y.intersp = 0.8)
  dev.off()
  
  # PCA plot
  filtPCA = prcomp(normVals)
  png("./QC_reporting/filtPCA.png",width=800,height = 800)
  plot(filtPCA$rotation[,1],filtPCA$rotation[,2],col=color[1:length(sampNames)],pch=16,
       xlab = paste("PC1, ",round(summary(filtPCA)$importance["Proportion of Variance",1]*100,digits = 1),"% of variance",sep=""),
       ylab = paste("PC2, ",round(summary(filtPCA)$importance["Proportion of Variance",2]*100,digits = 1),"% of variance",sep=""),
       main="PCA of normalized, filtered probes"
  )
  text(filtPCA$rotation[,1],filtPCA$rotation[,2],labels = sampNames, cex = 1,pos = 3)
  dev.off()
}

