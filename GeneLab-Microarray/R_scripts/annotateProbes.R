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
option_list=list(
  make_option(c("-i","--input"),type="character",help="Name of (or path to) the input file (\\t delimited .txt file)"),
  make_option(c("-a","--arrayInfo"),type="character",default="arrayInfo.txt",help="Name of (or path to) a file containing the array information [Line 1: Manufacturer, line 2: Array version]"),
  make_option(c("-o","--output"),type="character",default="annotExpValues.txt",help="Name of (or path to) file to write results to (default: annotExpValues.txt)"),
  make_option(c("-q","--QCoutput"),type="logical",default=TRUE,help="Output QC_reporting directory of QC plots (default = TRUE)"),
  make_option("--GLDS",type="character",help="GLDS accession number for plot outputs (ie '21' for GLDS-21)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

suppressPackageStartupMessages(library("genefilter"))

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
arrayNames = c("MoGene-1_0-st-v1",
               "MOE430A",
               "Drosophila_2",
               "HG-U133_Plus_2",
               "ATH1-121501",
               "HuGene-1_0-st-v1",
               "Yeast_2",
               "Rat230_2")

arrPackages = c("mogene10sttranscriptcluster.db",
                "moe430a.db",
                "drosophila2.db",
                "hgu133plus2.db",
                "ath1121501.db",
                "hugene10sttranscriptcluster.db",
                "yeast2.db",
                "rat2302.db")

# Call the appropriate annotation package
tryCatch({
    annotPack = arrPackages[grep(pattern = arrVer,x = arrayNames,ignore.case = T)] # Pick out appropriate package by the array version
    suppressPackageStartupMessages(library(annotPack,character.only = T)) # Load selected package
    packObjs = ls(paste("package:",as.character(annotPack),sep="")) # Stores a list of all the objects in the selected package
    if(any(grepl(pattern = "REFSEQ",x = packObjs, ignore.case = T))){
      annotEnv = packObjs[grepl(pattern = "REFSEQ",x = packObjs, ignore.case = T)] # Select the enivornment from the package to map probes to RefSeq IDs
    }else if(annotPack == "ath1121501.db"){
      annotEnv = packObjs[grepl(pattern = "ACCNUM",x = packObjs, ignore.case = T)] # Select the enivornment from the package to map probes to RefSeq IDs
    }else if(annotPack == "yeast2.db"){
      annotEnv = packObjs[grepl(pattern = "ORF",x = packObjs, ignore.case = T)] # Select the enivornment from the package to map probes to RefSeq IDs
    }
    cat("Annotating with R package",annotPack,"using object:",annotEnv,"\n")
  }, error=function(e){
    stop("Array version wasn't not recognized or the annotation package was unable to load.\n
         Check that the appropriate packages are installed and the array version is contained in the list of known arrays\n", call. = F)
  }
)

# inFH = "exprsValues.txt"
inFH = opt$input
tryCatch({
  eset = read.delim(inFH,header=T,sep = "\t",stringsAsFactors = F)
  # rownames(eset) = eset[,1]
  # eset[,1] = NULL
  neset = new("ExpressionSet",exprs = as.matrix(eset))
  neset@annotation = annotPack
}, error=function(e){
  stop("Input file was not recognized", call. = F)
})

cat("Filtering out unannotated probes...\n")
filt = nsFilter(neset, var.filter = F,require.entrez = T, remove.dupEntrez = T)
nDups = filt[[2]]$numDupsRemoved # Number of probes removed that map to non-unique gene IDs
# filt[[2]]

# Mapping probe IDs to RefSeq names from the imported library
mapFun = function(id, environ){ # Function to match the primary RefSeq ID for a given probe ID and return NA in all other cases
  return(tryCatch(get(id, env=environ)[1], error=function(e) NA))
}

filtID = featureNames(filt[[1]]) # Pulls out the probe IDs

cat("Mapping probes IDs to RefSeq IDs...\n")
filtRefSeq = lapply(filtID,FUN = mapFun, environ= eval(parse(text=annotEnv))) # Applying mapFun to all probe IDs

cat("\tDuplicated probes removed:",nDups,"\n")
cat("\tUnampped probes removed:",nrow(eset)-sum(!is.na(filtRefSeq))-nDups,"\n")
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
  sampNames = gsub(".CEL","",sampNames)
  if (is.null(opt$GLDS)){ # Include GLDS accession number in plot titles if provided
    print_help(opt_parser)
    glAn = ''
    cat("Warning: Generating plots without accession number")
  }else{
    glAn = paste('GLDS-',opt$GLDS,sep='')
  }
  if(!file.exists(paste('./',glAn,'_QC_reporting/',sep=''))) dir.create(paste('./',glAn,'_QC_reporting/',sep=''))
  # Post-normalization QC
  cat("Post annotation/filtering QC...\n")
  
  # Density distributions
  png(paste("./",glAn,"_QC_reporting/",glAn,"_microarray_filtDensityDistributions.png",sep=""),width=800,height=800)
  ylims = c(0,.8)
  xlims = c(0,16)
  for(i in 1:ncol(normVals)){
    if(i == 1){
      plot(density(normVals[,i]),ylim = ylims,xlim=xlims
           ,xlab='Normalized, annotated expression values[log2]'
           ,main=paste(glAn,' Normalized, filtered expression distributions',sep=''),col=color[i])
      par(new=T)
    }else{
      plot(density(normVals[,i]),ylim = ylims,xlim=xlims,axes=F,xlab='',ylab='',main='',col=color[i])
      par(new=T)
    }
  }
  legend(13,0.8,col=color[1:length(sampNames)],legend=sampNames
         ,pch=15,bty = "n",cex = 0.9,pt.cex = 0.8,y.intersp = 0.8)
  garbage <- dev.off()
  
  # PCA plot
  filtPCA = prcomp(normVals)
  png(paste("./",glAn,"_QC_reporting/",glAn,"_microarray_filtPCA.png",sep=""),width=800,height = 800)
  plot(filtPCA$rotation[,1],filtPCA$rotation[,2],col=color[1:length(sampNames)],pch=16,
       xlab = paste("PC1, ",round(summary(filtPCA)$importance["Proportion of Variance",1]*100,digits = 1),"% of variance",sep=""),
       ylab = paste("PC2, ",round(summary(filtPCA)$importance["Proportion of Variance",2]*100,digits = 1),"% of variance",sep=""),
       main=paste(glAn," PCA of normalized, filtered probes", sep = "")
  )
  text(filtPCA$rotation[,1],filtPCA$rotation[,2],labels = sampNames, cex = 1,pos = 3)
  garbage <- dev.off()
}

