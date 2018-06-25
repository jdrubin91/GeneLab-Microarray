#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("affy")
# biocLite("affyPLM")
# biocLite("oligo")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list=list(
  make_option(c("-i","--input"),type="character",help="Path to directory containing input .CEL files"),
  make_option(c("-n","--normalization"),type="character",default="rma",help="Normalization method [rma (default, full rma), quantile (no background correction), background (no quant. normalization), log2 (no quant. norm. or background correction)"),
  make_option(c("-o","--outFile"),type="character",default="expValues",help="Name of the output file [without extension!] (default: expValues)"),
  make_option(c("-t","--outType"),type="character",default="txt",help="Format of output data: R (Rdata eset object), txt (default, tab delimited file with identifiers and sample names), both"),
  make_option("--outputData",type="logical",default=TRUE,help="Output data at all (default TRUE)"),
  make_option(c("-a","--arrayInfoOnly"),type="logical",default=FALSE,help="Detect-affy-array-only mode. If true, script will exit after outputting the arrayInfo file. (Default: FALSE)"),
  make_option("--QCoutput",type="logical",default=TRUE,help="Output QC_reporting directory of QC plots (default = TRUE)"),
  make_option("--NUSEplot",type="logical",default=FALSE,help="Include a NUSE plot in the QC output, adds significantly to runtime (default = FALSE)"),
  make_option("--GLDS",type="character",help="GLDS accession number for plot outputs (ie '21' for GLDS-21)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

norm = opt$normalization
QCout = opt$QCoutput
NUSEplot = opt$NUSEplot

if (is.null(opt$GLDS)){ # Include GLDS accession number in outputs if provided
  print_help(opt_parser)
  glAn = ''
  cat("Warning: No GLDS accession number provided")
}else{
  glAn = paste('GLDS-',opt$GLDS,sep='')
}

if (is.null(opt$input)){ # Include GLDS accession number in outputs if provided
  print_help(opt_parser)
  stop("No path to input directory provided. Please look over the available options", call. = F)
}else{
  inPath = opt$input
  setwd(inPath)
}

detach_package = function(pkg, character.only = FALSE){
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

# Load initial libraries
suppressPackageStartupMessages(require(affy))

# setwd("~/Documents/genelab/rot1/GLDS4/microarray/")
celFiles <- list.celfiles(full.names=TRUE)
sampNames = gsub(".*_microarray_","",celFiles)
sampNames = gsub(".CEL","",sampNames)
sampNames = gsub("./","",sampNames) # Extract sample naems form the list of .CEL files

tryCatch({raw = ReadAffy()}, error=function(e){
  stop("No .CEL files detected in the current directory", call. = F)
  })

# Output array information to a separate file
write.table(c("Affymetrix",as.character(raw@cdfName)),file = paste(glAn,"_arrayInfo.txt",sep=""),quote = F,
            col.names = F, row.names = F)

# Exit script if arrayInfoOnly mode is True
if(opt$arrayInfoOnly == TRUE) stop("Detect-affy-array-only mode is on, exiting...", call. = F)

if (grepl("-st-",raw@cdfName,ignore.case = T)){
  detach_package(affy)
  rm(raw)
  suppressPackageStartupMessages(require(oligo))
  raw = read.celfiles(celFiles)
  st = T
}else{
  suppressPackageStartupMessages(require(affyPLM))
  st = F
}

## Raw QC
if(QCout == T){
  cat("Performing intial QC\n")
  # Prepare plotting options
  toMatch = c(8,183,31,45,51,100,101,118,128,139,147,183,254,421,467,477,
              483,493,498,503,508,535,552,575,635,655)
  color = grDevices::colors()[rep(toMatch,3)] # Create a library of colors for plotting
  if(!file.exists(paste('./',glAn,'_QC_reporting/',sep=''))) dir.create(paste('./',glAn,'_QC_reporting/',sep=''))
  
  #Images
  cat("\tGenerating raw images")
  if(st == T){
    for(i in 1:length(celFiles)){
      png(paste('./',glAn,'_QC_reporting/',glAn,'_',sampNames[i],'_image.png',sep=''),width=800, height = 800)
      image(raw, which = i)
      garbage <- dev.off()
      cat(".")
    }
  }else{
    nblines=length(celFiles)%/%4 + as.numeric((length(celFiles)%%4)!=0)
    png(paste('./',glAn,'_QC_reporting/',glAn,'_images.png',sep=''),width=800,height = 200*nblines)
    par(mfrow=c(nblines,4))
    image(raw)
    garbage <- dev.off()
  }
  cat("\n")

  #MA plot
  cat("\tGenerating raw data MA plots...\n")
  nblines=length(celFiles)%/%3 + as.numeric((length(celFiles)%%3)!=0)
  png(paste('./',glAn,'_QC_reporting/',glAn,'_rawPlotMA.png',sep=''),width=800, height = 300*nblines )
  par(mfrow=c(nblines,3))
  if(st == T){
    MAplot(raw)
  }else{
    MAplot(raw,type="pm")
  }
  garbage <- dev.off()
  
  # Intensity distributions of the pm probes from each microarray on the same graph
  cat("\tGenerating initial distribution plots")
  mypms = pm(raw)
  png(paste('./',glAn,'_QC_reporting/',glAn,'_rawDensityDistributions.png',sep=''),width=800,height=800 )
  ylims = c(0,.8)
  xlims = c(0,16)
  for(i in 1:ncol(mypms)){
    cat(".")
    if(i == 1){
      plot(density(log2(mypms[,i])),ylim = ylims,xlim=xlims,xlab='log2(Raw Intensities)',main=paste(glAn,' Raw intensity distributions',sep=''),col=color[i])
      par(new=T)
    }else{
      plot(density(log2(mypms[,i])),ylim = ylims,xlim=xlims,axes=F,xlab='',ylab='',main='',col=color[i])
      par(new=T)
    }
  }
  legend(13,0.8,col=color[1:length(celFiles)],legend=sampNames
         ,pch=15,bty = "n",cex = 0.9,pt.cex = 0.8,y.intersp = 0.8)
  garbage <- dev.off()
  cat("\n")
  
  # Boxplots
  png(paste('./',glAn,'_QC_reporting/',glAn,'_rawBoxplot.png',sep=''),width=800,height = 400)
  par(mar=c(7,5,1,1))
  if(st == T){
    boxplot(oligo::rma(raw, background=FALSE, normalize=FALSE, subset=NULL, target="core"), las=2,
            names = sampNames, main=paste(glAn," Raw intensities",sep=""),col=color[1:length(celFiles)]);
  }else{
    boxplot(raw,las=2,outline=FALSE,col=color[1:length(celFiles)],main = paste(glAn," Raw intensities",sep=""),names=sampNames)
  }
  mtext(text="log2 Intensity", side=2, line=2.5, las=0)
  garbage <- dev.off()
  
  # PCA
  cat("\tPerforming PCA of raw data...\n")
  rawPCA = prcomp(mypms)
  png(paste('./',glAn,'_QC_reporting/',glAn,'_rawPCA.png',sep=''),width=800,height = 800)
  plot(rawPCA$rotation[,1],rawPCA$rotation[,2],col=color[1:length(celFiles)],pch=16,
       xlab = paste("PC1, ",round(summary(rawPCA)$importance["Proportion of Variance",1]*100,digits = 1),"% of variance",sep=""),
       ylab = paste("PC2, ",round(summary(rawPCA)$importance["Proportion of Variance",2]*100,digits = 1),"% of variance",sep=""),
       main=paste(glAn," PCA of raw data",sep="")
  )
  text(rawPCA$rotation[,1],rawPCA$rotation[,2],labels = sampNames, cex = 1,pos = 3)
  garbage <- dev.off()
  
  #NUSE plot
  if(NUSEplot == T){
    cat("\tFitting probe-level model and generating RLE/NUSE plots...\n")
    if(st == T){
      Pset = fitProbeLevelModel(raw)
      # RLE plot
      png(paste('./',glAn,'_QC_reporting/',glAn,'_RLE.png',sep=''),width=800,height = 600)
      par(mar=c(7,5,1,1))
      RLE(Pset, col = color[1:length(sampNames)],
          names = sampNames, las=2, main=paste(glAn," Relative Log Expression (RLE) plot",sep=""))
      abline(h=0,lty=1,col="red")
      garbage <- dev.off()
      # NUSE plot
      png(paste('./',glAn,'_QC_reporting/',glAn,'_NUSE.png',sep=''),width=800,height = 600)
      par(mar=c(7,5,1,1))
      NUSE(Pset, col = color[1:length(sampNames)], las=2)
      title(main=paste(glAn," NUSE plot of microarray experiments",sep=""))
      abline(h=1.1,lty=1,col="red")
      garbage <- dev.off()
    }else{
      Pset=fitPLM(raw)
      # RLE plot
      png(paste('./',glAn,'_QC_reporting/',glAn,'_RLE.png',sep=''),width=800,height = 600)
      par(mar=c(7,5,1,1))
      RLE(Pset, col = color[1:length(sampNames)],
          names = sampNames, las=2, main="Relative Log Expression (RLE) plot")
      abline(h=0,lty=1,col="red")
      garbage <- dev.off()
      # NUSE plot
      png(paste('./',glAn,'_QC_reporting/',glAn,'_NUSE.png',sep=''),width=800,height = 600)
      par(mar=c(7,5,1,1))
      NUSE(Pset,col = color[1:length(sampNames)], las=2)
      title(main=paste(glAn," NUSE plot of microarray experiments",sep=""))
      abline(h=1.1,lty=1,col="red")
      garbage <- dev.off()
    }
  }
}

outFH = opt$outFile
if(opt$outputData == TRUE){
  
  ## Normalize
  cat("\nNormalizing with selected normalization technique...\n")
  if(norm=='rma'){
    eset = rma(raw)
  }else if(norm=='quantile'){
    eset = rma(raw, background = F, normalize = T)
  }else if(norm=='background'){
    eset = rma(raw, background = T, normalize = F)
  }else if(norm=='log2'){
    eset = rma(raw, background = F, normalize = F)
  }else{
    stop("Normalization did not occur, please examine script inputs and default values",call. = F)
  }
  
  if(opt$outType == "both"){
    save(eset,file=paste(outFH,".rda",sep=""))
    write.table(exprs(eset),file=paste(outFH,".txt",sep=""),sep="\t",quote = F)
  }else if(opt$outType == "R"){
    save(eset,file=paste(outFH,".rda",sep=""))
  }else if(opt$outType == "txt"){
    write.table(exprs(eset),file=paste(outFH,".txt",sep=""),sep="\t",quote = F)
  }else{
    print_help(opt_parser)
    stop("Help, I don't know how to save this data!",call. = F)
  }
}

if(QCout == T){
  cat("Post normalization QC steps...\n")
  # Post-normalization QC
  png(paste('./',glAn,'_QC_reporting/',glAn,'_normDensityDistributions.png',sep=''),width=800,height=800 )
  ylims = c(0,.8)
  xlims = c(0,16)
  normVals = exprs(eset)
  for(i in 1:ncol(normVals)){
    if(i == 1){
      plot(density(normVals[,i]),ylim = ylims,xlim=xlims,xlab='Normalized expression values[log2]',main=paste(glAn,' Normalized expression distributions',sep=''),col=color[i])
      par(new=T)
    }else{
      plot(density(normVals[,i]),ylim = ylims,xlim=xlims,axes=F,xlab='',ylab='',main='',col=color[i])
      par(new=T)
    }
  }
  legend(13,0.8,col=color[1:length(celFiles)],legend=sampNames
         ,pch=15,bty = "n",cex = 0.9,pt.cex = 0.8,y.intersp = 0.8)
  garbage <- dev.off()
  
  # Boxplots
  png(paste('./',glAn,'_QC_reporting/',glAn,'_normBoxplot.png',sep=''),width=800,height = 400)
  par(mar=c(7,5,1,1))
  if(st == T){
    boxplot(normVals,las=2,outline=FALSE,col=color[1:length(celFiles)],main=paste(glAn," Normalized intensities",sep=""),transfo='identity',names=sampNames)
    mtext(text="log2 Intensity", side=2, line=2.5, las=0) 
  }else{
    boxplot(normVals,las=2,outline=FALSE,col=color[1:length(celFiles)],main=paste(glAn," Normalized intensities",sep=""),names=sampNames)
    mtext(text="log2 Intensity", side=2, line=2.5, las=0)
  }
  garbage <- dev.off()
  
  #MA plot
  cat("\tGenerating MA plots from the normalized data...\n")
  nblines=length(celFiles)%/%3 + as.numeric((length(celFiles)%%3)!=0)
  png(paste('./',glAn,'_QC_reporting/',glAn,'_normPlotMA.png',sep=''),width=800, height = 300*nblines )
  par(mfrow=c(nblines,3))
  MAplot(eset)
  garbage <- dev.off()
  
  # PCA
  cat("\tPerforming PCA of normalized data...\n")
  normPCA = prcomp(normVals)
  png(paste('./',glAn,'_QC_reporting/',glAn,'_normPCA.png',sep=''),width=800,height = 800)
  plot(normPCA$rotation[,1],normPCA$rotation[,2],col=color[1:length(celFiles)],pch=16,
       xlab = paste("PC1, ",round(summary(normPCA)$importance["Proportion of Variance",1]*100,digits = 1),"% of variance",sep=""),
       ylab = paste("PC2, ",round(summary(normPCA)$importance["Proportion of Variance",2]*100,digits = 1),"% of variance",sep=""),
       main=paste(glAn," PCA of normalized data",sep="")
  )
  text(normPCA$rotation[,1],normPCA$rotation[,2],labels = sampNames, cex = 1,pos = 3)
  garbage <- dev.off()
}