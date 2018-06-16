#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("affy")
# biocLite("affyPLM")
# biocLite("oligo")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list=list(
  make_option(c("-n","--normalization"),type="character",default="rma",help="Normalization method [rma (default, full rma), quantile (no background correction), background (no quant. normalization), log2 (no quant. norm. or background correction)"),
  make_option(c("-o","--outFile"),type="character",default="expValues",help="Name of the output file [without extension!] (default: expValues"),
  make_option(c("-t","--outType"),type="character",default="both",help="Format of output data: R (Rdata eset object), txt (tab delimited file with identifiers and sample names), both (default)"),
  make_option("--outputData",type="logical",default=TRUE,help="Output data at all (default TRUE)"),
  make_option("--QCoutput",type="logical",default=TRUE,help="Output QC_reporting directory of QC plots (default = TRUE)"),
  make_option(c("-np","--NUSEplot"),type="logical",default=FALSE,help="Include a NUSE plot in the QC output, adds significantly to runtime (default = FALSE)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

norm = opt$normalization
QCout = opt$QCoutput
NUSEplot = opt$NUSEplot

detach_package <- function(pkg, character.only = FALSE){
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
sampNames = gsub("./","",celFiles)
sampNames = gsub(".CEL","",sampNames)

tryCatch({raw = ReadAffy()}, error=function(e){
  stop("No .CEL files detected in the current directory", call. = F)
  })

write.table(c("Affymetrix",as.character(raw@cdfName)),file = "arrayInfo.txt",quote = F,
            col.names = F, row.names = F)

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
#require(oligo)


## Raw QC
if(QCout == T){
  # Prepare plotting options
  toMatch = c(8,183,31,45,51,100,101,118,128,139,147,183,254,421,467,477,
              483,493,498,503,508,535,552,575,635,655)
  color = grDevices::colors()[rep(toMatch,3)] # Create a library of colors for plotting
  if(!file.exists('./QC_reporting/')) dir.create('./QC_reporting/')
  
  #Images
  cat("Generating raw images")
  if(st == T){
    for(i in 1:length(celFiles)){
      png(paste('./QC_reporting/image_',sampNames[i],'.png',sep=''),width=800, height = 800)
      image(raw, which = i)
      dev.off()
      cat(".")
    }
  }else{
    nblines=length(celFiles)%/%4 + as.numeric((length(celFiles)%%4)!=0)
    png('./QC_reporting/image.png',width=800,height = 200*nblines)
    par(mfrow=c(nblines,4))
    image(raw)
    dev.off()
  }
  cat("\n")

  #MA plot
  cat("Generating raw data MA plots...\n")
  nblines=length(celFiles)%/%3 + as.numeric((length(celFiles)%%3)!=0)
  png("./QC_reporting/rawPlotMA.png",width=800, height = 300*nblines )
  par(mfrow=c(nblines,3))
  if(st == T){
    MAplot(raw)
  }else{
    MAplot(raw,type="pm")
  }
  dev.off()
  
  # Intensity distributions of the pm probes from each microarray on the same graph
  cat("Generating initial distribution plots")
  mypms = pm(raw)
  png("./QC_reporting/rawDensityDistributions.png",width=800,height=800 )
  ylims = c(0,.8)
  xlims = c(0,16)
  for(i in 1:ncol(mypms)){
    cat(".")
    if(i == 1){
      plot(density(log2(mypms[,i])),ylim = ylims,xlim=xlims,xlab='log2(Raw Intensities)',main='Raw intensity distributions',col=color[i])
      par(new=T)
    }else{
      plot(density(log2(mypms[,i])),ylim = ylims,xlim=xlims,axes=F,xlab='',main='',col=color[i])
      par(new=T)
    }
  }
  legend(13,0.8,col=color[1:length(celFiles)],legend=sampNames
         ,pch=15,bty = "n",cex = 0.9,pt.cex = 0.8,y.intersp = 0.8)
  cat("\n")
  dev.off()
  
  # Boxplots
  png("./QC_reporting/rawBoxplot.png",width=800,height = 400)
  par(mar=c(7,5,1,1))
  if(st == T){
    boxplot(oligo::rma(raw, background=FALSE, normalize=FALSE, subset=NULL, target="core"), las=2,
            names = sampNames, main="Raw intensities",col=color[1:length(celFiles)])
  }else{
    boxplot(raw,las=2,outline=FALSE,col=color[1:length(celFiles)],main="Raw intensities",names=sampNames)
  }
  mtext(text="log2 Intensity", side=2, line=2.5, las=0)
  dev.off()
  
  # PCA
  cat("Performing PCA of raw data...\n")
  rawPCA = prcomp(mypms)
  png("./QC_reporting/rawPCA.png",width=800,height = 800)
  plot(rawPCA$rotation[,1],rawPCA$rotation[,2],col=color[1:length(celFiles)],pch=16,
       xlab = paste("PC1, ",round(summary(rawPCA)$importance["Proportion of Variance",1]*100,digits = 1),"% of variance",sep=""),
       ylab = paste("PC2, ",round(summary(rawPCA)$importance["Proportion of Variance",2]*100,digits = 1),"% of variance",sep=""),
       main="PCA of raw data"
  )
  text(rawPCA$rotation[,1],rawPCA$rotation[,2],labels = sampNames, cex = 1,pos = 3)
  dev.off()
  
  if(NUSEplot == T){
    #NUSE plot
    cat("Fitting probe-level model and generating NUSE plot...\n")
    png("./QC_reporting/NUSE.png",width=800,height = 600)
    if(st == T){
      Pset = fitProbeLevelModel(raw)
      NUSE(Pset, col = color[1:length(sampNames)], las=2)
    }else{
      Pset=fitPLM(raw)
      NUSE(Pset,col = color[1:length(sampNames)], las=2)
    }
    title(main="NUSE plot of microarray experiments")
    abline(h=1.1,lty=1,col="red")
    dev.off()
  }
}


## Normalize
cat("Normalizing with selected normalization technique...\n")
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

outFH = opt$outFile
if(opt$outputData == TRUE){
  if(opt$outType == "both"){
    save(eset,file=paste(outFH,".rda",sep=""))
    write.exprs(eset,file=paste(outFH,".txt",sep=""),sep="\t")
  }else if(opt$outType == "R"){
    save(eset,file=paste(outFH,".rda",sep=""))
  }else if(opt$outType == "txt"){
    write.exprs(eset,file=paste(outFH,".txt",sep=""),sep="\t")
  }else{
    stop("Help, I don't know how to save this data!")
  }
}

if(QCout == T){
  cat("Post normalization QC steps...\n")
  # Post-normalization QC
  png("./QC_reporting/normDensityDistributions.png",width=800,height=800 )
  ylims = c(0,.8)
  xlims = c(0,16)
  normVals = exprs(eset)
  for(i in 1:ncol(normVals)){
    if(i == 1){
      plot(density(normVals[,i]),ylim = ylims,xlim=xlims,xlab='Normalized expression values[log2]',main='Normalized expression distributions',col=color[i])
      par(new=T)
    }else{
      plot(density(normVals[,i]),ylim = ylims,xlim=xlims,axes=T,xlab='',main='',col=color[i])
      par(new=T)
    }
  }
  legend(13,0.8,col=color[1:length(celFiles)],legend=sampNames
         ,pch=15,bty = "n",cex = 0.9,pt.cex = 0.8,y.intersp = 0.8)
  dev.off()
  
  # Boxplots
  png("./QC_reporting/normBoxplot.png",width=800,height = 400)
  par(mar=c(7,5,1,1))
  if(st == T){
    boxplot(normVals,las=2,outline=FALSE,col=color[1:length(celFiles)],main="Normalized intensities",transfo='identity',names=sampNames)
    mtext(text="log2 Intensity", side=2, line=2.5, las=0) 
  }else{
    boxplot(normVals,las=2,outline=FALSE,col=color[1:length(celFiles)],main="Normalized intensities",names=sampNames)
    mtext(text="log2 Intensity", side=2, line=2.5, las=0)
  }
  dev.off()
  
  #MA plot
  cat("Generating MA plots from the normalized data...\n")
  nblines=length(celFiles)%/%3 + as.numeric((length(celFiles)%%3)!=0)
  png("./QC_reporting/normPlotMA.png",width=800, height = 300*nblines )
  par(mfrow=c(nblines,3))
  MAplot(eset)
  dev.off()
  
  # PCA
  cat("Performing PCA of normalized data...\n")
  normPCA = prcomp(normVals)
  png("./QC_reporting/normPCA.png",width=800,height = 800)
  plot(normPCA$rotation[,1],normPCA$rotation[,2],col=color[1:length(celFiles)],pch=16,
       xlab = paste("PC1, ",round(summary(normPCA)$importance["Proportion of Variance",1]*100,digits = 1),"% of variance",sep=""),
       ylab = paste("PC2, ",round(summary(normPCA)$importance["Proportion of Variance",2]*100,digits = 1),"% of variance",sep=""),
       main="PCA of normalized data"
  )
  text(normPCA$rotation[,1],normPCA$rotation[,2],labels = sampNames, cex = 1,pos = 3)
  dev.off()
}