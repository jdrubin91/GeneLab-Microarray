#!/usr/bin/env Rscript

library("optparse")
# Read options
option_list=list(
  make_option("--normalization",type="character",default="rma",help="Normalization method [rma (default, full rma), quantile (no background correction), background (no quant. normalization), log2 (no quant. norm. or background correction)"),
  make_option("--outFile",type="character",default="expValues",help="Name of the output file"),
  make_option("--outType",type="character",default="R",help="Format of output data: R (Rdata eset object), txt (tab delimited file with identifiers and sample names)"),
  make_option("--outputData",type="logical",default=TRUE,help="Output data at all"),
  make_option("--QCoutput",type="logical",default=TRUE,help="Output concatenated PDF file of QC plots (default = TRUE)"),
  make_option("--NUSEplot",type="logical",default=TRUE,help="Include a NUSE plot in the QC output, adds significantly to runtime (default = TRUE)")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

norm = opt$normalization
QCout = opt$QCoutput
NUSEplot = opt$NUSEplot

if(is.null(opt$input)){
  print_help(opt_parser)
  stop("Input required!", call.=FALSE)
}

# source("http://bioconductor.org/biocLite.R")
# biocLite("affy")
# biocLite("affyPLM")
# biocLite("oligo")

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
require(affy)

# setwd("~/Documents/genelab/rot1/GLDS4/microarray/")
celFiles <- list.celfiles(full.names=TRUE)
sampNames = gsub("./","",celFiles)
sampNames = gsub(".CEL","",sampNames)

tryCatch({raw = ReadAffy()}, error=function(e){
  stop("No .CEL files detected in the current directory", call. = F)
  })

if (grepl("-st-",raw@cdfName,ignore.case = T)){
  detach_package(affy)
  rm(raw)
  require(oligo)
  raw = read.celfiles(celFiles)
  st = T
}else{
  require(affyPLM)
  st = F
}
#require(oligo)


## Raw QC
if(QCout == T){
  # Prepare plotting options
  toMatch = c(8,183,31,45,51,100,101,118,128,139,147,183,254,421,467,477,
              483,493,498,503,508,535,552,575,635,655)
  color = grDevices::colors()[rep(toMatch,3)] # Create a library of colors for plotting
  if (NUSEplot == T){
    numPlts = 10
  }else{
    numPlts = 9
  }
  outPlts=vector(numPlts, mode = 'list')
  cntPlts=1
  
  #Images
  cat("Generating raw images\n")
  nblines=length(celFiles)%/%4 + as.numeric((length(celFiles)%%4)!=0) 
  par(mfrow=c(nblines,4))
  image(raw)
  outPlts[[cntPlts]] <- recordPlot()
  cntPlts = cntPlts + 1
  graphics.off()

    #MA plot
  cat("Generating raw data MA plots...\n")
  nblines=length(celFiles)%/%3 + as.numeric((length(celFiles)%%3)!=0) 
  par(mfrow=c(nblines,3))
  if(st == T){
    MAplot(raw)
  }else{
    MAplot(raw,type="pm")
  }
  outPlts[[cntPlts]] <- recordPlot()
  cntPlts = cntPlts + 1
  graphics.off()
  
  # Intensity distributions of the pm probes from each microarray on the same graph
  cat("Generating initial distribution plots")
  mypms = pm(raw)
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
         ,pch=15,bty = "n",cex = 0.8,pt.cex = 0.8,y.intersp = 0.6)
  cat("\n")
  outPlts[[cntPlts]] <- recordPlot()
  cntPlts = cntPlts + 1
  graphics.off()
  
  # Boxplots
  par(mar=c(7,5,1,1))
  if(st == T){
    boxplot(oligo::rma(raw, background=FALSE, normalize=FALSE, subset=NULL, target="core"), las=2,
            names = sampNames, main="Raw intensities",col=color[1:length(celFiles)])
    mtext(text="log2 Intensity", side=2, line=2.5, las=0) 
  }else{
    boxplot(raw,las=2,outline=FALSE,col=color[1:length(celFiles)],main="Raw intensities",names=sampNames)
    mtext(text="log2 Intensity", side=2, line=2.5, las=0)
  }
  outPlts[[cntPlts]] <- recordPlot()
  cntPlts = cntPlts + 1
  graphics.off()
  
  # PCA
  cat("Performing PCA of raw data...\n")
  rawPCA = prcomp(mypms)
  plot(rawPCA$rotation[,1],rawPCA$rotation[,2],col=color[1:length(celFiles)],pch=15,
       xlab = paste("PC1, ",round(summary(rawPCA)$importance["Proportion of Variance",1]*100,digits = 1),"% of variance",sep=""),
       ylab = paste("PC2, ",round(summary(rawPCA)$importance["Proportion of Variance",2]*100,digits = 1),"% of variance",sep=""),
       main="PCA of raw data"
  )
  text(rawPCA$rotation[,1],rawPCA$rotation[,2],labels = sampNames, cex = 0.7,pos = 3)
  outPlts[[cntPlts]] <- recordPlot()
  cntPlts = cntPlts + 1
  graphics.off()
  
  if(NUSEplot == T){
    #NUSE plot
    cat("Fitting probe-level model and generating NUSE plot...\n")
    if(st == T){
      Pset = fitProbeLevelModel(raw)
      NUSE(Pset, col = color[1:length(sampNames)], las=2)
    }else{
      Pset=fitPLM(raw)
      NUSE(Pset,col = color[1:length(sampNames)], las=2)
    }
    title(main="NUSE plot of microarray experiments")
    abline(h=1.1,lty=1,col="red")
    outPlts[[cntPlts]] <- recordPlot()
    cntPlts = cntPlts + 1
    graphics.off()
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
  if(opt$outType == "R"){
    save(eset,file=outFH)
  }else if(opt$outType == "txt"){
    write.exprs(eset,file=paste(outFH,".txt",sep=""),sep="\t")
  }else{
    stop("Help, I don't know how to save this data!")
  }
}

if(QCout == T){
  cat("Post normalization QC steps...\n")
  # Post-normalization QC step
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
         ,pch=15,bty = "n",cex = 0.8,pt.cex = 0.8,y.intersp = 0.6)
  outPlts[[cntPlts]] <- recordPlot()
  cntPlts = cntPlts + 1
  graphics.off()
  
  # Boxplots
  par(mar=c(7,5,1,1))
  if(st == T){
    boxplot(normVals,las=2,outline=FALSE,col=color[1:length(celFiles)],main="Normalized intensities",transfo='identity',names=sampNames)
    mtext(text="log2 Intensity", side=2, line=2.5, las=0) 
  }else{
    boxplot(normVals,las=2,outline=FALSE,col=color[1:length(celFiles)],main="Normalized intensities",names=sampNames)
    mtext(text="log2 Intensity", side=2, line=2.5, las=0)
  }
  outPlts[[cntPlts]] <- recordPlot()
  cntPlts = cntPlts + 1
  graphics.off()
  
  #MA plot
  cat("Generating MA plots from the normalized data...\n")
  nblines=length(celFiles)%/%3 + as.numeric((length(celFiles)%%3)!=0) 
  par(mfrow=c(nblines,3))
  MAplot(eset)
  outPlts[[cntPlts]] <- recordPlot()
  cntPlts = cntPlts + 1
  graphics.off()
  
  # PCA
  cat("Performing PCA of normalized data...\n")
  normPCA = prcomp(normVals)
  plot(normPCA$rotation[,1],normPCA$rotation[,2],col=color[1:length(celFiles)],pch=15,
       xlab = paste("PC1, ",round(summary(normPCA)$importance["Proportion of Variance",1]*100,digits = 1),"% of variance",sep=""),
       ylab = paste("PC2, ",round(summary(normPCA)$importance["Proportion of Variance",2]*100,digits = 1),"% of variance",sep=""),
       main="PCA of normalized data"
  )
  text(normPCA$rotation[,1],normPCA$rotation[,2],labels = sampNames, cex = 0.7,pos = 3)
  outPlts[[cntPlts]] <- recordPlot()
  cntPlts = cntPlts + 1
  graphics.off()
  
  pdf("qc_plots.pdf",onefile = T)
  for (myPlt in outPlts){
    replayPlot(myPlt)
  }
  dev.off()
}