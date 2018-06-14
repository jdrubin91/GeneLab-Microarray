rm(list=ls())

# source("http://bioconductor.org/biocLite.R")
# biocLite("affy")
# biocLite("mogene10sttranscriptcluster.db")
# biocLite("Risa")


library(affy)
library(limma)
library(mogene10sttranscriptcluster.db)
library(annotate)
library(Risa)

color = sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1:-7]) # Read in a randomly ordered library of colors for plotting

# ls("package:mogene10sttranscriptcluster.db") # List of R objects in the package
# mogene10sttranscriptcluster() # QC info

raw = ReadAffy() # Read in all of the raw .CEL files in the current directory

# Possible raw QC step, plotting the intensity distributions of the pm probes from each microarray on the same graph
mypms = pm(raw)
ylims = c(0,1.5)
xlims = c(0,5)
for(i in 1:ncol(mypms)){
  if(i == 1){
    plot(density(log10(mypms[,i])),ylim = ylims,xlim=xlims,xlab='log10(Raw Intensities)',main='All intensity distributions',col=color[i])
    par(new=T)
  }else{
    plot(density(log10(mypms[,i])),ylim = ylims,xlim=xlims,axes=T,xlab='',main='',col=color[i])
    par(new=T)
  }
  print(i)
}
par(new=F)

# Background correct, normalize, calculate expression values
eset <- expresso(raw, normalize.method="quantiles",bgcorrect.method="rma",pmcorrect.method="pmonly",summary.method="liwong")

# Post-normalization QC step
ylims = c(0,1.5)
xlims = c(0,5)
normVals = exprs(eset)
for(i in 1:ncol(normVals)){
  if(i == 1){
    plot(density(log10(normVals[,i])),ylim = ylims,xlim=xlims,xlab='log10(normalized expression values)',main='All expression distributions',col=color[i])
    par(new=T)
  }else{
    plot(density(log10(normVals[,i])),ylim = ylims,xlim=xlims,axes=T,xlab='',main='',col=color[i])
    par(new=T)
  }
  print(i)
}
par(new=F)

# Mapping Affy transcript cluster IDs to RefSeq names from the imported library
ID = featureNames(eset) # Pulls out the AffyIDs
mapFun = function(id){ # Function to match the primary RefSeq ID for a given AffyID and return NA in all other cases
  return(tryCatch(get(id, env=mogene10sttranscriptclusterREFSEQ)[1], error=function(e) NA))
}
RefSeq = lapply(ID,FUN = mapFun) # Applying mapFun to all AffyIDs

# Replace AffyIDs with RefSeq IDs, drop probes w/o RefSeq IDs?

write.exprs()

counts <- read.delim('"""+directory+"""affy_norm_annot.txt',header=T,sep="\t")
rownames(counts) <- counts$'Gene_Probe'
counts <- subset(counts,select=c("""+','.join([str(i) for i in range(2,samplenumber)])+"""))
summary(counts)
head(counts,10)
png('"""+directory+"""MDS.png')
plotMDS(counts)


