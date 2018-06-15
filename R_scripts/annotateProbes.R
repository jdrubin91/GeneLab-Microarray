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
