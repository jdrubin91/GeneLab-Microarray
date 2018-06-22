#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("Risa")
# biocLite("limma")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list=list(
  make_option(c("-d","--exprData"),type="character",help="Name of (or path to) the input file (\\t delimited .txt file)"),
  make_option(c("-i","--ISApath"),type="character",help="Path to the directory containing the dataset metadata"),
  make_option("--group1",type="character",help="'_'delimited list of factors to select samples for group 1 [ex: flight_geneKO]"),
  make_option("--group2",type="character",help="'_'delimited list of factors to select samples for group 2 [ex: ground_geneKO]"),
  make_option(c("-o","--output"),type="character",default="DGE.txt",help="Name of (or path to) file to write results to (default: DGE.txt)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$exprData)){
  print_help(opt_parser)
  stop("No expression data provided", call.=FALSE)
}else inFH = opt$exprData

if (is.null(opt$ISApath)){
  print_help(opt_parser)
  stop("No ISA directory provided", call.=FALSE)
} else isaFH = opt$ISApath

tryCatch({
  suppressPackageStartupMessages(library("Risa"))
  tabISA = FALSE
}, error=function(e){
  cat("Warning: Risa package could not be loaded. Metadata will be parsed as tab-delimited files.\n")
  tabISA = TRUE
})
suppressPackageStartupMessages(library("limma"))

# Read in ISA tab file and extract assay file
# isaFH = "../metadata/GLDS-4_metadata_GSE18388-ISA/"
if(tabISA==TRUE){
  tryCatch({
    isaFiles = dir(isaFH)
    sFile = isaFiles[grep("^s_*",isaFiles)]
    studyFactors = read.delim(paste(isaFH,sFile,sep=""),header=T, sep="",stringsAsFactors = F)
  }, error=function(e){
    stop("ISA files could not be read by parsing tab-delimited files", call. = F)
  })
  isaFiles = dir(isaFH)
  sFile = isaFiles[grep("^i_*",isaFiles)]
}else{
  tryCatch({
    isaIN = readISAtab(isaFH)
    studyFactors = isaIN@study.files[[1]]
  }, error=function(e){
    stop("ISA directory could not be read by Risa", call. = F)
  })
}


# Read in underscore-delimited factor levels and split into lists
if (!is.null(opt$group1) & !is.null(opt$group2)){
  fact1 = opt$group1
  fact1 = strsplit(fact1, split = "_")[[1]]
  fact2 = opt$group2
  fact2 = strsplit(fact2, split = "_")[[1]]
} else{
  print_help(opt_parser)
  stop("Factor levels not provided or improperly formated", call.=FALSE)
}

#From assay file, extract column containing 'Factor Value'
tryCatch({
  if(tabISA == TRUE){
    factorValues = studyFactors[,grepl("Factor.Value",colnames(studyFactors))]
    rownames(factorValues) = studyFactors[,grepl("Sample.Name",colnames(studyFactors))]
    
  }else{
    factorValues = studyFactors[,grepl("Factor Value",colnames(studyFactors))]
    rownames(factorValues) = studyFactors[,grepl("Sample Name",colnames(studyFactors))]
  }
}, error=function(e){
  stop("Error: Unable to pull sample names from the study level metadata", call. = F)
})

# Read in an expression value txt file
# inFH = "annotExpValues.txt"
eset <- read.delim(inFH, header = TRUE, stringsAsFactors = F)
#Modify matrix so that the first column becomes the row names
#rownames(eset) <- eset[,1]
#eset <- eset[,-1]
if(nrow(factorValues) != ncol(eset)){
  cat("\nWarning: Number of samples in the expression set not equal to the number of samples in the metadata\n")
}

#From the eset matrix, determine which columns correspond to which factor values
esetSampNames <- colnames(eset)
newOrder = rep(0,ncol(eset))
for(i in 1:nrow(factorValues)){ # Reorder the factorValues dataframe to match the order of sample names in the expression set
  newOrder[i] = grep(pattern = rownames(factorValues)[i], x = esetSampNames,ignore.case = T)
}
factorValues = factorValues[newOrder,] 

group <- rep(0,ncol(eset)) # Create list to hold group assignments for all
for(i in 1:nrow(factorValues)){ # Assign each sample to a group [3 = both groups 1 & 2, 0 = neither group]
  if(all(fact1 %in% factorValues[i,]) & all(fact2 %in% factorValues[i,])){
    group[i] = 3
  }else if(all(fact1 %in% factorValues[i,])){
    group[i] = 1
  }else if(all(fact2 %in% factorValues[i,])){
    group[i] = 2
  }
}
# Error handling
if(sum(group == 3) > 0){
  cat("The following samples belonged to both groups and were removed from further analysis:\n",rownames(factorValues)[group == 3],"\n")
}
if(sum(group == 3) == nrow(factorValues)){
  stop("Error: All of the samples belonged to both groups! Exiting.", call.=F)
}
cat(sum(group == 1),"sample(s) found in group 1:\n",rownames(factorValues)[group == 1],"\n")
cat(sum(group == 2),"sample(s) found in group 2:\n",rownames(factorValues)[group == 2],"\n")
if(sum(group == 0) > 0){
  cat("Warning:",sum(group == 0),"sample(s) not found in either group:\n",rownames(factorValues)[group == 0],"\nIf this is not expected, please ensure the provided factor levels match the factor levels in the study-level metadata exactly\n")
}
eset = eset[,!(group == 0 | group == 3)]
group = group[!(group == 0 | group == 3)]

#Create a design matrix based on the ordering of the columns within eset
group = as.factor(group)
design <- model.matrix(~0+group)

#This part of the script is straight from Limma documentation
fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(group1-group2,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Here we write the results to a tab delimited text file that is ordered by adjusted p-value
#coef refers to which column is of interest (1 is log2FC), adjust refers to multiple hypothesis testing method ("BH" = Benjamini & Hochberg)
table <- data.frame(topTable(fit2, coef=1, n=Inf, adjust="BH"))
write.table(table,file=opt$output,sep="\t")
cat("All done! Differential expression information saved to:",opt$output,"\n")

