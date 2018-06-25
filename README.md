# GeneLab-Microarray

## Table of Contents
1. <A href=#Installation>Installation</A>
2. <A href=#Running>Running GeneLab-Microarray</A>
   - <A href=#ProcessMode>Process Mode</A>
        + <A href=#BatchSubmode>Batch Submode</A>
   - <A href=#VisualizationMode>Visualization Mode</A>
3. <A href=#Directory>Directory Structure</A>
4. <A href=#BatchFile>Batch File Format</A>
5. <A href=#affyNormQC>Affy QC and Normalization</A>
6. <A href=#annotateProbe>Probe annotation</A>
7. <A href=#limmaDiffExp>Limma Differential Expression</A>


<H2 id="Installation">Installation</H2>
GeneLab-Microarray writes to itself so to ensure you can easily run it, it is recommended to pip install with the --user flag. Example:

```
pip install --user GeneLab-Microarray
```

GeneLab-Microarray should install mpld3 and scipy automatically however matplotlib will need to be installed manually if it is not installed. To do this, run the following pip command:

```
pip install --user matplotlib
```

And to finish installation, open an R session and run the following commands:
```
install.packages("optparse")
source("http://bioconductor.org/biocLite.R")
biocLite("zlibbioc")
biocLite("Biobase")
biocLite("affy")
biocLite("S4Vectors")
biocLite("IRanges")
biocLite("XVector")
biocLite("Biostrings")
biocLite("affyPLM")
biocLite("bit")
biocLite("ff")
biocLite("bitops")
biocLite("RCurl")
biocLite("GenomicRanges")
biocLite("matrixStats")
biocLite("Rcpp")
biocLite("bit64")
biocLite("digest")
biocLite("RSQLite")
biocLite("oligo")
biocLite("XML")
biocLite("genefilter")
biocLite("mogene10sttranscriptcluster.db")
biocLite("moe430a.db")
biocLite("drosophila2.db")
biocLite("hgu133plus2.db")
biocLite("ath1121501.db")
biocLite("yeast2.db")
biocLite("hugene10sttranscriptcluster.db")
biocLite("rat2302.db")
biocLite("limma")
biocLite("colorspace")
biocLite("plyr")
biocLite("scales")
biocLite("lazyeval")
biocLite("rlang")
biocLite("tibble")
biocLite("graph")
biocViews("biocViews")

biocLite("mzR") # I think this requires manual installation of ncdf first?
biocLite("MSnbase")

biocLite("RBGL")

biocLite("Risa")
```

Once the above steps are completed without error, you should be able to call GeneLab-Microarray from any directory. Try:

```
GeneLab-Microarray --help
```

<H2 id="Running">Running GeneLab-Microarray</H2>
Gene-Lab microarray has two modes: process and visualization. There are 3 possible flags to give GeneLab-Microarray (--process,--batch, and --visualize). In all cases, an output directory must be specified as the last argument (no flags) - this is a positional argument.

<H3 id="ProcessMode">Process Mode</H3>
Processing mode requires as input a path to a GLDS directory. GeneLab-Microarray will copy raw files from your specified directory into your output directory, rename them according to standard specifications, and perform QC and normalization.


Example:
```
GeneLab-Microarray --process /opt/genelab-genomespace-dev_mount_point/GLDS-4/ /opt/jdrubin/batch_out/
```

<H4 id="BatchSubmode">Batch Submode</H4>
Within the processing mode, there is a submode called batch (specified with -b/--batch). This submode will batch process all GLDS directories specified within a batch.txt file (see <A href=#BatchFile>batch file format</A>). If specified, input into process a batch.txt file instead of a GLDS directory. There are example batch files within the batch/ directory.



Example:
```
GeneLab-Microarray --batch --process /opt/jdrubin/GeneLab-Microarray/batch/batch.txt /opt/jdrubin/batch_out/
```

<H3 id="VisualizationMode">Visualization Mode</H3>
The visualization mode for GeneLab-Microarray is specified with the -v/--visualize flag and takes as input a comma separated list of factor values followed by an adjusted p-value cutoff (for plotting purposes). The visualization mode will output an html file with interactive graphs and a png file with identical graphs that can be used for publication. Similar to other modes, an output directory must be specified as a positional argument at the end (in this case this is also the input directory). Visualize mode assumes the input data has been processed using GeneLab-Microarray (looks for specific filenames within specific directories).



Example:
```
GeneLab-Microarray --visualize flight,ground,0.1 /opt/jdrubin/batch_out/GLDS-4/
```

The visualization mode can also do comparisons with multiple factor values which are specified with underscore delimiters.



Example:
```
GeneLab-Microarray --visualize flight_geneKO,flight_noKO,0.1 /opt/jdrubin/batch_out/GLDS-4/
```

<H2 id="Directory">Directory Structure</H2>
GeneLab-Microarray expects directories to be in a specific structure. A parent directory with GLDS-# followed by two subdirectories (where one is named metadata and the other microarray) each of which contain zipped archives with either raw microarray data or ISA formatted metadata. For example:

```
GLDS-#/
  metadata/
    metadata_ISA.zip
  microarray/
    microarray_raw.tar
```


<H2 id="BatchFile">Batch File Format</H2>
If `-b,--batch` option is desired. In addition to calling the flag, users must submit the full path to a batch.txt file (examples and a simple script to create this batch.txt file is located in the batch subdirectory). 

Briefly, the batch.txt file expects the first line to begin with '#' followed by 'Directory=' then a full path to a directory. The rest of the file is a tab delimited txt file with 3 columns (header is required). The first column is the name of a folder within the specified Directory. All subsequent columns are booleans (True or False) and are used to keep track of the progress of processing the desired data in batch. For example:

```
#Directory=/opt/genelab-genomespace-dev_mount_point/
GLDS#     Copied    Normalize/QC
GLDS-4    False     False
```

GeneLab-Microarray will overwrite the specified batch.txt file changing booleans to True or Skipped when the specific step is finished. An example of a batch.txt file can be found within the `batch/` folder

<H2 id="affyNormQC">Affy QC and Normalization</H2>
The `affyNormQC.R` script is to be called from the directory containing Affymetrix microarray files (with a .CEL extension). It can determine the version of the array and load the appropriate packages (ie "affy" for earlier microarrays and "oligo" for the newer arrays). No inputs are required to run it, but to view the available options, simply run the line below:

```
Rscript --vanilla affyNormQC.R --help
```

Before running this script, it may be necessary to run the commented out lines immediately below the shebang in an R session to be sure all of the necessary packages are installed

```
install.packages("optparse")
source("http://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("affyPLM")
biocLite("oligo")
```
An example run with all of the options explicitly set to the default or example options:

```
Rscript --vanilla affyNormQC.R -i path/to/input/files -n rma -o expValues --outType=txt --outputData=TRUE --arrayInfoOnly=FALSE --QCoutput=TRUE --NUSEplot=FALSE --GLDS=21
```

This script can also be used to detect the Affymetrix array information only, outputting a text file containing the manufacturer and the array version and quitting before normalizing the data or performing QC. This option can be accessed by setting `--arrayInfoOnly=TRUE`. However, the array information txt file will be output in the standard mode as well.

<H2 id="annotateProbe">Probe annotation</H2>
The `annotateProbes.R` script can be used to map probe IDs to RefSeq gene IDs using annotation packages. If an array type has not been seen before, the annotation package will need to be manually loaded into the array:annotation pseudo-dictonary. The available options for the script are viewable by the following command:

```
Rscript --vanilla annotateProbes.R --help
```

The required packages to be install prior to running shown here, as well as in the script immediately below the shebang.

```
install.packages("optparse")
source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
biocLite("mogene10sttranscriptcluster.db")
```

An example call with all of the default/recommended options explicitly defined:

```
Rscript --vanilla annotateProbes.R -i normalizedProbeExpressions.txt -a GLDS-4_arrayInfo.txt -o annotExpValues.txt --QCoutput=TRUE --GLDS=4
```

<H2 id="limmaDiffExp">Limma Differential Expression</H2>
The `limmaDiffExp.R` script can be used to calculate changes in expression between two groups of samples for a given dataset. The available options can be examined by calling:

```
Rscript --vanilla limmaDiffExp.R --help
```

This script requires the packages: `optparse`, `Risa`, and `limma`. If these packages are not installed, they can be by running the following line in an R session:

```
install.packages("optparse")
source("http://bioconductor.org/biocLite.R")
biocLite("Risa")
biocLite("limma")
```

An example call with options set for all parameters is shown below:

```
Rscript --vanilla limmaDiffExp.R -d GLDS-4_microarray_normalized_annotated.txt -i ../metadata/example_RSA_directory --group1=flight_KO --group2=ground_KO -o GLDS-4_microarray_DGE.txt
```

