# GeneLab-Microarray

## Table of Contents
1. <A href=#Installation>Installation</A>
2. <A href=#Running>Running GeneLab-Microarray</A>
3. <A href=#Directory>Directory Structure</A>
4. <A href=#BatchFile>Batch File Format</A>
5. <A href=#affyNormQC>Affy QC and Normalization</A>


<H2 id="Installation">Installation</H2>
To install GeneLab-Microarray, follow the steps below to clone this repository and add it to your path:

```
git clone https://github.com/jdrubin91/GeneLab-Microarray.git
cd GeneLab-Microarray/
pip install -e .
```

* Note: You should be in the topmost GeneLab-Microarray directory when pip installing (not the one that contains .py scripts)

If you are installing on the DP server, run this pip command instead:

```
pip --user install -e .
```

Once the above steps are completed without error, you should be able to call GeneLab-Microarray from any directory. Try:

```
GeneLab-Microarray --help
```

<H2 id="Running">Running GeneLab-Microarray</H2>
Gene-Lab microarray has two modes: process and visualization. 

<H3 id="ProcessMode">Process Mode</H3>
The default for this package is processing mode which requires as input a path to a GLDS directory. GeneLab-Microarray will copy raw files from your specified directory into your output directory (specified with the -o/--output flag - Note: the output option is required). 

Example:
```
GeneLab-Microarray --output /opt/jdrubin/batch_out/ /opt/genelab-genomespace-dev_mount_point/GLDS-4/
```
<H4 id="BatchSubmode">Batch Submode</H4>
Within the processing mode, there is a submode called batch (specified with -b/--batch). Similar to process mode, an output directory must be specified. However, instead of specifying a directory, batch mode takes as input a path to a batch file (see <A href=#BatchFile>batch file format</A> for guidelines on batch files. There are also example batch files within the batch/ directory).

Example:
```
GeneLab-Microarray --output /opt/jdrubin/batch_out/ --batch /opt/jdrubin/GeneLab-Microarray/batch/batch.txt
```

<H3 id="VisualizationMode">Visualization Mode</H3>
The visualization mode for GeneLab-Microarray is specified with the -v/--visualize flag and takes as input a comma separated list of factor values followed by an adjusted p-value cutoff (for plotting purposes). The visualization mode will output an html file with interactive graphs and a png file with identical graphs that can be used for publication. Similar to other modes, an output must be specified. As input, the visualization mode takes a path to a GLDS folder that has already been processed by GeneLab-Microarray.

Example:
```
GeneLab-Microarray --output /opt/jdrubin/batch_out/ --visualize flight,ground,0.1 /opt/jdrubin/batch_out/GLDS-4/
```

The visualization mode can also do comparisons with multiple factor values which are specified with underscore delimiters.
Example:
```
GeneLab-Microarray --output /opt/jdrubin/batch_out/ --visualize flight_geneKO,flight_noKO,0.1 /opt/jdrubin/batch_out/GLDS-4/
```

<H2 id="Directory">Directory Structure</H2>
GeneLab-Microarray expects directories to be in a specific structure. A parent directory with GLDS-# followed by two subdirectories (where one is named metadata and the other microarray) each of which contain zipped archives with either raw microarray data or ISA formatted metadata. For example:

```
GLDS-#/
  metadata/
    Metadata_ISA.zip
  microarray/
    microarray_raw.tar
```


<H2 id="BatchFile">Batch File Format</H2>
If `-b,--batch` option is desired. In addition to calling the flag, users must submit the full path to a batch.txt file (examples and a simple script to create this batch.txt file is located in the batch subdirectory). Briefly, the batch.txt file expects the first line to begin with '#' followed by 'Directory=' then a full path to a directory. For example:

```
#Directory=/opt/genelab-genomespace-dev_mount_point/
```

The rest of the file is a tab delimited txt file with 3 columns:

```
GLDS#     Copied    Normalize/QC
GLDS-4    False     False
```

The first column is the name of a folder within the specified Directory. All subsequent columns are booleans (True or False) and are used to keep track of the progress of processing the desired data in batch. GeneLab-Microarray will overwrite the specified batch.txt file changing booleans to True or Skipped when the specific step is finished. An example of a batch.txt file can be found within the `batch/` folder

<H2 id="affyNormQC">Affy QC and Normalization</H2>
This script is to be called from the directory containing Affymetrix microarray files (with a .CEL extension). It can determine the version of the array and load the appropriate packages (ie "affy" for earlier microarrays and "oligo" for the newer arrays). No inputs are required to run it, but to view the available options, simply run the line below:

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
An example run with specified options for all available parameters is given below:

```
Rscript --vanilla affyNormQC.R -n rma --outFile=GLDS-n_FLT-GC_microarray_exprsValues --outType=txt --outputData=TRUE --QCoutput=TRUE --NUSEplot=TRUE
```

