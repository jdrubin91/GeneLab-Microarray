# GeneLab-Microarray

## Table of Contents
1. <A href=#Installation>Installation</A>


<H2 id="Installation">Installation</H2>
To install GeneLab-Microarray, follow the steps below to clone this repository and add it to your path:

```
git clone https://github.com/jdrubin91/GeneLab-Microarray.git
cd GeneLab-Microarray/
pip install -e .
```

* Note: You should be in the topmost GeneLab-Microarray directory (not the one that contains .py scripts)


Once the above steps are completed without error, you should be able to call GeneLab-Microarray from any directory. Try:

```
GeneLab-Microarray --help
```

## Directory structure
GeneLab-Microarray expects directories to be in a specific structure. A parent directory with GLDS-# followed by two subdirectories (where one is named metadata and the other microarray) each of which contain zipped archives with either raw microarray data or ISA formatted metadata. For example:

```
GLDS-#/
  metadata/
    Metadata_ISA.zip
  microarray/
    microarray_raw.tar
```


## Batch file format
If `-b,--batch` option is desired. In addition to calling the flag, users must submit the full path to a batch.txt file (examples and a simple script to create this batch.txt file is located in the batch subdirectory). Briefly, the batch.txt file expects the first line to begin with '#' followed by 'Directory=' then a full path to a directory. For example:

`#Directory=/opt/genelab-genomespace-dev_mount_point/`

The rest of the file is a tab delimited txt file with 6 columns:

`GLDS# Copied  Chip  1.QC/QA 2.Normalization 3.Normalized_QC/QA`

The first column is the name of a folder within the specified Directory. All subsequent columns are booleans (True or False) and are used to keep track of the progress of processing the desired data in batch. GeneLab-Microarray will overwrite the specified batch.txt file changing booleans to True or Skipped when the specific step is finished.
