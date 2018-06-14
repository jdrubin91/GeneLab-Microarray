# GeneLab-Microarray

## Batch file format
If -b,--batch option is desired. In addition to calling the flag, users must submit the full path to a batch.txt file (examples and a simple script to create this batch.txt file is located in the batch subdirectory). Briefly, the batch.txt file expects the first line to begin with '#' followed by 'Directory=' then a full path to a directory. For example:

`#Directory=/opt/genelab-genomespace-dev_mount_point/`

This directory is expected to be in a specific format:
```
GLDS-#/
metadata/ microarray/
zipped files
```

