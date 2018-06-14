__author__ = 'Jonathan Rubin'

from stat import ST_MTIME
import os, sys, time, config

def clean(metadata_directory):
    #Path to the directory (absolute)
    dirpath = metadata_directory

    #Get all entries in the directory w/ stats
    entries = (os.path.join(dirpath, fn) for fn in os.listdir(dirpath))
    entries = ((os.stat(path), path) for path in entries)

    #Insert creation date
    entries = ((stat[ST_MTIME], path) for stat, path in entries)

    #Find name of GLDS number
    GLDS = os.path.basename(os.path.dirname(metadata_directory))
    metadata_out = os.path.join(config.outdir,GLDS,'metadata')

    #Make appropriate output directory
    if not os.path.exists(metadata_out):
        os.makedirs(metadata_out)

    #Get last modified metadata zip file, copy to the output directory, and unzip it
    i = 0
    for cdate, path in sorted(entries,reverse=True):
        if 'zip' in path and i == 0:
            metadata_zip = os.path.join(metadata_directory,os.path.basename(path))
            command1 = "rsync -r " + metadata_zip + " " + metadata_out
            command2 = "unzip " + os.path.join(metadata_out,os.path.basename(metadata_zip))
            os.system(command1)
            os.system(command2)
            i += 1