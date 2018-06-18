__author__ = 'Jonathan Rubin'

from stat import ST_MTIME
import os, sys, time, config

def clean(metadata_directory):
    #Path to the directory (absolute)
    dirpath = metadata_directory

    #Get all entries in the directory w/ stats
    entries = (os.path.join(dirpath, fn) for fn in os.listdir(dirpath))
    entries = ((os.stat(path), path) for path in entries)

    #Insert creation date so we can get the last modified metadata
    entries = ((stat[ST_MTIME], path) for stat, path in entries)

    #Find name of GLDS number
    GLDS = os.path.basename(os.path.dirname(metadata_directory))
    #metadata_out is the path to the output metadata
    metadata_out = os.path.join(config.outdir,GLDS,'metadata')

    #Make appropriate output directory
    if not os.path.exists(metadata_out):
        os.makedirs(metadata_out)

    #Get last modified metadata zip file, copy to the output directory, unzip it, remove the zipped directory, and finally bring all files within folders to the top metadata directory
    i = 0
    for cdate, path in sorted(entries,reverse=True):
        if 'zip' in path and i == 0:
            metadata_zip = os.path.join(metadata_directory,os.path.basename(path))

            #Copy the last modified metadata
            cp_command = "cp -r " + metadata_zip + " " + metadata_out
            #Unzip it into the metadata_out directory
            unzip_command = "unzip -o -qq " + os.path.join(metadata_out,os.path.basename(metadata_zip)) + " -d " + metadata_out
            #Remove the .zip compressed file to avoid confusion and save space
            remove_zip_command = "rm " + os.path.join(metadata_out,os.path.basename(metadata_zip))

            #Execute commands
            os.system(cp_command)
            os.system(unzip_command)
            os.system(remove_zip_command)

            i += 1

    #Loop through the metadata_out directory in case the unzipping produces a folder. If so, mv contents of folder up one directory and remove folder
    for filename in os.listdir(metadata_out):
        if os.path.isdir(os.path.join(metadata_out,filename)):
            move_command = "mv " + os.path.join(metadata_out,filename,"*") + " " + metadata_out
            remove_folder_command = "rm -r " + os.path.join(metadata_out,filename)
            os.system(move_command)
            os.system(remove_folder_command)

def read_assay(metadata_out):
    #Loop through metadata files, find the assay file (starts with 'a_')
    for filename in os.listdir(metadata_out):
        if 'a_' in filename[:2]:
            assay_file = os.path.join(metadata_out,filename)

    #Create an assay dictionary where the key is the name of the sample file
    assay_dict = dict()
    try:
        with open(assay_file) as F:
            F.readline()
            for line in F:
                linelist = [item.strip('"') for item in line.strip('\n').split()]
                assay_dict[linelist[0]] = linelist[1:]

        return assay_dict
    except:
        print "File Error: No assay file found in ISA metadata. Exiting..."
        sys.exit(1)
