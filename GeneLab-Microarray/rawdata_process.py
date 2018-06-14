__author__ = 'Jonathan Rubin'

import os, config

def copy(rawdata_directory):
    #Find name of GLDS number
    GLDS = os.path.basename(os.path.dirname(rawdata_directory))
    rawdata_out = os.path.join(config.outdir,GLDS,'microarray')

    #Make appropriate output directory
    if not os.path.exists(rawdata_out):
        os.makedirs(rawdata_out)


    for file1 in os.listdir(rawdata_directory):
        if 'raw' in file1 or 'RAW' in file1 or 'Raw' in file1 or 'CEL' in file1 or not 'processed' in file1:
            out_file_path = os.path.join(rawdata_out,file1)
            rsync_command = "rsync " + os.path.join(rawdata_directory,file1) + " " + os.path.join(rawdata_out,file1)
            os.system(rsync_command)
            if 'zip' in file1:
                unzip_command = "unzip -o " + os.path.join(rawdata_out,file1) + " -d " + rawdata_out
                os.system(unzip_command)
            if 'gz' in file1:
                gunzip_command = "gunzip " + os.path.join(rawdata_out,file1)
                os.system(gunzip_command)
            if 'tar' in file1:
                untar_command = "tar -xvfz " + os.path.join(rawdata_out,file1)
                os.system(untar_command)
    