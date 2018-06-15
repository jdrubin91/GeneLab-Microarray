__author__ = 'Jonathan Rubin'

import os, time, subprocess, config

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
            rsync_command = "rsync " + os.path.join(rawdata_directory,file1) + " " + out_file_path
            os.system(rsync_command)
            if '.zip' in file1:
                unzip_command = "unzip -o " + out_file_path + " -d " + rawdata_out
                os.system(unzip_command)
            if '.gz' in file1:
                gunzip_command = "gunzip -f " + out_file_path
                os.system(gunzip_command)
            if '.tar' in file1:
                untar_command = "tar -xf " + out_file_path + " -C " + rawdata_out
                os.system(untar_command)

    for file2 in os.listdir(rawdata_out):
        out_file_path = os.path.join(rawdata_out,file2)
        if '.zip' in file2:
            unzip_command = "unzip -o " + out_file_path + " -d " + rawdata_out
            os.system(unzip_command)
        if '.gz' in file2:
            gunzip_command = "gunzip -f " + out_file_path
            print gunzip_command
            os.system(gunzip_command)
            time.sleep(1)
        if '.tar' in file2:
            untar_command = "tar -xf " + out_file_path + " -C " + rawdata_out
            os.system(untar_command)

    