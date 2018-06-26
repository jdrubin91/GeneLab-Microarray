__author__ = 'Jonathan Rubin'

import os, sys, subprocess, config, metadata_process

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
            cp_command = ["cp", os.path.join(rawdata_directory,file1),out_file_path]
            subprocess.call(cp_command)
            if '.zip' in file1:
                unzip_command = ["unzip", "-o", "-qq",out_file_path,"-d",rawdata_out]
                subprocess.call(unzip_command)
                remove_command = ["rm", out_file_path]
                subprocess.call(remove_command)
            if '.gz' in file1:
                gunzip_command = ["gunzip", "-f",out_file_path]
                subprocess.call(gunzip_command)
            if '.tar' in file1:
                untar_command = ["tar", "-xf", out_file_path, "-C", rawdata_out]
                subprocess.call(untar_command)
                remove_command = ["rm", out_file_path]
                subprocess.call(remove_command)

    for file2 in os.listdir(rawdata_out):
        out_file_path = os.path.join(rawdata_out,file2)
        if '.zip' in file2:
            unzip_command = ["unzip", "-o", "-qq", out_file_path, "-d", rawdata_out]
            subprocess.call(unzip_command)
            remove_command = ["rm", out_file_path]
            subprocess.call(remove_command)
        if '.gz' in file2:
            gunzip_command = ["gunzip", "-f", out_file_path]
            subprocess.call(gunzip_command)
        if '.tar' in file2:
            untar_command = ["tar", "-xf", out_file_path, "-C", rawdata_out]
            subprocess.call(untar_command)
            remove_command = ["rm", out_file_path]
            subprocess.call(remove_command)

    
def detect_array(GLDS_path):
    try:
        os.chdir(rawdata_out)
        R_script = os.path.join(config.R_dir,'affyNormQC.R')
        R_command = ["Rscript", R_script, "--arrayInfoOnly=TRUE"]
        subprocess.call(R_command)
        if not 'exprsValues.txt' in os.listdir(rawdata_out):
            print "Warning: Normalized expression file missing, some processing steps may have failed"
        os.chdir(config.srcdir)
    except OSError:
        print "Error: Microarray raw data directory missing. Exiting..."
        sys.exit(1)

def rename(GLDS_path):
    metadata_out = os.path.join(GLDS_path,'metadata')
    rawdata_out = os.path.join(GLDS_path,'microarray')
    assay_dict = metadata_process.read_assay(metadata_out)
    GLDS = os.path.basename(GLDS_path)
    for key in assay_dict:
        for filename in os.listdir(rawdata_out):
            if key in filename:
                extension = filename.split('.')[-1]
                move_command = ["mv", os.path.join(rawdata_out,filename), os.path.join(rawdata_out,GLDS+'_microarray_'+key+'.'+extension)]
                try:
                    with open(os.devnull,'w') as FNULL:
                        subprocess.check_call(move_command,stdout=FNULL, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError:
                    pass


def qc_and_normalize(rawdata_out,GLDS):
    R_script = os.path.join(config.R_dir,'affyNormQC.R')
    R_command = ["Rscript", R_script, "--normalization", "rma", "-o", GLDS+"_microarray_normalized", "-i", rawdata_out,"--outType=txt", "--outputData=TRUE" ,"--QCoutput=TRUE", "--NUSEplot=FALSE", "--GLDS="+GLDS.split('-')[1]]
    subprocess.call(R_command)
    if not GLDS+'_microarray_normalized.txt' in os.listdir(rawdata_out):
        print "Warning: Normalized expression file missing, some processing steps may have failed"
    
def limma_differential(rawdata_out,metadata_out,GLDS):
    condition1,condition2,pval_cut = config.visualize.split(',')
    limma_script = os.path.join(config.R_dir,'limmaDiffExp.R')
    if GLDS + "_microarray_normalized_annotated.txt" in os.listdir(rawdata_out):
        d_path = os.path.join(rawdata_out,GLDS + "_microarray_normalized_annotated.txt")
    elif GLDS + "_microarray_normalized.txt" in os.listdir(rawdata_out):
        print "Warning: No annotated expression file detected, using unannotated file instead - Labelling will be done with probe IDs and will need to be converted to gene names manually."
        d_path = os.path.join(rawdata_out,GLDS + "_microarray_normalized.txt")
    else:
        print "Error: No expression count file detected, exiting..."
        sys.exit(1)
    group1_option = "--group1=" + condition1
    group2_option = "--group2=" + condition2
    limma_differential_command = ["Rscript", "--vanilla", limma_script, "-d", d_path, "-i", metadata_out, group1_option, group2_option, "-o", os.path.join(rawdata_out,GLDS + "_microarray_DGE.txt")]
    subprocess.call(limma_differential_command)



