__author__ = 'Jonathan Rubin'

import os, sys, subprocess, config, metadata_process
from difflib import SequenceMatcher

#Finds similarity between two strings a, b. Not currently used but here in case it's needed (JDR 7/5/18)
def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

#Copies most recent metadata and all raw data and then unzips, and removes compressed files
def copy(rawdata_directory):
    #Find name of GLDS number
    GLDS = os.path.basename(os.path.dirname(rawdata_directory))
    rawdata_out = os.path.join(config.outdir,GLDS,'microarray')

    #Make appropriate output directory
    if not os.path.exists(rawdata_out):
        os.makedirs(rawdata_out)

    #Now search the microarray folder for raw data files (this part could be done in a smarter way...)
    for file1 in os.listdir(rawdata_directory):
        if 'raw' in file1 or 'RAW' in file1 or 'Raw' in file1 or 'CEL' in file1 or not 'processed' in file1:
            out_file_path = os.path.join(rawdata_out,file1)

            #Command for copying the raw files to desired output
            cp_command = ["cp", os.path.join(rawdata_directory,file1),out_file_path]

            #md5sum command to check original files
            config.get_md5sum(cp_command[1],'original',action='copy')

            #Then execute the copy command to copy raw files to output directory
            subprocess.call(cp_command)

            #md5sum command to check copied files
            config.get_md5sum(cp_command[2],'new')

            #Once copied, unzip/untar/gunzip compressed directories (if there are any)
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

    #Sometimes compressed files spit out more compressed files so loop through the files once again to catch those and uncompress them
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

#This function renames all raw data files in a specified GLDS_path
def rename(GLDS_path):
    #First get all the proper paths according to specifications
    metadata_out = os.path.join(GLDS_path,'metadata')
    rawdata_out = os.path.join(GLDS_path,'microarray')
    assay_dict = metadata_process.read_assay(metadata_out)
    GLDS = os.path.basename(GLDS_path)
    final_rawdata_out = os.path.join(rawdata_out,'raw')

    #Make the 'raw' directory if it doesn't already exist
    if not os.path.isdir(final_rawdata_out):
        os.makedirs(final_rawdata_out)

    #Loop through the raw data files 
    for filename in os.listdir(rawdata_out):
        if not os.path.isdir(os.path.join(rawdata_out,filename)):

            #Boolean to detect whether the first column corresponds well to the filenames
            sample_in_first_column = False
            sample_in_other_column = False
            extension = filename.split('.')[-1]
            for key in assay_dict:

                #If the first column does correspond well, then assume the first column is the sample name and rename accordingly
                if key in filename:
                    sample_in_first_column = True
                    sample_name = key.replace(' ','-').replace('_','-').replace('(','-').replace(')','-').strip('-')
                    move_command = ["mv", os.path.join(rawdata_out,filename), os.path.join(final_rawdata_out,GLDS+'_'+sample_name+'_microarray_raw.'+extension)]
                    # config.get_md5sum(os.path.join(rawdata_out,filename),'original',action='rename')
                    new_md5sum_file = os.path.join(final_rawdata_out,GLDS+'_'+sample_name+'_microarray_raw.'+extension)

            #If the first column does not correspond to any filenames, simply remove '_', '(', and ')' characters from filename and append appropriate naming conventions
            if not sample_in_first_column:
                for key in assay_dict:
                    for item in assay_dict[key]:
                        if item in filename and item != '':
                            sample_in_other_column = True
                            new_filename = filename.split('.')[0].replace('_','-').replace('(','-').replace(')','-').replace(' ','-').strip('-')
                            move_command = ["mv", os.path.join(rawdata_out,filename), os.path.join(final_rawdata_out,GLDS+'_'+new_filename+'_microarray_raw.'+extension)]
                            # config.get_md5sum(os.path.join(rawdata_out,filename),'original',action='rename')
                            new_md5sum_file = os.path.join(final_rawdata_out,GLDS+'_'+new_filename+'_microarray_raw.'+extension)

            #Execute the command if the file was in metadata - catch whether the file already exists and don't output an error
            if sample_in_first_column or sample_in_other_column:
                try:
                    config.get_md5sum(move_command[1],'original',action='rename')
                    with open(os.devnull,'w') as FNULL:
                        subprocess.check_call(move_command,stdout=FNULL, stderr=subprocess.STDOUT)
                    config.get_md5sum(new_md5sum_file,'new')
                except subprocess.CalledProcessError:
                    config.md5sum['new'].append(('Move Error','N/A'))
            elif not os.path.isdir(os.path.join(rawdata_out,filename)):
                remove_command = ["rm",os.path.join(rawdata_out,filename)]
                config.get_md5sum(os.path.join(rawdata_out,filename),'original',action='remove')
                subprocess.call(remove_command)
                config.md5sum['new'].append(('Removed','N/A'))

#A function to simply detect the array type. Outputs into the 'QC_reporting' directory. Assumes the input is already within 'raw'. This means that the rename
#function has already been called.
def detect_array(GLDS_path):
    GLDS = os.path.basename(GLDS_path)
    rawdata_out = os.path.join(GLDS_path,'microarray')
    R_script = os.path.join(config.R_dir,'affyNormQC.R')
    R_command = ["Rscript", R_script, 
                    "--arrayInfoOnly=TRUE",
                    "--outDir="+rawdata_out, 
                    "--QCDir="+os.path.join(rawdata_out,'QC_reporting'), 
                    "-i", os.path.join(rawdata_out,'raw'),
                    "--GLDS="+GLDS.split('-')[1]]
    subprocess.call(R_command)

    array_info = os.path.join(rawdata_out,'QC_reporting',GLDS+'_arrayInfo.txt')
    if os.path.exists(array_info):
        with open(array_info) as F:
            array = F.readline().strip('\n')
    else:
        print "Warning: Array detection failed for " + GLDS + " - arrayInfo.txt file missing. Skipping..."
        array = 'Skipped'

    return array

#This function simply runs the R script affyNormQC.R specifying the correct inputs
def qc_and_normalize(rawdata_out,GLDS):
    R_script = os.path.join(config.R_dir,'affyNormQC.R')
    R_command = ["Rscript", R_script, 
                    "--normalization", "rma", 
                    "-o", GLDS+"_microarray_normalized",
                    "--outDir="+rawdata_out, 
                    "--QCDir="+os.path.join(rawdata_out,'QC_reporting'), 
                    "-i", os.path.join(rawdata_out,'raw'),
                    "--outType=txt", 
                    "--outputData=TRUE",
                    "--QCoutput=TRUE", 
                    "--NUSEplot=FALSE", 
                    "--GLDS="+GLDS.split('-')[1]]
    subprocess.call(R_command)
    if not GLDS+'_microarray_normalized.txt' in os.listdir(rawdata_out):
        print "Warning: Normalized expression file missing, some processing steps may have failed"

def annotate(rawdata_out,GLDS):
    R_script = os.path.join(config.R_dir,'annotateProbes.R')
    normalized_expression = os.path.join(rawdata_out,GLDS+"_microarray_normalized.txt")
    array_info = os.path.join(rawdata_out,'QC_reporting',GLDS+"_arrayInfo.txt")
    R_command = ["Rscript", R_script, 
                    "-i", normalized_expression,
                    "-a", array_info, 
                    "-o", os.path.join(rawdata_out,GLDS+"_microarray_normalized-annotated.txt"),
                    "--GLDS="+GLDS.split('-')[1],
                    "--QCDir="+os.path.join(rawdata_out,'QC_reporting')]
    subprocess.call(R_command)

    
#This function runs limma differential script limmaDiffExp.R
def limma_differential(rawdata_out,metadata_out,GLDS):
    condition1,condition2,pval_cut = config.visualize.split(',')
    limma_script = os.path.join(config.R_dir,'limmaDiffExp.R')
    if GLDS + "_microarray_normalized_annotated.txt" in os.listdir(rawdata_out):
        d_path = os.path.join(rawdata_out,GLDS + "_microarray_normalized-annotated.txt")
    elif GLDS + "_microarray_normalized.txt" in os.listdir(rawdata_out):
        print "Warning: No annotated expression file detected, using unannotated file instead - Labelling will be done with probe IDs and will need to be converted to gene names manually."
        d_path = os.path.join(rawdata_out,GLDS + "_microarray_normalized.txt")
    else:
        print "Error: No expression count file detected, exiting..."
        sys.exit(1)
    limma_differential_command = ["Rscript", "--vanilla", limma_script, 
                                    "-d", d_path, 
                                    "-i", metadata_out, 
                                    "--group1=" + condition1, 
                                    "--group2=" + condition2, 
                                    "-o", os.path.join(rawdata_out,GLDS + "_microarray_DGE.txt")]
    subprocess.call(limma_differential_command)



