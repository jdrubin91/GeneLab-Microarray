__author__ = 'Jonathan Rubin'

import os, config, metadata_process, rawdata_process

def run(batch_file):
    compatible_arrays = ['Affymetrix']
    batch_list = list()
    with open(batch_file) as F:
        parent_dir = F.readline().strip('\n').split('=')[1]
        header = F.readline()
        for line in F:
            linelist = line.strip('\n').split('\t')
            if len(linelist) != 5:
                print "Error, batch file line not formatted properly: " + line + " skipping..."
            else:
                batch_list.append(linelist)


    for i in range(len(batch_list)):
        if 'False' in batch_list[i]:
            GLDS, copy, array, norm_qc, annotate = batch_list[i]
            GLDS_path = os.path.join(config.outdir,GLDS)
            rawdata_out = os.path.join(config.outdir,GLDS,'microarray')
            metadata_out = rawdata_out = os.path.join(config.outdir,GLDS,'metadata')

            #Copy module, copies and unzips both metadata and raw data. If precise directories are not found,
            #that GLDS is skipped.
            if copy == 'False':
                print "Copying files for " + GLDS + "..."
                #Process metadata
                metadata_in = os.path.join(parent_dir,GLDS,'metadata')
                if os.path.isdir(metadata_in):
                    metadata_process.clean(metadata_in)
                else:
                    print "metadata directory within " + GLDS + " not found, skipping..."
                    copy, array, norm_qc, annotate = ['Skipped' for j in range(4)]
                    batch_list[i] = [GLDS, copy, array, norm_qc, annotate]

                #Copy rawdata into output
                rawdata_in = os.path.join(parent_dir,GLDS,'microarray')
                if os.path.isdir(rawdata_in):
                    rawdata_process.copy(rawdata_in)
                    rawdata_process.rename(os.path.join(config.outdir,GLDS))
                else:
                    print "microarray directory within " + GLDS + " not found, skipping..."
                    copy, array, norm_qc, annotate = ['Skipped' for j in range(4)]
                    batch_list[i] = [GLDS, copy, array, norm_qc, annotate]

                metadata_process.create_md5sum_out(rawdata_out,GLDS)
                batch_list[i][1] = 'True'
                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"
            elif copy != 'True':
                print "Warning: Files were not copied for " + GLDS + ". If this was not desired, check batch file and make sure this GLDS was set to 'False'."

            #Array module, this part simply generates an arrayInfo.txt file and reads it in. If the array is part of a list of arrays that we can process then
            #continue, otherwise skip the GLDS
            if array == 'False':
                print "Detecting array type for " + GLDS + "..."
                array = rawdata_process.detect_array(GLDS_path)
                if array != 'Skipped':
                    batch_list[i][2] = array
                    if array in compatible_arrays:
                        update_batch(parent_dir,header,batch_file,batch_list)
                        print "done"
                    else:
                        print "Warning: " + GLDS + " " + array + " arrays not currently supported, skipping..."
                        norm_qc, annotate = ['Skipped' for j in range(2)]
                        batch_list[i] = batch_list[i][:2]+[array,norm_qc,annotate]
                        update_batch(parent_dir,header,batch_file,batch_list)
                else:
                    norm_qc, annotate = ['Skipped' for j in range(2)]
                    batch_list[i] = batch_list[i][:2]+[array,norm_qc,annotate]
                    update_batch(parent_dir,header,batch_file,batch_list)
            elif array != 'True' or array != 'Skipped':
                print "Warning: Array was not detected for " + GLDS + ". If this was not desired, check batch file and make sure this GLDS was set to 'False'."

            #Performs normalization and qc pre- and post-normalization
            if norm_qc == 'False':
                print "Performing QC, normalization, and post-normalization QC on data for " + GLDS + "..."
                rawdata_process.qc_and_normalize(rawdata_out,GLDS)
                batch_list[i][3] = 'True'
                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"
            elif norm_qc != 'True' or norm_qc != 'Skipped':
                print "Warning: QC and normalization not performed for " + GLDS + ". If this was not desired, check batch file and make sure this GLDS was set to 'False'."

            #Annotates probeIDs with gene names. Autodetection of array annotation package is attempted but if it fails then return 'Skipped'.
            if annotate == 'False':
                print "Annotating probe IDs with gene names for " + GLDS + "..."
                rawdata_process.annotate(rawdata_out,GLDS)
                batch_list[i][4] = 'True'
                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"


#This function updates the batch file
def update_batch(parent_dir,header,batch_file,batch_list):
    with open(batch_file,'w') as outfile:
        outfile.write('#Directory='+parent_dir+'\n')
        outfile.write(header)
        for linelist in batch_list:
            outfile.write('\t'.join(linelist)+'\n')
        