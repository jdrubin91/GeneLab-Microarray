__author__ = 'Jonathan Rubin'

import os, config, metadata_process, rawdata_process

def run(batch_file):
    batch_list = list()
    with open(batch_file) as F:
        parent_dir = F.readline().strip('\n').split('=')[1]
        header = F.readline()
        for line in F:
            linelist = line.strip('\n').split()
            batch_list.append(linelist)


    for i in range(len(batch_list)):
        if 'False' in batch_list[i]:
            GLDS, copy, norm_qc = batch_list[i]

            #Copy module, copies and unzips both metadata and raw data. If precise directories are not found,
            #that GLDS is skipped.
            if copy == 'False':
                print "Copying files for " + GLDS + "..."
                #Process metadata
                metadata_dir = os.path.join(parent_dir,GLDS,'metadata')
                if os.path.isdir(metadata_dir):
                    metadata_process.clean(metadata_dir)
                else:
                    print "metadata directory within " + GLDS + " not found, skipping..."
                    copy, norm_qc = ['Skipped' for j in range(2)]
                    batch_list[i] = [GLDS, copy, norm_qc]

                #Copy rawdata into output
                rawdata_dir = os.path.join(parent_dir,GLDS,'microarray')
                if os.path.isdir(rawdata_dir):
                    rawdata_process.copy(rawdata_dir)
                    rawdata_process.rename(os.path.join(config.outdir,GLDS))
                else:
                    print "microarray directory within " + GLDS + " not found, skipping..."
                    copy, norm_qc = ['Skipped' for j in range(2)]
                    batch_list[i] = [GLDS, copy, norm_qc]

                batch_list[i][1] = 'True'
                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"

            if norm_qc == 'False':
                print "Performing QC, normalization, and post-normalization QC on data for " + GLDS + "..."
                rawdata_out = os.path.join(config.outdir,GLDS,'microarray')
                rawdata_process.qc_and_normalize(rawdata_out,GLDS)
                batch_list[i][2] = 'True'
                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"


def update_batch(parent_dir,header,batch_file,batch_list):
    outfile=open(batch_file,'w')
    outfile.write('#Directory='+parent_dir+'\n')
    outfile.write(header)
    for linelist in batch_list:
        outfile.write('\t'.join(linelist)+'\n')
    outfile.close()