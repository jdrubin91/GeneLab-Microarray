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
            GLDS, copy, chip, qc, norm, norm_qc = batch_list[i]

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
                    copy, chip, qc, norm, norm_qc = ['Skipped' for j in range(5)]
                    batch_list[i] = [GLDS, copy, chip, qc, norm, norm_qc]

                #Copy rawdata into output
                rawdata_dir = os.path.join(parent_dir,GLDS,'microarray')
                if os.path.isdir(rawdata_dir):
                    rawdata_process.copy(rawdata_dir)
                    rawdata_process.rename(os.path.join(config.outdir,GLDS))
                else:
                    print "microarray directory within " + GLDS + " not found, skipping..."
                    copy, chip, qc, norm, norm_qc = ['Skipped' for j in range(5)]
                    batch_list[i] = [GLDS, copy, chip, qc, norm, norm_qc]

                batch_list[i][1] = 'True'
                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"

            if chip == 'False':
                print "Detecting array type for " + GLDS + "..."
                batch_list[i][2] = 'True'
                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"

            if qc == 'False':
                print "Performing initial QC for " + GLDS + "..."
                batch_list[i][3] = 'True'
                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"

            if norm == 'False':
                print "Normalizing data for " + GLDS + "..."
                batch_list[i][4] = 'True'
                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"

            if norm_qc == 'False':
                print "Performing QC on normalized data for " + GLDS + "..."
                batch_list[i][5] = 'True'
                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"


def update_batch(parent_dir,header,batch_file,batch_list):
    outfile=open(batch_file,'w')
    outfile.write('#Directory='+parent_dir+'\n')
    outfile.write(header)
    for linelist in batch_list:
        outfile.write('\t'.join(linelist)+'\n')
    outfile.close()