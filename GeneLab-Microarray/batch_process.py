__author__ = 'Jonathan Rubin'

import config, metadata_process, rawdata_process

def run(batch_file):
    batch_list = list()
    with open(batch_file) as F:
        parent_dir = F.readline().strip('\n').split('=')[1]
        header = F.readline()
        for line in F:
            linelist = line.strip('\n').split('\t')
            batch_list.append(linelist)

    new_batch = open(batch_file,'w')
    for i in range(len(batch_list)):
        if 'False' in batch_list[i]:
            GLDS, copy, chip, qc, norm, qc_norm = batch_list[i]

            #Copy module, copies and unzips both metadata and raw data. If precise directories are not found,
            #that GLDS is skipped.
            if copy == 'False':
                #Process metadata
                metadata_dir = os.path.join(config.indir,'metadata')
                if os.path.isdir(metadata_dir):
                    metadata_process.clean(metadata_dir)
                else:
                    print "metadata directory within " + GLDS + " not found, skipping..."
                    GLDS, copy, chip, qc, norm, qc_norm = ['Skipped' for i in range(6)]
                    batch_list[i] = ['Skipped' for i in range(6)]

                #Copy rawdata into output
                rawdata_dir = os.path.join(config.indir,'microarray')
                if os.path.isdir(rawdata_dir):
                    rawdata_process.copy(rawdata_dir)
                else:
                    print "microarray directory within " + GLDS + " not found, skipping..."
                    GLDS, copy, chip, qc, norm, qc_norm = ['Skipped' for i in range(6)]
                    batch_list[i] = ['Skipped' for i in range(6)]

                batch_list[i][1] = 'True'
                update_batch(parent_dir,header,batch_file,batch_list)


def update_batch(parent_dir,header,batch_file,batch_list):
    outfile=open(batch_file,'w')
    for linelist in batch_list:
        outfile.write('\t'.join(linelist)+'\n')
    outfile.close()