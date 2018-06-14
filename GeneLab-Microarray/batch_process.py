__author__ = 'Jonathan Rubin'

import config, metadata_process, rawdata_process

def run(batch_file):
    batch_list = list()
    with open(batch_file) as F:
        parent_dir = F.readline().strip('\n').split('=')[1]
        for line in F:
            linelist = line.strip('\n').split('\t')
            batch_list.append(linelist)

    for linelist in batch_list:
        if 'False' in linelist:
            GLDS, copy, chip, qc, norm, qc_norm = linelist
            if not boolean(copy):
                
                #Process metadata
                metadata_dir = os.path.join(config.indir,'metadata')
                if os.path.isdir(metadata_dir):
                    metadata_process.clean(metadata_dir)
                else:
                    raise IOError('metadata directory within input not found. See README for expected directory structure.')

                #Copy rawdata into output
                rawdata_dir = os.path.join(config.indir,'microarray')
                if os.path.isdir(rawdata_dir):
                    rawdata_process.copy(rawdata_dir)
                else:
                    raise IOError('microarray directory within input not found. See README for expected directory structure.')


