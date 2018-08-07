indir = "False"
outdir = "../priority_batch/GLDS-111"
srcdir = "/Users/jonathanrubin/Google_Drive/NASA/home/GeneLab-Microarray/GeneLab-Microarray"
wrkdir = "/Users/jonathanrubin/Google_Drive/NASA/home/GeneLab-Microarray"
tempdir = "/Users/jonathanrubin/Google_Drive/NASA/home/GeneLab-Microarray/temp"
R_dir = "/Users/jonathanrubin/Google_Drive/NASA/home/GeneLab-Microarray/GeneLab-Microarray/R_scripts"
md5sum = {"original": [], "new": []}
batch = "False"
visualize = "Spaceflight,Ground Control,0.1"
output = str()
microarray_out="microarray"
def get_md5sum(filepath,key,action=False):
    import os, subprocess
    if not action:
        md5sum_command = ["md5sum",filepath]
        md5sum_character = subprocess.check_output(md5sum_command).split(' ')[0].encode("utf-8")
        md5sum[key].append((os.path.basename(os.path.normpath(filepath)),md5sum_character))
    else:
        md5sum_command = ["md5sum",filepath]
        md5sum_character = subprocess.check_output(md5sum_command).split(' ')[0].encode("utf-8")
        md5sum[key].append((action,os.path.basename(os.path.normpath(filepath)),md5sum_character))
def write_output(output_file,output):
    with open(output_file,'a') as F:
        F.write(output)
def detect_2channel(infile):
    is2channel = False
    with open(infile,'r') as F:
        for line in F:
            if 'rMedianSignal' in line or 'rBGMedianSignal' in line:
                is2channel = True
                return is2channel
    return is2channel
