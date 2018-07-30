indir = "./genelab-genomespace-dev_mount_point-local/GLDS-4/"
outdir = "test"
srcdir = "/Users/jonathanrubin/Google_Drive/NASA/home/GeneLab-Microarray/GeneLab-Microarray"
wrkdir = "/Users/jonathanrubin/Google_Drive/NASA/home"
tempdir = "/Users/jonathanrubin/Google_Drive/NASA/home/GeneLab-Microarray/temp"
R_dir = "/Users/jonathanrubin/Google_Drive/NASA/home/GeneLab-Microarray/GeneLab-Microarray/R_scripts"
md5sum = {"original": [], "new": []}
batch = "False"
visualize = "False"
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
