__author__ = 'Jonathan Rubin'

import os,sys,argparse

def run():
    parser = argparse.ArgumentParser(prog='GeneLab-Microarray',usage='%(prog)s [options] Directory',description='Standardized processing pipeline for microarray data on GeneLab.')
    parser.add_argument('Directory',
        help='The full path to a directory containing either a single dataset structured according to specifications or a batch.txt file. See README for more information.')
    parser.add_argument('-o','--output',help='Required. Takes as input an output directory.',metavar='',required=True)
    parser.add_argument('-b','--batch',help='If batch processing is desired, provide a full path to a batch.txt file as the /Directory/ (see README for format guidelines).'
     ,default=False,action='store_const',const=True,metavar='')
    parser.add_argument('-v','--visualize',help='Specify for visualization mode (default is processing mode). If selected, must input a comma-separated list of factor values and an adjusted p-value cutoff (ex. --visualize flight,ground,0.1) to compare. Outputs an interactive scatter plot of given conditions.',
        ,default=False,metavar='')


    #If user does not provide any arguments, simply display help message
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)


    #Parse user-provided arguments 
    args = parser.parse_args()
    batch = args.batch
    indir = args.Directory
    outdir = args.output
    condition1,condition2,pval_cut = args.visualize.split(',')


    #Get full paths to locations within this package
    srcdir = os.path.dirname(os.path.realpath(__file__))
    tempdir = os.path.join(os.path.dirname(srcdir),'temp')
    R_dir = os.path.join(os.path.dirname(srcdir),'R_scripts')


    #Write full paths to locations to a config.py file to be used by other scripts in this package
    with open(os.path.join(srcdir,'config.py'),'w') as outfile:
        outfile.write('indir = "' + indir + '"\n')
        outfile.write('outdir = "' + outdir + '"\n')
        outfile.write('srcdir = "' + srcdir + '"\n')
        outfile.write('tempdir = "' + tempdir + '"\n')
        outfile.write('R_dir = "' + R_dir + '"\n')
        outfile.write('batch = "' + str(batch) + '"\n')
        outfile.write('condition1 = "' + condition1 + '"\n')
        outfile.write('condition2 = "' + condition2 + '"\n')
        outfile.write('pval_cut = "' + pval_cut + '"\n')


    #Either run batch module or just run the processing steps on a single dataset
    if visualize == False:
        if batch:
            print "Batch option specified.\nUsing batch file: " + indir + "\nWriting output to: " + outdir
            import batch_process
            batch_process.run(indir)
        else:
            import metadata_process, rawdata_process
            print "Processing " + indir + "\nWriting output to: " + outdir
            GLDS = os.path.basename(indir)
            rawdata_out = os.path.join(outdir,GLDS,'microarray')
            metadata_out = os.path.join(outdir,GLDS,'metadata')
            metadata_in = os.path.join(indir,'metadata')
            rawdata_in = os.path.join(indir,'microarray')
            if os.path.isdir(metadata_dir):
                metadata_process.clean(metadata_in)
            else:
                raise IOError('metadata directory within input not found. See README for expected directory structure.')

            #Copy rawdata into output
            if os.path.isdir(rawdata_dir):
                rawdata_process.copy(rawdata_in)
                rawdata_process.rename(os.path.join(outdir,GLDS))
                rawdata_process.qc_and_normalize(rawdata_out,GLDS)
            else:
                raise IOError('microarray directory within input not found. See README for expected directory structure.')

            print "done."
    else:
        condition1,condition2,pval_cut = visualize.split(',')
        print "Visualization mode specified.\nComparing: " + condition1 + " vs. " + condition2 + "\nAdjusted p-value cutoff set at: " + pval_cut
        import differential_plot
        rawdata_out = os.path.join(outdir,'microarray')
        metadata_out = os.path.join(outdir,'metadata')
        GLDS = os.path.basename(outdir)
        rawdata_process.limma_differential(rawdata_out,metadata_out,GLDS)
        differential_plot.differential_visualize(rawdata_out,GLDS)
        print "done. Output in: " + rawdata_out


