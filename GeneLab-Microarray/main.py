__author__ = 'Jonathan Rubin'

import argparse
import sys

def run():
    parser = argparse.ArgumentParser(prog='GeneLab-Microarray',usage='%(prog)s [-c ChipType]',description='Standardized processing pipeline for microarray data on GeneLab.')
    parser.add_argument('-c','--chip',help='Specify the type of microarray chip. Options: Affymetrix, Agilent, Nimblegen, Illumina, Custom',
        choices=['Affymetrix', 'Agilent', 'Nimblegen', 'Illumina', 'Custom'],metavar='',required=True)
    if len(sys.argv)==1:
        # display help message when no args are passed. 
        parser.print_help()
        sys.exit(1)

    chip = parser.parse_args().chip
    print chip