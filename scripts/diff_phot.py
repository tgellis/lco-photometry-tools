# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 10:18:38 2016

@author: rstreet
"""

import sys
from os import path
import dataset

def run_diff_phot():
    """Function to perform difference photometry for selected stars within
    a given dataset"""
    #Fixed input mode for debug
    params = {}
    params['input_dir'] = "/mnt/c/Users/tellis/Documents/Data/lcodata"
    params['output_dir'] = "/mnt/c/Users/tellis/Documents/Data/output"
    params['star_file'] = "/mnt/c/Users/tellis/Documents/Projects/targets.txt"

#    params = get_params()
    check_sanity(params)

    phot_data = dataset.DataSet(params=params)
    phot_data.load_star_list(params)
    phot_data.print_star_list()
    phot_data.load_data(params)

    phot_data.diff_photometry()

    phot_data.quality_control()

    phot_data.plot_lightcurves(phot_data.star_list)
    phot_data.output_lightcurves(phot_data.star_list)

    phot_data.plot_ensemble_lightcurve()

def get_params():
    """Function to gather the required commandline arguments for different
    photometry"""

    params = {}
    if len(sys.argv) != 4:
        params['input_dir'] = raw_input('Please enter the input directory path: ')
        params['output_dir'] = raw_input('Please enter the output directory path: ')
        params['star_file'] = raw_input('Please enter the path to the starlist file: ')
    else:
        params['input_dir'] = sys.argv[1]
        params['output_dir'] = sys.argv[2]
        params['star_file'] = sys.argv[3]

    return params

def check_sanity(params):
    """Function to perform basic checks to see if the reduction can continue"""

    for dpath in ['input_dir', 'output_dir']:
        if not path.isdir(params[dpath]):
            print "ERROR: Cannot find directory: %s"%params[dpath]
            sys.exit(0)

    if not path.isfile(params['star_file']):
        print "ERROR: Cannot find star file: %s"%params['star_file']
        sys.exit(0)


if __name__ == '__main__':
    run_diff_phot()
