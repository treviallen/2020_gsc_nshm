# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 11:40:40 2015

@author: tallen
"""

from misc_tools import listdir_extension
from os import system, path, remove
from sys import argv

'''
# first remove hdf5 files
hdf5files = listdir_extension('hdf5', 'hdf5')
if len(hdf5files) > 0:
    for hf in hdf5files:
        remove(hf)
'''
# get xls gmm files
#tabin = argv[1] #'tables' #path.join('2015-08_Electronic Appendices','AppendixA_AdjustedGMMs')
#dist = argv[2]


extension = 'txt'
gmpetables = listdir_extension(folder, extension)

'''
usage: gmpe_table_builder.py [-h] --input-file INPUT_FILE --output-file
                             OUTPUT_FILE --distance-key DISTANCE_KEY
                             [--is-log IS_LOG] [--skip-rows SKIP_ROWS]
                             [--accel-units ACCEL_UNITS]
                             [--sigma-row SIGMA_ROW]
'''
#gmpetables = ['A.10_2CCSP_adjusted_spec.txt']

for tab in gmpetables:

    '''
    # get params for converting to hdf5
    if tab.startswith('ENA') or tab.startswith('Wcrust_') \
       or tab.startswith('Winslab') or tab.startswith('Woffshore'):
        dist = 'Rhypo'
    elif tab.startswith('WcrustFRjb'):
        dist = 'Rjb'
    elif tab.startswith('WinterfaceCombo'):
        dist = 'Rrup'
    else:
        dist = 'Rrup' # for NGA-E
    '''
        
    tabin  = path.join(folder, tab)
    
    hdf5out = path.join('gmm_hdf5_tables', path.split(tab)[-1].strip('txt')+'hdf5')
    
    try:
        remove(hdf5out)
    except:
        print 'Continue: hdf5 file does not exist'
    
    # read in table and get diatance metric
    lines = open(tab).readlines()
    distType = lines[0].split('.')[0].split()[-1].capitalize()
    print distType
    
    system(''.join(('python /Users/tallen/Documents/Code/my_codes/gmpe_table_builder.py --input-file=', \
                    tabin, ' --output-file=', hdf5out, ' --distance-key=', distType, ' --is-log True')))
    
    

