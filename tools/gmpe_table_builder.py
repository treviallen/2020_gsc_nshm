#!/usr/bin/env/python

"""
Canada GMPE tables to hdf5 - modified from code written by Graeme Weatherill
"""
import argparse
import h5py
import numpy as np
#import pandas as pd
from linecache import getlines
from collections import OrderedDict
from scipy.constants import g
from openquake.hazardlib.gsim.base import CoeffsTable
from openquake.hazardlib.imt import PGA, PGV, SA
from openquake.hazardlib import const
RLIST = ["Rjb", "Repi", "Rrup", "Rhypo"]


sigma = np.array([[0.05, 0.53],
                  [1.0 / 3.3, 0.53],
                  [0.5, 0.576],
                  [1.0, 0.622],
                  [10.0, 0.622]])
pga_sigma = 0.53
pgv_sigma = 0.622

SIGMA_TABLE = CoeffsTable(sa_damping=5, table="""
    IMT     Sigma
    pgv     0.622
    pga     0.530
    0.050   0.530
    0.100   0.530
    0.200   0.530
    0.3003  0.553
    0.500   0.576
    1.000   0.622
    2.000   0.622
    5.000   0.622
    10.00   0.622
    """)

def convert_data(data, data_type, units, is_log):
    """
    Converts data from log units to ground motion values, and from specified
    units to "g" if not input in "g"
    """
    if is_log:
        data = 10.0 ** data
    if data_type == "PGV":
        return data
        
    if (units == "cm/s/s") or (units == "cgs"):
        data /= (100.0 * g)
    elif units == "m/s":
        data /= g
    else:
        pass
    return data


def read_in_tables_simple(filename, rtypes, acc_units="cm/s/s", is_log=False, 
        skiprows=6,sigma=False):
    """
    Reads in GMPE tables in the EQHaz format
    """
    raw = getlines(filename)[skiprows:]
    # First row is header row
    headers = (raw[0].rstrip("\n")).split()
    data = []
    for row in raw[1:]:
        data.append(map(float, (row.rstrip("\n")).split('\t')))
    data = np.array(data)
    predictors = {"Mw": data[:, 0]}
    spectra = []
    spectra_idx = []
    pga = []
    pgv = []
    #print headers
    for iloc, header in enumerate(headers):
        if "Mw" in header:
            continue
        elif header in RLIST:
            if header in rtypes:
                if header == "Rrup":
                    predictors["rrup"] = data[:, iloc]
                else:
                    predictors[header.lower()] = data[:, iloc]
            else:
                continue
        elif int(float(header)) == 99:
            # PGA
            #pga = data[:, iloc] / (100.0 * g)
            pga = convert_data(data[:, iloc], "PGA", acc_units, is_log)
        elif int(float(header)) == 89:
            # PGV
            #pgv = data[:, iloc]
            pgv = convert_data(data[:, iloc], "PGV", acc_units, is_log)
        else:
            spectra.append(1.0 / float(header))
            spectra_idx.append(iloc)
    spectra.reverse()
    spectra_idx.reverse()
    imtls = {
        "Periods": spectra,
        "IMLs": convert_data(data[:, spectra_idx], "SA", acc_units, is_log)}
        #     "IMLs": data[:, spectra_idx] / (100.0 * g)}
    if len(pga):
        imtls["PGA"] = pga
    if len(pgv):
        imtls["PGV"] = pgv

    return predictors, imtls


def read_in_tables_gsc(filename, rtype):
    """
    Reads in GMPE tables in the format provided by GSC
    """
    print 'Reading tables in GSC format...'
    
    raw = getlines(filename)[2:]
    
    headers = (raw[0].rstrip("\n").replace('#', ' ')).split()

    sigmas = map(float, (raw[1].rstrip("\n").replace('#', ' ').split()))
    #print 'sigmas', len(sigmas), sigmas
    
    data = []
    for row in raw[2:]:
        data.append(map(float, (row.rstrip("\n")).split()))
    data = np.array(data)
    predictors = {"Mw": data[:, 0], rtype: data[:, 1]}
    
    # Reverse columns
    headers.reverse()
    data = data[:, ::-1]
    	
    # reorder sigma and re-append PGA/PGV at end
    oldSigmas = sigmas
    sigmas = sigmas[::-1][2:]
    print sigmas
    sigmas.append(oldSigmas[-2])
    sigmas.append(oldSigmas[-1])
    
    '''
    np.savetxt('hdf5_example.txt',(10.0 ** data[:, 2:-2]) / (100.0 * g), \
               delimiter='\t', fmt='%0.5e')
    '''
    imtls = {"Periods": np.array(map(float, headers[2:])),
             "IMLs": (10.0 ** data[:, 2:-2]) / (100.0 * g),
             "PGA": (10.0 ** data[:, 1]) / (100.0 * g),
             "PGV": 10.0 ** data[:, 0],
             "Sigma": sigmas}
    
    return predictors, imtls, sigmas



def build_hdf5_tables(filename, predictors, imtls):
    """
    Constructs the output hdf5 table
    """
    #print len(imtls["Periods"])
    #print imtls["IMLs"][0]
    #print imtls.keys()
    fle = h5py.File(filename, "w-")
    unique_mags = np.unique(predictors["Mw"])
    nmags = len(unique_mags)
    ndists = len(predictors["Mw"]) / nmags
    npers = imtls["IMLs"].shape[1]
    mag_dset = fle.create_dataset("Mw", (nmags,), dtype="f")
    mag_dset[:] = unique_mags
    #dist_grp = fle.create_group("Distances")
    for key in predictors.keys():
        if key == "Mw":
            continue
        else:
            dist_dset = fle.create_dataset("Distances",
                                           (ndists, 1, nmags),
                                           dtype="f")
            dist_dset.attrs["metric"] = key.lower()
            for iloc, mag in enumerate(unique_mags):
                idx = np.fabs(predictors["Mw"] - mag) < 1.0E-7
                dist_dset[:, :, iloc] = np.resize(predictors[key][idx],
                                                  [ndists, 1])
    iml_grp = fle.create_group("IMLs")
    spectra_dset = iml_grp.create_dataset("SA",
                                          (ndists, npers, nmags),
                                          dtype="f")
    if len(imtls["PGA"]):
        pga_dset = iml_grp.create_dataset("PGA",
                                          (ndists, 1, nmags),
                                          dtype="f")
    if len(imtls["PGV"]):
        pgv_dset = iml_grp.create_dataset("PGV",
                                          (ndists, 1, nmags),
                                          dtype="f")
    for iloc, mag in enumerate(unique_mags):
        idx = np.fabs(predictors["Mw"] - mag) < 1.0E-7
        spectra_dset[:, :, iloc] = imtls["IMLs"][idx, :]
        if len(imtls["PGA"]):
            pga_dset[:, :, iloc] = np.resize(imtls["PGA"][idx],
                                             [ndists, 1])
        if len(imtls["PGV"]):
            pgv_dset[:, :, iloc] = np.resize(imtls["PGV"][idx],
                                             [ndists, 1])
    per_dset = iml_grp.create_dataset("T", (npers,), dtype="f")
    per_dset[:] = np.array(imtls["Periods"])
    
    # Sigma Group
    sigma_grp = fle.create_group(const.StdDev.TOTAL)
    per_dset = sigma_grp.create_dataset("T", (npers,), dtype="f")
    per_dset[:] = np.array(imtls["Periods"])
    
    if len(imtls["PGA"]):
        pga_sigma_dset = sigma_grp.create_dataset("PGA",
                                                  (ndists, 1, nmags),
                                                  dtype="f")
    if len(imtls["PGV"]):
        pgv_sigma_dset = sigma_grp.create_dataset("PGV", 
                                                  (ndists, 1, nmags),
                                                  dtype="f")
    sa_sigma_dset = sigma_grp.create_dataset("SA", (ndists, npers, nmags))
    
    #sa_sigmas = np.array([SIGMA_TABLE[SA(period)]["Sigma"]
    #                     for period in imtls["Periods"]])
                         
    sa_sigmas = imtls["Sigma"][0:-2]
    
    for iloc in xrange(nmags):
        if len(imtls["PGA"]):
            #pga_sigma_dset[:, :, iloc] = SIGMA_TABLE[PGA()]["Sigma"] *\
            #    np.ones([ndists, 1])
                
            pga_sigma_dset[:, :, iloc] = imtls["Sigma"][-2] *\
                np.ones([ndists, 1])
            
        if len(imtls["PGV"]):
            #pgv_sigma_dset[:, :, iloc] = SIGMA_TABLE[PGV()]["Sigma"] *\
            #    np.ones([ndists, 1])
                
            pgv_sigma_dset[:, :, iloc] = imtls["Sigma"][-1] *\
                np.ones([ndists, 1])
                
        for jloc in xrange(ndists):
            sa_sigma_dset[jloc, :, iloc] = sa_sigmas
    fle.close()

def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert NRML stochastic event set file to tab delimited '
            ' .txt files. Inside the specified output directory, create a .txt '
            'file for each stochastic event set.'
            'To run just type: python eventset_converter.py '
            '--input-file=PATH_TO_INPUT_TABLE '
            '--output-dir=PATH_TO_OUTPUT_BINARY', add_help=False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--input-file',
        help='path to ses NRML file (Required)',
        default=None,
        required=True)
    flags.add_argument('--output-file',
        help='path to output directory (Required, raise an error if it already exists)',
        default=None,
        required=True)
    flags.add_argument('--distance-key',
        help='Preferred distance type',
        default=None,
        required=True)
    flags.add_argument('--is-log',
        help='Ground motion values givne in logarithmic scale (True) or not (False)',
        default=False,
        required=False)
    flags.add_argument('--skip-rows',
        help='Number of header rows to skip',
        default=6,
        required=False)
    flags.add_argument('--accel-units',
        help='Units of acceleration (cgs, cm/s/s, m/s, g)',
        default='cgs',
        required=False)
    flags.add_argument('--sigma-row',
        help="File has sigma row (True/False)",
        default=False,
        required=False)
    return parser


if __name__=="__main__":
    
    parser = set_up_arg_parser()
    args = parser.parse_args()
    print args.input_file
    predictors, imtls, sigmas = read_in_tables_gsc(args.input_file, 
                                           args.distance_key)
    #predictors, imtls = read_in_tables_simple(args.input_file, 
    #                                          [args.distance_key],
    #                                          args.accel_units,
    ##                                          args.is_log,
    #                                          int(args.skip_rows),
    #                                            args.sigma_row)
    print args.output_file
    build_hdf5_tables(args.output_file, predictors, imtls)
