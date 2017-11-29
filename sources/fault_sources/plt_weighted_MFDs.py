# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:57:46 2013

@author: tallen
"""


    
'''
START MAIN CODE NOW
'''

from numpy import array, arange, log10, exp, array, around, sin, pi, nan, radians, where
from fault_tools import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import operator
from os import path, sep
import shapefile
from mapping_tools import get_field_data
#from matplotlib import mpl


# read shapefile with fault length and sliprate
root = path.split('shapefiles//DMF_LRF.shp')[0]
shpfile = 'shapefiles//DMF_LRF.shp'
sf = shapefile.Reader(shpfile)
fname = get_field_data(sf, 'FAULT_NAME', 'str')
fcode = get_field_data(sf, 'CODE', 'str')
flen  = get_field_data(sf, 'LENGTH', 'float')
fdip  = get_field_data(sf, 'DIP', 'float')
fwt  = get_field_data(sf, 'WEIGHT', 'float')
slipbest = get_field_data(sf, 'SLIP_RATE', 'float')
slipmin = get_field_data(sf, 'SLIP_MIN', 'float')
slipmax = get_field_data(sf, 'SLIP_MAX', 'float')

# get characteristic magnitude from Wells & Coppersmith 1994
mchar = []
sig = []
area = []
seis_thickness = 15.
for i in range(0,len(flen)):
    
    '''
    Note: for hazard model, reduce mchar by 0.25 and assume mchar = mmax for
    imput to slip inversion
    '''
    ftype = 'ss'
    area.append(flen[i] * seis_thickness / sin(radians(fdip[i])))
    mchar.append(return_area2mag_for_ss(area[i])[1] - 0.25)
    sig.append(0.15) # assumed from various distributions
    

# make colourmap
#ncolours = len(sliprates) * len(mmax)
ncolours = 8
cmap = cm.get_cmap('jet', ncolours)
cs = cmap(arange(ncolours))

# Assumptions
bval_char = .8
beta = bval2beta(bval_char)

#mmin = [7.0, 7.0, 7.0, 7.0, 7.0, 7.0]
mu = 3.3E10 # Pa

# M range for plotting
bin_width = 0.01

# now loop thru the fault sources
mmin = 6.3
legend_txt = []
figure = plt.figure(1,figsize=(16,8))
for j, fault in enumerate(fname):
    
    plthandle = []
    
    
    N0rng = []
    
    # use best slip only
    sr = slipbest[j]

    mc = mchar[j]
            
    # set mrange
    mrange = arange(around(mmin, decimals=2), 7.61, bin_width)
    
#        for mc in mchar_range_fix:

    # get characteristic log10 A value
    cA0, mx, n_min_mag, n_char_mag = characteristic_mag_rec_from_slip( \
                                    sr, mu, bval_char, mmin, mc, \
                                    area=area[j], ftype='ss')
                                        
    cN0 = 10**cA0
    
    cmrate, cslip_per_mag, ccumrate = get_cummulative_stats(cA0, bval_char, mc, \
                                      'characteristic', mrange, n_char_mag, fault)
    
    # now plot non-weighted
    if j > 0: # don't plt full length
        col = cs[j]
        plt.subplot(121)
        plt.semilogy(mrange, ccumrate, linewidth=2, color=col)
        
        # plot weighted
        plt.subplot(122)
        plt.semilogy(mrange, array(ccumrate)*fwt[j], linewidth=2, color=col)
        
        # add leg txt
        legend_txt.append(' '.join((fcode[j], str('%0.0f' % flen[j]), 'km')))

# now plot
plt.subplot(121)
plt.title('Non-Weighted Rates (Best Slip & Mmax)', fontsize=18)
plt.ylabel('Cumulative Rate (/yr)', fontsize=16)
plt.xlabel('Moment Magnitude', fontsize=16)
plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
plt.xlim([6.0, 7.5])
plt.ylim([10**-6, 10**-3])

plt.subplot(122)
plt.title('Weighted Rates (Best Slip & Mmax)', fontsize=18)
plt.ylabel('Weighted Cumulative Rate (/yr)', fontsize=16)
plt.xlabel('Moment Magnitude', fontsize=16)
plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
plt.xlim([6.0, 7.5])
plt.ylim([10**-6, 10**-3])

plt.legend(legend_txt, loc='upper right', fontsize=12)

# make file name
#outfile = root + sep + 'DMF_LRF_weighted_rates.png'
outfile = 'DMF_LRF_weighted_rates.png'

plt.savefig(outfile,format='png', bbox_inches='tight', dpi=300)
plt.show()