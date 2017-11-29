# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:57:46 2013

@author: tallen
"""


    
'''
START MAIN CODE NOW
'''

from numpy import arange, log10, exp, array, around, sin, pi, nan, radians, where
from fault_tools import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import operator
from os import path, sep
import shapefile
from mapping_tools import get_field_data
#from matplotlib import mpl


# read shapefile with fault length and sliprate
root = path.split('DMF\\DMF_LRF.shp')[0]
shpfile = 'DMF\\DMF_LRF.shp'
sf = shapefile.Reader(shpfile)
fname = get_field_data(sf, 'FAULT_NAME', 'str')
fcode = get_field_data(sf, 'CODE', 'str')
flen  = get_field_data(sf, 'LENGTH', 'float')
fdip  = get_field_data(sf, 'DIP', 'float')
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
ncolours = 4
cmap = cm.get_cmap('hsv', ncolours)
cs = cmap(arange(ncolours))

# Assumptions
bval_char = .8
beta = bval2beta(bval_char)

#mmin = [7.0, 7.0, 7.0, 7.0, 7.0, 7.0]
mu = 3.3E10 # Pa

# M range for plotting
bin_width = 0.01

# now loop thru the fault sources
j = 0
#for fault in [fname[j]]:
mmin = 6.3
for j, fault in enumerate(fname):
    
    # get upper and lower mchar    
    mchar_range = [mchar[j]-sig[j], mchar[j], mchar[j]+sig[j]]
    
    plthandle = []
    legend_txt = []
    
    figure = plt.figure(j,figsize=(10,10))
    N0rng = []
    
    # compile min-max-best slip rates
    slipr = [slipmin[j], slipbest[j], slipmax[j]]
    for i, sr in enumerate(slipr):
#        ii = 0
        # !!!!! to use just extreme values - temp fix only!!!!
        mc = mchar_range[i]
                
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
        
        '''# get bounded curve (sliprate, mu, bval, mchar, mmin, **kwargs)
        bA0, mx = bounded_mag_rec_from_slip(sr, mu, bval_bound, mc, mmin, area=area[j])
        bN0 = 10**bA0
        bN  = 10**bA0 * exp(-beta*mrange) * (1 - exp(-beta*(mx - mrange)))
        
        bn_char_mag = 10**bA0 * exp(-beta*mc) * (1 - exp(-beta*(mx - mc)))
        
        bmrate, bslip_per_mag, bcumrate = get_cummulative_stats(bA0, bval_bound, mc, \
                                          'bounded', mrange, nan, fault)
        # get N0 range for GSCFRISK
        bN0frisk = bN0 * (1-exp(-beta*mx))
        N0rng.append(bN0frisk)
        '''
                
        # now plot
        col = [cs[i][0],cs[i][1],cs[i][2]]
        if i == 0:
            h1 = plt.semilogy(mrange, ccumrate, linewidth=1., color=col)
        elif i == 1:
            h2 = plt.semilogy(mrange, ccumrate, linewidth=1., color=col)
        elif i == 2:
            h3 = plt.semilogy(mrange, ccumrate, linewidth=1., color=col)
            
        #plt.semilogy(mrange, ccumrate, '--', linewidth=1., color=col)
        
        chidx = where((mrange > mc-bin_width/2) & (mrange < mc+bin_width/2))[0]
        
        legend_txt.append('\nMmin = '+str("%0.2f" % mmin)+'; Mchar = '+str("%0.2f" % mc) \
                          +'; Mmax = '+str("%0.2f" % mx) \
                          +';\nbeta = '+str("%0.2f" % beta)+'; SR = '+str("%0.2f" % sr)+ \
                          ' mm/yr; WC94 disp = '+str("%0.2f" % sum(cslip_per_mag))+' mm/yr;' \
                          +'\nMchar RP = '+str("%0.0f" % (1./ccumrate[chidx]))+' yrs')
    
#            i += 1
#            ii += 1

    # now plot
    plt.title(fault+' ('+fcode[j]+')\nSR = ' + str(slipr[1]) + ' (' + str(slipr[0]) + '-' \
              + str(slipr[2]) +') mm/yr; LEN = '+str('%0.0f' % flen[j])+' km')
    plt.ylabel('Cumulative Rate (/yr)')
    plt.xlabel('Moment Magnitude')
    plt.legend((h1[0],h2[0],h3[0]), legend_txt, loc='lower left')
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize='small')
    plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
    plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
    plt.xlim([mmin-0.5, mx+0.25])
    plt.ylim([10**-6, 10**-2])
    
    # make file name
    outfile = root + sep + reduce(operator.add, fault.strip().split(' ')) \
              + '_bval.eq.' + str("%0.1f" % bval_char) + '.png'
    
    plt.savefig(outfile,format='png')
    plt.show()
    
    '''
    # output parameters for GSCFRISK
    header = fcode[j] + ' - ' + fault
    mxrange = array(mchar_range)  + 0.25
    magline = ' '.join((str("%0.2f" % mmin), str("%0.2f" % mxrange[1]), \
                        str("%0.2f" % mxrange[0]), str("%0.2f" % mxrange[2])))
                        
    betaline = '  '+'  '.join((str("%0.4f" % N0rng[1]), str("%0.4f" % beta), \
                               str("%0.4f" % N0rng[0]), str("%0.4f" % beta), \
                               str("%0.4f" % N0rng[2]), str("%0.4f" % beta)))
    
    outtxt = '\n'.join((header, magline, betaline))
    outfile = path.join(shpfile.split('.')[0], fcode[j], fcode[j] + '.beta0')

    f = open(outfile, 'wb')
    f.write(outtxt)    
    f.close()
    '''
    
    j += 1