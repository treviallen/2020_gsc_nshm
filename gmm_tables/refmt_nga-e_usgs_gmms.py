from __future__ import print_function
import xlrd
from os import path, system
from numpy import array, hstack, ones_like, interp, unique, loadtxt, savetxt, log10, log, zeros_like
from scipy.constants import g
from misc_tools import listdir_extension

vsAdjust = True # adjust from 3000 - 450
weightTab = False # output weighted tables
altSigma = False # do TA alternative sigma - see ngae_working/ngae_test_sigma.png

'''
params for VS30 adjustments
'''
vsPeriod = array([0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 5, 10])

bcVsFact = array([1.02, 1.36, 1.65, 1.67, 1.57, 1.37, 1.27, 1.17, 1.08]) # boore-campbell factors to convert 3000-760 m/s

aaVs30Fact = array([1.164, 1.14, 1.176, 1.259, 1.369, 1.443, 1.466, 1.481, 1.406]) # 2015 NBCC factors to convert 750-450 m/s

totalVsFact = bcVsFact * aaVs30Fact

'''
# USGS selected weights P 120 from PEER 2017/03
'''
# load period dependent weights for USGS NGA-East GMMs
wts = loadtxt(path.join('ngae_working', 'usgs_ngae_wgts.txt'), skiprows=1)
wthead = open(path.join('ngae_working', 'usgs_ngae_wgts.txt')).readlines()[0].split()
wthead = [x.strip('f=') for x in wthead]


# load xls gmm files for 13 USGS NGA-E tables
folder =  path.join('ngae_working', 'peer_usgs_xls')
extension = 'xls'

xlstables = listdir_extension(folder, extension)
for xls in xlstables:
    
    # get model num
    modnum = int(xls[:-4].split('_')[-1])
    
    # get model weights
    for wt in wts:
        if wt[0] == modnum:
            modwts = wt[1:]
    
    xlsfile = path.join(folder, xls)
    
    # Open the workbook
    xl_workbook = xlrd.open_workbook(xlsfile)
    
    # List sheet names, and pull a sheet by name
    sheet_names = xl_workbook.sheet_names()
    
    # get mag distance columns
    xl_sheet = xl_workbook.sheet_by_name(sheet_names[0])
    
    mag = []
    rrup = []
    for row_idx in range(1, xl_sheet.nrows):    # Iterate through rows - skip header
        mag.append(xl_sheet.cell(row_idx, 0).value)
        rrup.append(xl_sheet.cell(row_idx, 1).value)
        
    mag = array(mag).reshape((len(mag), 1))
    rrup = array(rrup).reshape((len(rrup), 1))
    
    # stack mag & dist
    allsa = hstack((mag, rrup))
    allsa_wt = hstack((mag, rrup))
    
    # now loop thru tabs and get period data
    perstr = []

    for sheet, fwt in zip(sheet_names, modwts):
        sa = []
        
        # get periods
        perstr.append(sheet.strip('F'))
        
        # period to get amp factors
        if sheet.strip('F') == 'PGA':
        	# assume T = 0.05 for amp factor
        	period = 0.05
        
        elif sheet.strip('F') == 'PGV':
            # assume T = 1.0 for amp factor
            period = 1.0
        
        else:
            period = float(sheet.strip('F'))
        
        # read each sheet
        xl_sheet = xl_workbook.sheet_by_name(sheet)
        
        # loop thru rows
        for row_idx in range(1, xl_sheet.nrows):    # Iterate through rows - skip header
            # get sa value
            saValue = float(xl_sheet.cell(row_idx, 2).value)
            
            if vsAdjust == True:
                # interpolate to get amp factor
                ampFact = interp(log(period), log(vsPeriod), totalVsFact)
                
                # now apply amp factor
                sa.append(saValue * ampFact)
            else:
                sa.append(saValue)
        
        # convert to cgs
        if savetxt == 'PGV':
            #sawt = log10(array(sa) * fwt) # arithmetic
            sawt = log10(array(sa)) * fwt # geometric
            sa = log10(array(sa))
            #sawt = ones_like(sa) * fwt # keep for testing
        else:
            #sawt = log10(array(sa) * 100 * g * fwt) # arithmetic
            sawt = log10(array(sa) * 100 * g) * fwt # geometric
            sa = log10(array(sa) * 100 * g)
            #sawt = ones_like(sa) * fwt # keep for testing
        
        # add sa to big array
        allsa = hstack((allsa, sa.reshape((len(sa), 1))))
        allsa_wt = hstack((allsa_wt, sawt.reshape((len(sawt), 1))))
    
    ####################################################################
    # get sigma data from AA13
    ####################################################################
    if altSigma == False:
        pernum = array([1./float(x) for x in perstr[0:-2]])
        
        shortT = 0.25
        longT  = 1.0
        
        # set periods in between based on AA13
        log10sig = interp(log10(pernum[::-1]), [log10(shortT), log10(longT)], [0.23, 0.27], left=0.23, right=0.27)[::-1]
        
        # append pga & pgv
        log10sig = hstack((log10sig, [0.23, 0.27]))
        lnsig = log(10**log10sig)
    
    ####################################################################
    #  crude NGA-E sigma
    ####################################################################
    elif altSigma == True:
    
        # get periods
        pernum = array([1./float(x) for x in perstr[0:-2]])
        
        shortT = 0.01
        midT = 0.1
        longT  = 10.0
        
        shortTsig = 0.70
        midTsig = 0.81
        longTsig = 0.57
        
        lnsig = zeros_like(pernum)
        #print lnsig
        
        idx = pernum <= midT
        lnsig[idx] = interp(log10(pernum[idx][::-1]), [log10(shortT), log10(midT)], [shortTsig, midTsig])[::-1]
        
        idx = pernum > midT
        lnsig[idx] = interp(log10(pernum[idx][::-1]), [log10(midT), log10(longT)], [midTsig, longTsig])[::-1]
        
        # add PGA & PGV sigmas
        lnsig = hstack((lnsig, [0.700, 0.810]))
    
    ####################################################################
    # write to table
    ####################################################################
    basename = path.split(xls)[1].strip('.xls')
    header = basename+' NGA-East USGS GMM as published in PEER Report No. 2017/03, distance is Rrup. Log10 hazard values in cgs units. Uses NRCan simplified, magnitude-independent sigma'
    config = '          ' + ' '.join((str(len(unique(mag))), str(len(unique(rrup))), str(len(sheet_names)), ': nmag, ndist, nperiod'))
    pertxt = '          ' + ' '.join([str('%0.3f' % x) for x in pernum]) + ' PGA PGV'
    sigtxt = '          ' + ' '.join([str('%0.2f' % x) for x in lnsig])
    
    # combine headers
    headers = '\n'.join((header, config, pertxt, sigtxt))
    
    
    # write headers
    if vsAdjust == True:
        outfile = path.join('ngae_usgs_txt_tables', basename+'.vs450.txt')
        
        if weightTab == True:
            outfile_wt = path.join('weighted_usgs_ngae_tables', basename+'_weighted.450mps.txt')
            
        if altSigma == False:
            outfile = path.join('ngae_usgs_txt_tables', basename+'_AA13_sigma.vs450.txt')
        
    else:
        outfile = path.join('ngae_usgs_txt_tables', basename+'.vs3000.txt')
        
        if weightTab == True:
            outfile_wt = path.join('weighted_usgs_ngae_tables', basename+'_weighted.txt')
            
        if altSigma == False:
            outfile = path.join('ngae_usgs_txt_tables', basename+'_AA13_sigma.vs3000.txt')
        #outfile_wt = basename+'_weighted.txt'
    
    # now append numpy array
    if weightTab == False:
        savetxt(outfile, allsa, fmt='%.5f', delimiter=' ', header=headers, comments='')
    else:
        savetxt(outfile_wt, allsa_wt, fmt='%.5f', delimiter=' ', header=headers, comments='')
    
