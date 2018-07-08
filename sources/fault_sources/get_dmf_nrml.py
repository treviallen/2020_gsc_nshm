from numpy import array, arange, log10, exp, array, around, sin, pi, nan, radians, where
from fault_tools import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import operator
from os import path, sep
import shapefile
from mapping_tools import get_field_data
#from matplotlib import mpl

def make_collapse_occurrence_text(mmin, mc, binwid, bval_char, area, fwt):
    from numpy import zeros   
    from fault_tools import characteristic_mag_rec_from_slip 
    
    sr_wt    = [0.16, 0.68, 0.16]
    max_mag_wt = [0.60, 0.30, 0.10]
    srate = [0.15, 0.25, 0.35]
    wtd_list  = []
    maglen = 0

    for i, swt in enumerate(sr_wt):
        
        mchar = [mc-0.15, mc, mc+0.15]
        for j, mwt in enumerate(max_mag_wt):
            
            # get characteristic log10 A value
            cA0, mx, n_min_mag, n_char_mag = characteristic_mag_rec_from_slip( \
                                            srate[i], mu, bval_char, mmin, mchar[j], \
                                            area=area, ftype='ss')
                                                
            cN0 = 10**cA0
            
            cmrate, cslip_per_mag, ccumrate = get_cummulative_stats(cA0, bval_char, mchar[j], \
                                              'characteristic', mrange, n_char_mag, fault)
                                      
            wtd_list.append(array(cmrate)	 * swt * mwt *fwt)  
            
            # get max length of arrays
            if len(cmrate) > maglen:
                maglen = len(cmrate)

    # sum rates for each branch
    wtd_rates = zeros(maglen)
    for rates in wtd_list:
        # go element by element
        for r, rate in enumerate(rates):
            wtd_rates[r] += rate
    
    """            
    # convert cummulative rates to annual occurrence rates
    occ_rates = []
    for b in range(0, len(wtd_rates[0:-1])):
        occ_rates.append(wtd_rates[b] - wtd_rates[b+1])
    occ_rates.append(wtd_rates[-1])
    """
    
    # make text object                        
    octxt = str('%0.5e' % wtd_rates[0])
    for bc in wtd_rates[1:]:
        octxt += ' ' + str('%0.5e' % bc)
        
    return octxt



# read shapefile with fault length and sliprate
root = path.split('shapefiles\\DMF_LRF.shp')[0]
shpfile = path.join('shapefiles', 'DMF_LRF.shp')
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
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
tectonic_reg = 'Active Shallow Fault'
for i in range(0,len(flen)):
    
    '''
    Note: for hazard model, reduce mchar by 0.25 and assume mchar = mmax for
    imput to slip inversion
    '''
    ftype = 'ss'
    area.append(flen[i] * seis_thickness / sin(radians(fdip[i])))
    mchar.append(return_area2mag_for_ss(area[i])[1] - 0.25)
    sig.append(0.15) # assumed from various distributions


# Assumptions
bval_char = .8
beta = bval2beta(bval_char)

#mmin = [7.0, 7.0, 7.0, 7.0, 7.0, 7.0]
mu = 3.3E10 # Pa

# M range for plotting
bin_width = 0.1

# now loop thru the fault sources
mmin = 6.3
legend_txt = []
newxml = ''
for j, fault in enumerate(fname):
    
    N0rng = []
    
    # use best slip only
    sr = slipbest[j]

    mc = mchar[j]
    a = area[j]
    
    shape = shapes[j].points
            
    # set mrange
    mrange = arange(around(mmin, decimals=2), mc+0.25, bin_width)
    
    octxt = make_collapse_occurrence_text(mmin, mc, bin_width, bval_char, a, fwt[j])
    
    print octxt
    
    ##################################################################################
    # start making xml
    ##################################################################################
    
    if j > 0:
        newxml += '        <simpleFaultSource id="'+fcode[j]+'" name="'+ \
                             fname[j]+'" tectonicRegion="'+tectonic_reg+'">\n'
        newxml += '            <simpleFaultGeometry>\n'
        newxml += '                <gml:LineString>\n'
        newxml += '                    <gml:posList>\n'
        
        # simple faults use surface projection!
        '''
        # calculate lat lons from surface projection
        # get upper h-dist
        upperhdist = m['src_dep'][0] / tan(radians(m['fault_dip'][0]))
        upperxy = get_line_parallels(m['src_shape'], upperhdist)[0]
        '''
        
        xytxt = ''
        # reverse coords
        for xy in shape[::-1]: 
            xytxt += '                            ' + \
                     ' '.join((str('%0.4f' % xy[0]), str('%0.4f' % xy[1])))+'\n'
        newxml += xytxt
        
        newxml += '                    </gml:posList>\n'
        newxml += '                </gml:LineString>\n'
        newxml += '                <dip>'+'70'+'</dip>\n'
        newxml += '                <upperSeismoDepth>0.0</upperSeismoDepth>\n'
        newxml += '                <lowerSeismoDepth>15.0</lowerSeismoDepth>\n'
        newxml += '            </simpleFaultGeometry>\n'
        
        '''
        # get fault area scaling model
        '''
        src_code = fcode[j]
        if src_code == 'CIS':
            newxml += '            <magScaleRel>GSCCascadia</magScaleRel>\n'
        elif src_code.startswith('WIN'):
            newxml += '            <magScaleRel>GSCOffshoreThrustsWIN</magScaleRel>\n'
        elif src_code.startswith('HGT'):
            newxml += '            <magScaleRel>GSCOffshoreThrustsHGT</magScaleRel>\n'
        elif src_code.startswith('QCSS') or src_code.startswith('FWF'):
            newxml += '            <magScaleRel>WC1994_QCSS</magScaleRel>\n'
        elif src_code.startswith('EISO'):
            newxml += '            <magScaleRel>GSCEISO</magScaleRel>\n'
        elif src_code.startswith('EISB'):
            newxml += '            <magScaleRel>GSCEISB</magScaleRel>\n'
        elif src_code.startswith('EISI'):
            newxml += '            <magScaleRel>GSCEISI</magScaleRel>\n'
        else:
            newxml += '            <magScaleRel>WC1994</magScaleRel>\n'
        
        newxml += '            <ruptAspectRatio>1.0</ruptAspectRatio>\n'
        #newxml += '            <ruptAspectRatio>2.0</ruptAspectRatio>\n'
        '''
        # now get appropriate MFD
        '''
        # do incremental MFD
        if bval_char > -99:
            
                        
            # make text
            newxml += '            <incrementalMFD minMag="'+str('%0.2f' % (mmin+0.5*bin_width))+'" binWidth="'+str(bin_width)+'">\n'
            newxml += '                <occurRates>'+octxt+'</occurRates>\n'
            newxml += '            </incrementalMFD>\n'
                        
        if fdip[j] != 90.:
            newxml += '            <rake>90.0</rake>\n'
        else:
            newxml += '            <rake>0.0</rake>\n'
        
        newxml += '        </simpleFaultSource>\n\n'
    
# write xml to file
f = open('dmf_nrml.xml', 'wb')
f.write(newxml)
f.close()