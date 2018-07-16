from calc_oq_gmpes import hdf5_gsim, inslab_gsims, interface_gsims
from numpy import logspace, sqrt, array, exp, interp, log, arange
import matplotlib.pyplot as plt
from os import path, getcwd
from scipy.constants import g
import matplotlib as mpl
mpl.style.use('classic')


mag  = 8.
dep = 15.
ztor = 10. 
rake = 90. 
dip  = 30.

# set site details
rjb = array([30]) # only use dist 1
rrup = sqrt(rjb**2 + dep**2) # assume point source; i.e. repi = rjb
rhypo = rrup

fig = plt.figure(1, figsize=(12, 10))
wdir = getcwd()

mags = arange(7.5, 9.1, 0.1)
periods = [0.0, 0.2, 1.0, 2.0]

for i, period in enumerate(periods):
    magAmpHDF = []
    magAmpOQ = []
    for mag in mags:

        ##########################################################################
        # do Site Class B/C
        ##########################################################################

        ax = plt.subplot(2,2,i+1)
        vs30 = 760.
        
        
        # plt AtkinsonMacias2009
        hdf5file = path.join(wdir, 'gmm_hdf5_tables', 'AtkinsonMacias2009.vs760.h15.hdf5')
        table_hdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[0], rjb[0], rhypo[0], vs30, hdf5file)
        
        
        # plt OQ implementations
        Yea97imt, AB03imt, Zea06imt, Zea06CISimt, AM09imt, MP10imt, GA14imt, GA14CISimt, Aea15imt \
            = interface_gsims(mag, dep, ztor, dip, rake, rrup[0], rjb[0], vs30)
        
        # interp to period of interest
        if period == 0.0:
            magAmpHDF.append(table_hdf5imt['pga'][0])
            magAmpOQ.append(AM09imt['pga'][0])
        else:
            magAmpHDF.append(interp(log(period), log(table_hdf5imt['per']), table_hdf5imt['sa'])) # sa already log
            magAmpOQ.append(interp(log(period), log(AM09imt['per']), AM09imt['sa'])) # sa already log
        
    plt.semilogy(mags, exp(magAmpHDF), 'b-', lw=2, label='Atkinson & Macias (2009) - HDF')
    plt.semilogy(mags, exp(magAmpOQ), '--', c='r', lw=3., label='Atkinson & Macias (2009) - OQ')
        
    plt.xlabel('Period (sec)')
    plt.ylabel('Spectral Acceleration (g)')
    plt.suptitle('Rrup = '+str('%0.1f' % rrup[0]))
    plt.title('SA'+str(period))
    plt.grid(which='both')
    #plt.xlim([0.01, 10.])
    
    if i == 0:
        plt.legend(loc=3, fontsize=11)


plt.savefig('am09_hdf5_test.rrup_'+str('%0.1f' % rrup[0])+'.png', fmt='png', bbox_inches='tight')


plt.show()

