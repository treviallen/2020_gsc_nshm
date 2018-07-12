from calc_oq_gmpes import hdf5_gsim, inslab_gsims
from numpy import logspace, sqrt, array, exp
import matplotlib.pyplot as plt
from os import path, getcwd
from scipy.constants import g
import matplotlib as mpl
mpl.style.use('classic')


mag  = 6.0
depths = [50., 60.]
ztor = 50. # guess
rake = 90. # USGS CMT
dip  = 30.

# set site details
vs30 = 760.
rjb = logspace(1,2.7,10)
rjb = array([55, 65])
rrup = rjb #sqrt(rjb**2 + dep**2) # assume point source; i.e. repi = rjb
rhypo = rjb

fig = plt.figure(1, figsize=(18, 9))
wdir = getcwd()

for i, dep in enumerate(depths):

    # get ground motion estimates from GMPEs
    #for i in range(0, len(rrup)):
    #Yea97imt, AB03imt, AB03CISimt, Gea05imt, Zea06imt, Zea06CISimt, MP10imt, AA13imt, Aea15imt, Zea16imt \
    Yea97imt, AB03imt, AB03CISimt, Gea05imt, Zea06imt, Zea06CISimt, MP10imt, Aea15imt \
        = inslab_gsims(mag, dep, ztor, dip, rake, rrup[1], rjb[1], vs30)
    
    # plot Garcia
    plt.subplot(1,2,i+1)    
    hdf5file = path.join(wdir, 'gmm_hdf5_tables', 'GarciaEtAl2005SSlab.vs1100.h50.hdf5') # OQ will not adjust using VS30, so default use VS30 table
    Geahdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], rhypo[1], vs30, hdf5file)
    # M6
    
    plt.loglog(array(Geahdf5imt['per']), exp(Geahdf5imt['sa']), 'r-', lw=3, label='G05 GSIM Table')
    plt.loglog(array(Gea05imt['per']), exp(Gea05imt['sa']), 'b-', lw=1.5, label='G05 Parametric')
    
    # plt AB03imt
    hdf5file = path.join(wdir, 'gmm_hdf5_tables', 'AtkinsonBoore2003SSlabCascadia.vs760.h50.hdf5')
    AB03hdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], rhypo[1], vs30, hdf5file)
    # M6
    
    plt.loglog(array(AB03hdf5imt['per']), exp(AB03hdf5imt['sa']), 'g-', lw=3, label='AB03 GSIM Table')
    plt.loglog(array(AB03CISimt['per']), exp(AB03CISimt['sa']), 'k-', lw=1.5, label='AB03 Parametric')
    
    # plot from table
    '''
    if i == 1:
       tabPer = Geahdf5imt['per'][::-1]
       tabSA = array([2.477, 2.541, 2.633, 2.667, 2.503, 2.348, 2.159, 2.019, 1.748, 1.560, 1.262, 1.030, 0.697, 0.445, 0.250, -0.294])
       plt.loglog(tabPer, (10**tabSA)/(g*100.), 'g--', lw=2, label='ASCII Table')
    '''
    plt.legend()
    
    plt.xlabel('Period (sec)')
    plt.ylabel('Spectral Acceleration (g)')
    plt.title(' '.join(('Garcia et al (2005); MW =', str(mag) + '; Rrup =',str('%0.1f' % rrup[1]), 'km; h =',str('%0.0f' % dep), 'km')))
    plt.grid(which='both')
    plt.xlim([0.01, 10.])
    
    """
    plt.subplot(122)    
    hdf5file = path.join(wdir, 'hdf5', 'A.15_SP15_adjusted_spec.hdf5')
    SPhdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], vs30, hdf5file)
    tabtxt = '-0.278 0.011 0.360 0.541 0.773 1.086 1.298 1.584 1.776 2.027 2.154 2.304 2.391 2.488 2.595 2.707 2.754 2.769 2.752 2.710 2.671 2.624 2.481'
    # M5
    #tabtxt = '-1.374 -0.966 -0.596 -0.418 -0.177 0.190 0.463 0.850 1.113 1.455 1.625 1.824 1.938 2.066 2.211 2.371 2.450 2.502 2.500 2.469 2.431 2.385 2.230'
    tabvals = (10**array([float(x) for x in tabtxt.split()]) / (100 * g))[::-1]
    
    plt.loglog(array(SP16imt['per']), exp(SP16imt['sa']), 'b-', lw=2, label='Parametric')
    plt.loglog(array(SPhdf5imt['per']), exp(SPhdf5imt['sa']), 'r-', lw=2, label='GSIM Table')
    plt.loglog(array(SPhdf5imt['per']), tabvals, 'g--', lw=2, label='Tables')
    plt.legend()
    
    plt.xlabel('Period (sec)')
    plt.ylabel('Spectral Acceleration (g)')
    plt.title(' '.join(('Shahjouei & Pezeshk (2016); MW =', str(mag) + '; Rrup =',str('%0.1f' % rrup[1]), 'km')))
    plt.grid(which='both')
    plt.xlim([0.01, 10.])
    """
plt.savefig('hdf5_test_M'+str(mag)+'_vs'+str('%0.0f' % vs30)+'.png', fmt='png', bbox_inches='tight')

plt.show()

