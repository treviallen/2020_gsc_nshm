from calc_oq_gmpes import hdf5_gsim, inslab_gsims, interface_gsims
from numpy import logspace, sqrt, array, exp
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
rjb = array([50, 100])
rrup = sqrt(rjb**2 + dep**2) # assume point source; i.e. repi = rjb
rhypo = rrup

fig = plt.figure(1, figsize=(18, 9))
wdir = getcwd()

#for i, dep in enumerate(depths):


##########################################################################
# do Site Class C
##########################################################################

ax = plt.subplot(121)    
vs30 = 450.

# plt ZhaoEtAl2006SInterCascadia
hdf5file = path.join(wdir, 'gmm_hdf5_tables', 'ZhaoEtAl2006SInterCascadia.vs450.h15.hdf5')
table_hdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], rhypo[1], vs30, hdf5file)
plt.loglog(array(table_hdf5imt['per']), exp(table_hdf5imt['sa']), 'g-', lw=2, label='Zhao et al (2006)')

# plt AtkinsonMacias2009
hdf5file = path.join(wdir, 'gmm_hdf5_tables', 'AtkinsonMacias2009.vs450.h15.hdf5')
table_hdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], rhypo[1], vs30, hdf5file)
plt.loglog(array(table_hdf5imt['per']), exp(table_hdf5imt['sa']), 'b-', lw=2, label='Atkinson & Macias (2009)')

# plt GhofraniAtkinson2014Cascadia
hdf5file = path.join(wdir, 'gmm_hdf5_tables', 'GhofraniAtkinson2014Cascadia.vs450.h15.hdf5')
table_hdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], rhypo[1], vs30, hdf5file)
plt.loglog(array(table_hdf5imt['per']), exp(table_hdf5imt['sa']), '-', c='orange', lw=2, label='Ghofrani & Atkinson (2014)')

# plot AbrahamsonEtAl2015SInter
hdf5file = path.join(wdir, 'gmm_hdf5_tables', 'AbrahamsonEtAl2015SInter.vs450.h15.hdf5')
table_hdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], rhypo[1], vs30, hdf5file)
plt.loglog(array(table_hdf5imt['per']), exp(table_hdf5imt['sa']), 'r-', lw=2, label='Abrahamson et al (2015)')

# plot 2015 NBCC
#hdf5file = path.join(wdir,'..', '..', '..', '2015_National_Hazard', '2015_gsc_nshm', 'gm_tables', 'WinterfaceCombo_medclC.h15.hdf5')
hdf5file = '/Users/tallen/Documents/2015_National_Hazard/2015_gsc_nshm/gm_tables/WinterfaceCombo_medclC.hdf5'
table_hdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], rhypo[1], vs30, hdf5file)
plt.loglog(array(table_hdf5imt['per']), exp(table_hdf5imt['sa']), 'k-', lw=2, label='2015 NBCC')

# plt OQ implementations
Yea97imt, AB03imt, Zea06imt, Zea06CISimt, AM09imt, MP10imt, GA14imt, GA14CISimt, Aea15imt \
    = interface_gsims(mag, dep, ztor, dip, rake, rrup[1], rjb[1], 450.)

plt.loglog(array(Zea06CISimt['per']), exp(Zea06CISimt['sa']), '--', c='m', lw=3., label='Zhao et al (2006) - OQ')
plt.loglog(array(AM09imt['per']), exp(AM09imt['sa']), '--', c='y', lw=3., label='Atkinson & Macias (2009) - OQ (B/C only)')


plt.xlabel('Period (sec)')
plt.ylabel('Spectral Acceleration (g)')
plt.title('Test 2020 NBCC Interface GMMs on Site Class C')
plt.grid(which='both')
plt.xlim([0.01, 10.])

plt.legend(loc=3, fontsize=11)

##########################################################################
# do Site Class B/C
##########################################################################

ax = plt.subplot(122)    
vs30 = 760.

# plt ZhaoEtAl2006SInterCascadia
hdf5file = path.join(wdir, 'gmm_hdf5_tables', 'ZhaoEtAl2006SInterCascadia.vs760.h15.hdf5')
table_hdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], rhypo[1], vs30, hdf5file)
plt.loglog(array(table_hdf5imt['per']), exp(table_hdf5imt['sa']), 'g-', lw=2, label='Zhao et al (2006)')

# plt AtkinsonMacias2009
hdf5file = path.join(wdir, 'gmm_hdf5_tables', 'AtkinsonMacias2009.vs760.h15.hdf5')
table_hdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], rhypo[1], vs30, hdf5file)
plt.loglog(array(table_hdf5imt['per']), exp(table_hdf5imt['sa']), 'b-', lw=2, label='Atkinson & Macias (2009)')

# plt GhofraniAtkinson2014Cascadia
hdf5file = path.join(wdir, 'gmm_hdf5_tables', 'GhofraniAtkinson2014Cascadia.vs760.h15.hdf5')
table_hdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], rhypo[1], vs30, hdf5file)
plt.loglog(array(table_hdf5imt['per']), exp(table_hdf5imt['sa']), '-', c='orange', lw=2, label='Ghofrani & Atkinson (2014)')

# plot AbrahamsonEtAl2015SInter
hdf5file = path.join(wdir, 'gmm_hdf5_tables', 'AbrahamsonEtAl2015SInter.vs760.h15.hdf5')
table_hdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], rhypo[1], vs30, hdf5file)
plt.loglog(array(table_hdf5imt['per']), exp(table_hdf5imt['sa']), 'r-', lw=2, label='Abrahamson et al (2015)')

# plot 2015 NBCC
#hdf5file = path.join(wdir,'..', '..', '..', '2015_National_Hazard', '2015_gsc_nshm', 'gm_tables', 'WinterfaceCombo_medclC.h15.hdf5')
hdf5file = '/Users/tallen/Documents/2015_National_Hazard/2015_gsc_nshm/gm_tables/WinterfaceCombo_medclC.hdf5'
table_hdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], rhypo[1], vs30, hdf5file)
plt.loglog(array(table_hdf5imt['per']), exp(table_hdf5imt['sa']), 'k-', lw=2, label='2015 NBCC')

# plt OQ implementations
Yea97imt, AB03imt, Zea06imt, Zea06CISimt, AM09imt, MP10imt, GA14imt, GA14CISimt, Aea15imt \
    = interface_gsims(mag, dep, ztor, dip, rake, rrup[1], rjb[1], 760.)

plt.loglog(array(Zea06CISimt['per']), exp(Zea06CISimt['sa']), '--', c='m', lw=3., label='Zhao et al (2006) - OQ')
plt.loglog(array(AM09imt['per']), exp(AM09imt['sa']), '--', c='y', lw=3., label='Atkinson & Macias (2009) - OQ (B/C only)')


plt.xlabel('Period (sec)')
plt.ylabel('Spectral Acceleration (g)')
plt.title('Test 2020 NBCC Interface GMMs on Site Class B/C')
plt.grid(which='both')
plt.xlim([0.01, 10.])

#plt.legend(loc=3, fontsize=11)


plt.savefig('interface_hdf5_test.png', fmt='png', bbox_inches='tight')



plt.show()

