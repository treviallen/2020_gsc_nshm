'''
Code makes table text for various GMMs using OQ-hazardlib implementation

Various updates are made to the models in the process:
	- extrapolation of long-period accelerations to 10s
	- addition of PGV for models with no PGV parameterisation
	- application of SS14/AB06 amplification factors for models with no Vs30 parameterisation
	- depth-specific tables produced
	
'''

from calc_oq_gmpes import gsim2table # from https://github.com/treviallen/my_codes/blob/master/calc_oq_gmpes.py
from numpy import logspace, arange, array

dists = logspace(0,3, 31)
vs30rng = [115, 250, 450, 760., 1100, 1600] # in m/s
extrapPeriod = [5., 10.]
folder = 'gmm_txt_tables' # output folder

# function to make table
def make_tables(gmmClass, gmmName, vs30rng, depths, modVs30, vs30ref, interpPeriods):
    for vs30 in vs30rng:
        for depth in depths:
            
            # if vs30 parameterisation
            if modVs30 == True:
                # use AB06 anyway if vs30 > Vs30Ref
                if vs30 > vs30ref:
                    vs30refPass = vs30ref
                else:
                    vs30refPass = vs30
            
            # else if no vs30 parameterisation
            else:
                vs30refPass = vs30ref
            
            tabtxt = gsim2table(gmmClass, gmmName, mags, dists, depth, vs30, vs30refPass, extrapPeriod, interpPeriods, rtype, folder)



'''
# make slab models
'''
mags = arange(4.75, 8.1, 0.25)
"""
from openquake.hazardlib.gsim.garcia_2005 import GarciaEtAl2005SSlab
gmmClass = GarciaEtAl2005SSlab()
gmmName = 'GarciaEtAl2005SSlab'
rtype = 'rhypo'
modVs30 = False
vs30ref = 1100.
depths = [20., 30., 50.]
#depths = [50.]
interpPeriods = False
make_tables(gmmClass, gmmName, vs30rng, depths, modVs30, vs30ref, interpPeriods)
"""
# no PGV
from openquake.hazardlib.gsim.atkinson_boore_2003 import AtkinsonBoore2003SSlabCascadia
gmmClass = AtkinsonBoore2003SSlabCascadia()
gmmName = 'AtkinsonBoore2003SSlabCascadia'
modVs30 = True
vs30ref = 760.
rtype = 'rrup'
depths = [20., 30., 50.]
interpPeriods = False
make_tables(gmmClass, gmmName, vs30rng, depths, modVs30, vs30ref, interpPeriods)

# no PGV
from openquake.hazardlib.gsim.zhao_2006 import ZhaoEtAl2006SSlabCascadia
gmmClass = ZhaoEtAl2006SSlabCascadia()
gmmName = 'ZhaoEtAl2006SSlabCascadia'
modVs30 = True
vs30ref = 760.
depths = [20., 30., 50.]
rtype = 'rrup'
interpPeriods = False
make_tables(gmmClass, gmmName, vs30rng, depths, modVs30, vs30ref, interpPeriods)

from openquake.hazardlib.gsim.abrahamson_2015 import AbrahamsonEtAl2015SSlab
gmmClass = AbrahamsonEtAl2015SSlab()
gmmName = 'AbrahamsonEtAl2015SSlab'
modVs30 = True
vs30ref = 1000.
depths = [20., 30., 50.]
rtype = 'rhypo'
interpPeriods = False
make_tables(gmmClass, gmmName, vs30rng, depths, modVs30, vs30ref, interpPeriods)


'''
# make interface models
'''
mags = arange(4.75, 9.3, 0.25) # go to 9.25
"""
from openquake.hazardlib.gsim.atkinson_macias_2009 import AtkinsonMacias2009
gmmClass = AtkinsonMacias2009()
gmmName = 'AtkinsonMacias2009'
rtype = 'rrup'
modVs30 = False
vs30ref = 760.
depths = [30.]
interpPeriods = False
make_tables(gmmClass, gmmName, vs30rng, depths, modVs30, vs30ref, interpPeriods)



from openquake.hazardlib.gsim.ghofrani_atkinson_2014 import GhofraniAtkinson2014, GhofraniAtkinson2014Cascadia
gmmClass = GhofraniAtkinson2014Cascadia()
gmmName = 'GhofraniAtkinson2014Cascadia'
modVs30 = True
vs30ref = 760.
rtype = 'rrup'
depths = [30.]
interpPeriods = False
make_tables(gmmClass, gmmName, vs30rng, depths, modVs30, vs30ref, interpPeriods)


from openquake.hazardlib.gsim.zhao_2006 import ZhaoEtAl2006SInterCascadia
gmmClass = ZhaoEtAl2006SInterCascadia()
gmmName = 'ZhaoEtAl2006SInterCascadia'
modVs30 = True
vs30ref = 760.
depths = [30.]
rtype = 'rrup'
interpPeriods = False
make_tables(gmmClass, gmmName, vs30rng, depths, modVs30, vs30ref, interpPeriods)

from openquake.hazardlib.gsim.abrahamson_2015 import AbrahamsonEtAl2015SInter
gmmClass = AbrahamsonEtAl2015SInter()
gmmName = 'AbrahamsonEtAl2015SInter'
modVs30 = True
vs30ref = 1000.
depths = [30.]
rtype = 'rrup'
interpPeriods = False
make_tables(gmmClass, gmmName, vs30rng, depths, modVs30, vs30ref, interpPeriods)


"""