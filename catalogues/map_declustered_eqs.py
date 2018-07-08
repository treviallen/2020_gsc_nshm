from mpl_toolkits.basemap import Basemap
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from numpy import array, arange, mean, percentile, array, unique, where, ones_like, zeros_like
import matplotlib.pyplot as plt
from mapping_tools import drawoneshapepoly, distance, reckon
import shapefile
from os import path

suffix = 'swbc'
cnrs = [-132, -120, 46.5, 53] # SWBC

#suffix = 'arctic'
#cnrs = [-148, -118, 57, 69] # arctic

suffix = 'foothills'
cnrs = [-123, -112, 47, 57] # foothills

suffix = 'stlaurence'
cnrs = [-80, -60, 38, 50] # st laurence


mapres = 'i'
grdsize = 2.

fig = plt.figure(1, figsize=(20,10))

def make_basemap(cnrs):
    llcrnrlon = cnrs[0]
    urcrnrlon = cnrs[1]
    llcrnrlat = cnrs[2]
    urcrnrlat = cnrs[3]
    
    lon_0 = mean([llcrnrlon, urcrnrlon])
    lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
    lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
    
    plt.tick_params(labelsize=12)
    
    m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
                llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
                rsphere=6371200.,resolution=mapres,area_thresh=300)
    
    # draw coastlines, state and country boundaries, edge of map.
    #m.shadedrelief()
    m.drawcoastlines(linewidth=0.4) 
    m.drawstates()
    m.drawcountries()
    m.drawmapboundary(fill_color='0.8', zorder=100)
    m.fillcontinents(color='w',lake_color='0.8')
    m.drawparallels(arange(-90.,90.,grdsize), labels=[1,0,0,0],fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
    m.drawmeridians(arange(0.,360.,grdsize), labels=[0,0,0,1], fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)

    return m

#############################################################
# parse catalogue & plot
#############################################################

sheef_full = path.join('2010SHEEF','SHEEF2010Mw2_hmtk.csv')
sheef_decl = path.join('2010SHEEF','SHEEF2010Mw2_hmtk_declustered.csv')

parser = CsvCatalogueParser(sheef_full)
cat_full = parser.read_file()

lonf = cat_full.data['longitude']
latf = cat_full.data['latitude']
magf = cat_full.data['magnitude']

###############################################################
# plt full catalogue
###############################################################
plt.subplot(121)

# map earthquakes that pass completeness
m = make_basemap(cnrs)

# get index of events
idx = where((lonf >= cnrs[0]) & (lonf <= cnrs[1]) \
            & (latf >= cnrs[2]) & (latf <= cnrs[3]))[0]

for lo, la, ma in zip(lonf[idx], latf[idx], magf[idx]):
    ms = 2. * ma
    x, y = m(lo, la)
    m.plot(x, y, 'o', mec='k', mfc='none', mew=1.0, ms=ms)
    
plt.title('Full Catalogue', fontsize=18)

# make legend - set dummy x, y params
x, y= -100000, -100000

# get marker sizes
ms = 2. * 3.0
l1 = m.plot(x, y, 'ko', ms=ms)
ms = 2. * 5.0
l2 = m.plot(x, y, 'ko', ms=ms)
ms = 2. * 7.0
l3 = m.plot(x, y, 'ko', ms=ms)
plt.legend([l1[0], l2[0], l3[0]],['3.0', '5.0', '7.0'], fontsize=10, numpoints=1)

###############################################################
# plt declustered catalogue
###############################################################

parser = CsvCatalogueParser(sheef_decl)
cat_decl = parser.read_file()

lond = cat_decl.data['longitude']
latd = cat_decl.data['latitude']
magd = cat_decl.data['magnitude']

plt.subplot(122)

# map earthquakes that pass completeness
m = make_basemap(cnrs)

# get index of events
idx = where((lond >= cnrs[0]) & (lond <= cnrs[1]) \
            & (latd >= cnrs[2]) & (latd <= cnrs[3]))[0]

for lo, la, ma in zip(lond[idx], latd[idx], magd[idx]):
    ms = 2. * ma
    x, y = m(lo, la)
    m.plot(x, y, 'o', mec='k', mfc='none', mew=1.0, ms=ms)

plt.title('Declustered Catalogue', fontsize=18)        
###############################################################
# save figure
###############################################################
plt.savefig('map_declustering_eqs_'+suffix+'.png', fmt='png', bbox_inches='tight')
plt.show()