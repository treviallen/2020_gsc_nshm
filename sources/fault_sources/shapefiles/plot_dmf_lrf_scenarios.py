from mpl_toolkits.basemap import Basemap
from numpy import array, arange, mean, percentile, array, unique, where, ones_like, zeros_like, sin, radians
import matplotlib.pyplot as plt
from mapping_tools import drawoneshapepoly, distance, reckon
import shapefile
from mapping_tools import get_field_data
from fault_tools import *


cnrs = [-124.4, -121.6, 48, 49]
mapres = 'h'
grdsize = .4

# raed shp
shpfile = '//Users//tallen//Documents//2020_National_Hazard//2020_gsc_nshm//sources//fault_sources//shapefiles//DMF_LRF.shp'
sf = shapefile.Reader(shpfile)
fcode = get_field_data(sf, 'CODE', 'str')
flen  = get_field_data(sf, 'LENGTH', 'float')

seis_thickness = 15.
dip = 70.


fig = plt.figure(1, figsize=(15,15))


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

# loop thru 7 scenarios - ignore full rupture
for i in range(1,8):
    plt.subplot(4,2,i)


    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.drawmapboundary(fill_color='0.8', zorder=100)
    m.fillcontinents(color='w',lake_color='0.8')
    m.drawparallels(arange(-90.,90.,grdsize), labels=[1,0,0,0],fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
    m.drawmeridians(arange(0.,360.,grdsize*2), labels=[0,0,0,1], fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
    
    # plt fault
    points = sf.shapes()[i].points
    
    la = []
    lo = []
    for p in points:
        lo.append(p[0])
        la.append(p[1])
    
    x, y = m(lo, la)
    
    m.plot(x, y, 'r-', lw=2)
    
    # get fault stats
    ftype = 'ss'
    area = flen[i] * seis_thickness / sin(radians(dip))
    mchar = return_area2mag_for_ss(area)[1] - 0.25
    
    # add title
    plt.title('; '.join((fcode[i], 'L = '+str('%0.0f' % flen[i])+' km', 'Mc = '+str('%0.2f' % mchar))), fontsize=14)
    

"""
###############################################################
# schematic evs
evla = [49, 47] 
evlo = [-130, -125]
cols = ['g', 'b']
xoff = [-0.1, 0.1]

minperpdist = ones_like(evla) * 99999
minparldist = zeros_like(evla)
cumdist = 0

# loop thru records and get evdist to boundary
for k in range(0, len(evla)):
    pla = []
    plo = []
    
    for i in range(0, len(points)-1):
        # get rng and az between points
        fault_rngkm, fault_az, baz = distance(points[i][1], points[i][0], points[i+1][1], points[i+1][0])
        
        # now discretize between points
        for j in range(0, 51):
            fdist = j * fault_rngkm / 50.
            cumdist += fault_rngkm / 50.
            
            lo, la = reckon(points[i][1], points[i][0], fdist, fault_az)
            
            rngkm, az, baz = distance(la, lo, evla[k], evlo[k])
            
            # fill mindist
            if rngkm < abs(minperpdist[k]):
            	
                # get diff in azimuth
                az_diff = fault_az - az
                
                if az_diff >= 0 and az_diff <= 180:
                    minperpdist[k] = -1 * rngkm
                else:
                    minperpdist[k] = rngkm
                
                minparldist[k] = cumdist
                
                plo.append(lo)
                pla.append(la)
            
    x, y = m(array(plo)+xoff[k], array(pla))
    m.plot(x, y, '--', c=cols[k], lw=3)
    
    x, y = m(evlo[k], evla[k])
    m.plot(x, y, 'o', c=cols[k], ms=20)
    
    # link eqs to along trench lines
    llo = [evlo[k], plo[-1]+xoff[k]]
    lla = [evla[k], pla[-1]]
    x, y = m(array(llo), array(lla))
    m.plot(x, y, '--', c=cols[k], lw=3)

"""
###############################################################
plt.savefig('dmf_lrf_scenarios.png', fmt='png', bbox_inches='tight')
plt.show()