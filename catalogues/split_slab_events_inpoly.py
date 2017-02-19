# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 11:30:47 2013

This version of the "split slab events" code only considers earthquakes 
located within the bounds of a given shapefile

Use deep source zones:
    
run split_slab_events_inpoly.py 2010SHEEF/SHEEF2010Mw2.0.full.gmtdat ../sources/area_sources/2020_H_model_2016-12-30/2020_H_model.shp

@author: tallen
"""
def checkfloat(floatstr):
    from numpy import nan
    try:
        return float(floatstr)
    except:
        return nan
        
'''
start main
'''
#try:        
# set filenames
from os import system, path
from sys import argv
from numpy import nan, isnan
from shapely.geometry import Point, Polygon
import shapefile
from gmt_tools import get_grd_extent

sheeffile = argv[1]
shpfile = argv[2]

# now parse sheef
print 'Reading SHEEF...'
data = open(sheeffile).readlines()

# get Cascadia & Alaska extents
casgrd = path.sep + path.join('Users', 'tallen', 'Documents', 'DATA', 'SLAB1.0', 'cas_slab1.0_clip.grd')
alugrd = path.sep + path.join('Users', 'tallen', 'Documents', 'DATA', 'SLAB1.0', 'alu_slab1.0_clip.grd')
xcas, ycas, zcas = get_grd_extent(grdfile=casgrd)
#xcas = 360 + xcas
xalu, yalu, zalu = get_grd_extent(grdfile=alugrd)

# set output files
slabfile = sheeffile.strip('full.gmtdat') + '.slab.gmtdat'
crustfile = sheeffile.strip('full.gmtdat') + '.crust.gmtdat'
f = open(slabfile, 'wb')
f.close()
f = open(crustfile, 'wb')
f.close()

# first parse shapefile and make shapely objects
print 'Reading shapefile...'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))

# now loop through events
print 'Looping thru events...'
for line in data:
    # reset variables
    inpoly = False
    slabdiff = nan 
    slabdep = nan

    lat = line[31:37]
    lon = line[20:28].strip()
    
    # convert lon for Alaska grid file
    if float(lon) < 0.0:
        alon = str(360 + float(lon))
        
    # now loop through source zone polygons to see if we need to bother
    pt = Point(float(lon), float(lat))
    for poly in polygons:
        if pt.within(poly):
            inpoly = True
        
    # if the point is in the regions of interest, continue
    if inpoly == True:
        # check depth
        if len(line) < 75:
            eqdep = line[44:51].strip()
            fixdep = line[52:53].strip()
        else:
            eqdep = line[47:53].strip()
            fixdep = line[54:55].strip()
            
        # call gmt grdtrack on Cascadia
        if float(lon) >= xcas[0] and float(lon) <= xcas[1] \
           and float(lat) >= ycas[0] and float(lat) <= ycas[1]:
           
           # write coords to temp file
           f = open('temp.xy','wb')
           f.write('\t'.join((lon, lat)))
           f.close()
           
           system('gmt grdtrack temp.xy -G'+casgrd+' > temp.xyz')
           slabstr = open('temp.xyz').read().strip().split('\t')[-1]
           slabdep = abs(float(slabstr))
           
        # call gmt grdtrack on Alaska
        elif float(alon) >= xalu[0] and float(alon) <= xalu[1] \
           and float(lat) >= yalu[0] and float(lat) <= yalu[1]:
           
           # write coords to temp file       
           f = open('temp.xy','wb')      
           f.write('\t'.join((alon, lat)))
           f.close()                     
           
           system('gmt grdtrack temp.xy -G'+alugrd+' > temp.xyz')
           slabstr = open('temp.xyz').read().strip().split('\t')[-1]
           slabdep = abs(checkfloat(slabstr))
              
        # check to see if related to slab or crust
        if not eqdep.strip() == '':
            slabdiff = slabdep - float(eqdep)
        else:
            slabdiff = nan
        
    # now prepare to write files
    #print 'Writing to files...'
    if inpoly == True:
        #if isnan(slabdiff) == False:        
            if checkfloat(eqdep) >= 25.0 or slabdiff <= 10.0:
                # now write files
                try:
                    fs.write(line)
                except:
                    fs = open(slabfile, 'a')
                    fs.write(line)
            # do crust events
            else:
                try:
                    fc.write(line)
                except:
                    fc = open(crustfile, 'a')
                    fc.write(line)
                    '''
        # do nan events
        else:
            try:
                fc.write(line)
            except:
                fc = open(crustfile, 'a')
                fc.write(line)
                '''
    else:
        try:
            fc.write(line)
        except:
            fc = open(crustfile, 'a')
            fc.write(line)
   
fs.close()
fc.close()

#except:
#    print '\nUsage: python split_slab_events_inpoly.py <sheeffile> <shpfile>\n'
