from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from hmtk.seismicity.declusterer.dec_gardner_knopoff import GardnerKnopoffType1
from hmtk.seismicity.declusterer.distance_time_windows import GardnerKnopoffWindow, GruenthalWindow
from io_catalogues import sheef2hmtk_csv
from os import path, remove
from copy import deepcopy
from numpy import isnan

'''
got basics from here:

http://seiscode.iag.usp.br/gitlab/hazard/pshab_source_models/raw/646202c6c5a38426783b4851b188280a1441e032/notes/01_workflow_decluster.py
'''

#####################################################
# parse catalogues and prep declusterer
#####################################################

# reformat SHEEF
#sheef2hmtk_csv(path.join('2010SHEEF','SHEEF2010Mw2.0.full.gmtdat')) # only need to do once

# parse HMTK catalogue
parser = CsvCatalogueParser(path.join('2010SHEEF','SHEEF2010Mw2_hmtk.csv'))
catalogue = parser.read_file()

decluster_config = {'time_distance_window': GruenthalWindow(),
                    'fs_time_prop': 1.0}

#####################################################
# decluster here
#####################################################

print 'Running GK declustering...'
decluster_method = GardnerKnopoffType1()

#---------------------------------------------
# declustering
cluster_index_gk, cluster_flags_gk = decluster_method.decluster(catalogue, decluster_config)
#---------------------------------------------

# adding to the catalog
# The cluster flag (main shock or after/foreshock) and cluster index to the catalogue keys
catalogue.data['cluster_index_gk'] = cluster_index_gk
catalogue.data['cluster_flags_gk'] = cluster_flags_gk

#####################################################
# purge remove non-poissonian events
#####################################################

# create a copy from the catalogue object to preserve it
catalogue_gk = deepcopy(catalogue)

catalogue_gk.purge_catalogue(cluster_flags_gk == 0) # cluster_flags == 0: mainshocks

print 'Gardner-Knopoff\tbefore: ', catalogue.get_number_events(), " after: ", catalogue_gk.get_number_events()

#####################################################
# write declustered catalogue
#####################################################

# setup the writer
declustered_catalog_file = path.join('2010SHEEF','SHEEF2010Mw2_hmtk_declustered.csv')

# if it exists, delete previous file
remove(declustered_catalog_file)

# set-up writer
writer = CsvCatalogueWriter(declustered_catalog_file) 

# write
writer.write_file(catalogue_gk)
#writer.write_file(catalogue_af)
print 'Declustered catalogue: ok\n'

#####################################################
# write hmtk catalogue to SHEEF fmt
#####################################################

print 'Writing to SHEEF...'
newsheef = ''
cat = catalogue_gk.data #just shorten variable name
for i in range(0, catalogue_gk.get_number_events()):
    mag = str('%0.1f' % cat['magnitude'][i])
    
    lat = '  ' + str('%0.3f' % cat['latitude'][i])
    
    lon = str('%0.3f' % cat['longitude'][i])
    if cat['longitude'][i] > -100.:
        lon = ' ' + lon
    
    dep = str('%0.2f' % cat['depth'][i])
    if isnan(cat['depth'][i]):
        dep = '         '
    elif cat['depth'][i] < 10.:
        dep = '     ' + dep
    elif cat['depth'][i] < 100.:
        dep = '    ' + dep
    else:
        dep = '   ' + dep
    
    src = '  ' + cat['Agency'][i]
    if len(src) == 3:
        src += '  '
    elif len(src) == 3:
        src += ' '
        
    mtype = cat['magnitudeType'][i]
    if len(mtype) == 2:
        mtype += ' '
        
    mag2 = str('%0.2f' % cat['magnitude'][i])
    
    newsheef += ' '.join(('', str(cat['eventID'][i]), '', mag, '', lon, lat, src, dep, ' ', mtype, '', mag, 'UK    ', mag2)) + '\n'
    
newsheeffile = declustered_catalog_file.strip('_hmtk_declustered.csv') + '_declustered.gmtdat'
f = open(newsheeffile, 'wb')
f.write(newsheef)
f.close()