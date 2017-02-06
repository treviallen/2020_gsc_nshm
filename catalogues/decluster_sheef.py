from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from hmtk.seismicity.declusterer.dec_gardner_knopoff import GardnerKnopoffType1
from hmtk.seismicity.declusterer.distance_time_windows import GardnerKnopoffWindow, GruenthalWindow
from io_catalogues import sheef2hmtk_csv
from os import path
from copy import deepcopy

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
writer = CsvCatalogueWriter(declustered_catalog_file) 

# write
writer.write_file(catalogue_gk)
#writer.write_file(catalogue_af)
print 'Declustered catalogue: ok\n'
