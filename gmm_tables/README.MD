Steps for generating hdf5 tables
================================

1. Run "gsim2table.py":
	- Calls [oq-engine hazardlib](https://github.com/gem/oq-engine/tree/master/openquake/hazardlib) functions and outputs gmround motion text tables for a given Vs30
	- Outputs to the folder "gmm_txt_tables"

* Run "gmpe_txt2hdf5.py":
	- Converts text tables to hdf5 format for use with the oq-engine
	- Outputs to the folder "gmm_hdf5_tables"
	- Useage:
		> python gmpe_txt2hdf5.py [in folder] [out folder]
		
* Test hdf5 tables by running "plt_interface_hdf5_tables.py":
	- calls hazardlib function "openquake.hazardlib.gsim.gmpe_table.GMPETable"
	- source and site parameters hardwired in code