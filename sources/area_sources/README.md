To run for ONE source zone:

`run get_src_mfds.py <param file> <zone code>`

To run for ALL source zones:

`run get_src_mfds.py <param file>`

The parameter file (e.g. h_model.param) is of the following format:

- SHEEF FOLDER     # path to SHEEF-formatted files.  Expects 3 files: full, deep & crustal (i.e. shallow)
- DECLUSTERED      # tell the code if the catalogue is declustered
- SHAPEFILE INPUT  # path to zone shapefile
- OUTPUT FOLDER    # this is not currently used
- SHAPEFILE OUTPUT # name of output shapefile - saves to source model folder (e.g. 2020_H_model_YY-MM-DD)
- MFD BIN WIDTH    # MFD bin width

Note, the output shapefile should be a direct copy of the input shapefile, UNLESS any of the following are changed in the attribute table:

- source boundaries
- completeness periods
- Mmax/Mmin
- Min regression magnitude
- b-value & std are fixed

