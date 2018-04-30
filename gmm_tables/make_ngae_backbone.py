from misc_tools import listdir_extension
from os import system, path, remove
from numpy import zeros, log10, savetxt, zeros_like, interp, log, log10, array

# set variables

nmag = 11
ndist = 34
nperiod = 25
nrows = nmag * ndist
ncols = nperiod+2

altSigma = False
siteC = True

wght_sum = zeros((nrows, ncols))

# get pre-weighted gmm files
folder =  'weighted_usgs_ngae_tables_450' 
extension = '450mps.txt'
gmpetables = listdir_extension(folder, extension)

# read tables and append data
for table in gmpetables:
    tablepath = path.join(folder, table)
    lines = open(tablepath).readlines()[4:]
    
    for i, line in enumerate(lines):
        dat = line.strip().split()
        
        # first fill mag and dist - overwrite is ugly but does the job
        wght_sum[i, 0] = float(dat[0])
        wght_sum[i, 1] = float(dat[1])
        
        # sum amplitude data - already weighted, so just need to add
        for j, d in enumerate(dat[2:]):
            #wght_sum[i, j+2] += 10**float(d) # if arithmetic
            wght_sum[i, j+2] += float(d) # if geometric

# take log of amps - if arithmetic
#wght_sum[:, 2:] = log10(wght_sum[:, 2:])
	
# now put back in table
if altSigma == True: # use TA alternate sigma model
    headers = open(tablepath).readlines()[0:3]
    	
    ####################################################################
    #  crude NGA-E sigma
    ####################################################################
    
    # get periods
    pertxt = headers[-1].strip().split()[:-2]
    periods = array([float(x) for x in pertxt])
    
    shortT = 0.01
    midT = 0.1
    longT  = 10.0
    
    shortTsig = 0.70
    midTsig = 0.81
    longTsig = 0.57
    
    ngaesig = zeros_like(periods)
    
    idx = periods <= midT
    ngaesig[idx] = interp(log10(periods[idx][::-1]), [log10(shortT), log10(midT)], [shortTsig, midTsig])[::-1]
    
    idx = periods > midT
    ngaesig[idx] = interp(log10(periods[idx][::-1]), [log10(midT), log10(longT)], [midTsig, longTsig])[::-1]
    	
    ngaesigtxt = '          ' + ' '.join([str('%0.3f' % x) for x in ngaesig]) + ' 0.700 0.810'
    headers += ngaesigtxt

    
else: # use AA13 sigma
    headers = open(tablepath).readlines()[0:4] 
header = ''
for h in headers:
    header += h

# write file
#outfile = 'NGA-East_Backbone_Model.arithmetic.txt'
if altSigma == True:
    outfile = 'NGA-East_Backbone_Model.geometric.altSigma.txt'
else:
    outfile = 'NGA-East_Backbone_Model.geometric.450mps.txt'
savetxt(outfile, wght_sum, fmt='%.3f', delimiter=' ', header=header, comments='')

print '\n!!!!MAKE MINOR EDITS TO BACKBONE HEADER!!!!\n'

	