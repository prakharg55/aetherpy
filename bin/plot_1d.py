#!/usr/bin/env python
# Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
"""Block-based model visualization routines."""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from glob import glob
import os
import re
import argparse

from aetherpy import logger
from aetherpy.io import read_routines
from aetherpy.plot import data_prep
from aetherpy.plot import movie_routines
from aetherpy.utils import inputs
from aetherpy.utils import time_conversion

# ----------------------------------------------------------------------------
# parse arguments and get list of files
# ----------------------------------------------------------------------------

def get_args():

    parser = argparse.ArgumentParser(
        description = 'Plot Aether / GITM model results')

    parser.add_argument('-IsLog', default =False,
                        help='plot the log of the variable', 
                        action="store_true")
    
    parser.add_argument('-list', default =False,
                        help='list variable names', 
                        action="store_true")    
        
    parser.add_argument('-lat', metavar = 'lat',  default = 0.0, type = float, \
                        help = 'latitude : latitude in degrees (closest)')
    parser.add_argument('-lon', metavar = 'lon',  default = 0.0, type = float, \
                        help = 'longitude in degrees (closest)')

    parser.add_argument('-var',  \
                        default = 3, type = int, \
                        help = 'variable to plot (number)')
    
    parser.add_argument('filelist', nargs='+', \
                        help = 'list files to use for generating plots')
    
    args = parser.parse_args()

    return args


# ----------------------------------------------------------------------------
# change variable names and remove some weird characters
# ----------------------------------------------------------------------------

def fix_vars(vars):
    newvars = []
    for v in vars:
        nv = re.sub('!U', '', v)
        nv = re.sub('!N', '', nv)
        nv = re.sub('!D', '', nv)
        newvars.append(nv)

    return newvars

# ----------------------------------------------------------------------------
# Determine if these are GITM or Aether files
# ----------------------------------------------------------------------------

def determine_file_type(file):

    IsGitm = False
    HasHeader = False
    m = re.match(r'(.*)bin', file)
    if m:
        IsGitm = True
        # check for a header file:
        checkFile = glob(m.group(1)+"header")
        if (len(checkFile) > 0):
            if (len(checkFile[0]) > 1):
                HasHeader = True

    return IsGitm, HasHeader



# Needed to run main script as the default executable from the command line
if __name__ == '__main__':

    # Get the input arguments
    args = get_args()

    # determine what kind of files we are dealing with
    IsGitm, HasHeader = determine_file_type(args.filelist[0])

    # assume Aether NC files....
    header = read_routines.read_aether_headers(args.filelist)
    header['vars'] = fix_vars(header['vars'])
    
    # Read headers for input files (assumes all files have same header)
    header = read_routines.read_blocked_netcdf_header(args.filelist[0])

    if (args.list):
        for i, v in enumerate(header['vars']):
            print(i, '-->', v)
        exit()
        
    if ('z' in header['vars']):
        alt_var = 'z'
    else:
        alt_var = 'alt'
    
    # Determine variables to plot (currently hardcoded)
    # TODO: handle winds correctly

    loc_vars = ['lon', 'lat', alt_var]
    data = read_routines.read_blocked_netcdf_file(args.filelist[0], loc_vars)

    # these are block arrays, so first element is block number
    alts = data[alt_var][0, 0, 0, :] / 1000.0  # Convert from m to km
    lons = data['lon'][0, :, 0, 0]
    lats = data['lat'][0, 0, :, 0]

    print(args.lon)
    
    dLon = np.abs(lons - args.lon)
    iLon = np.argmin(dLon)
    dLat = np.abs(lats - args.lat)
    iLat = np.argmin(dLat)

    print(iLon, iLat)
    print('Closest Lon : ', lons[iLon])
    print('Closest lat : ', lats[iLat])

    varName = header['vars'][args.var]
    var_list = [varName]

    fig = plt.figure(figsize=(10, 8.5))
    ax = fig.add_subplot(111)
        
    for filename in args.filelist:
        # Retrieve data
        data = read_routines.read_blocked_netcdf_file(filename, var_list)

        x = data[varName][0, iLon, iLat, :]
        y = alts
        ax.plot(x, y)
        if (args.IsLog):
            ax.set_xscale('log')

    fig.savefig('test.png')
        
        
    #file_vars = ['lon', 'lat', alt_var, args['var']] if args['var'] else None
