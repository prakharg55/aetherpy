#!/usr/bin/env python
# Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
"""Routines to read Aether files."""

import datetime as dt
from netCDF4 import Dataset
import numpy as np
from struct import unpack
import re

from aetherpy.utils.time_conversion import epoch_to_datetime
from aetherpy import logger


def read_aether_headers(filelist, finds=-1):
    """Obtain ancillary information from the netCDF file.

    Parameters
    ----------
    filelist : array-like
        Array-like object of names for netCDF files
    finds : int, list, or slice
        Index(es) for file(s) from which the header information will be read.
        (default=-1)

    Returns
    -------
    header : dict
        A dictionary containing header information from the netCDF files,
        including:
        nfiles - number of files in the list (not the number read)
        nlons - number of longitude grids
        nlats - number of latitude grids
        nalts - number of altitude grids
        vars - list of data variable names
        time - list of datetimes for the processed file start times
        filename - list of the processed filenames

    Raises
    ------
    IOError
        If there are an unexpected number or name of variables or dimensions
        for any file in the file list.

    Notes
    -----
    This routine obtains the same info as `read_gitm_headers`

    """
    # Initialize the output
    header = {"nfiles": len(filelist), "version": 0.1, "nlons": 0, "nlats": 0,
              "nalts": 0, "vars": [], "time": []}

    # Ensure the filelist is array-like, allowing slicing of input
    file_list = np.asarray(filelist)
    header['filename'] = file_list

    # Read the header info from the desired files
    for filename in header['filename']:
        with Dataset(filename, 'r') as ncfile:
            ncvars = list()
            for var in ncfile.variables.values():
                if len(var.shape) == 3:
                    nlons = var.shape[0]
                    nlats = var.shape[1]
                    nalts = var.shape[2]
                    ncvars.append(var.name)

                    # Test the dimensions
                    if(header['nlons'] == 0 and header['nlats'] == 0
                       and header['nalts'] == 0):
                        header["nlons"] = nlons
                        header["nlats"] = nlats
                        header["nalts"] = nalts
                    elif(header['nlons'] != nlons or header['nlats'] != nlats
                         or header['nalts'] != nalts):
                        raise IOError(''.join(['unexpected dimensions for ',
                                               'variable ', var.name,
                                               ' in file ', filename]))

            # Test the variable names
            if len(header["vars"]) == 0:
                header["vars"] = list(ncvars)
            elif header["vars"] != ncvars:
                raise IOError(''.join(['unexpected number or name of ',
                                       'variables in file ', filename]))

            # Add the time for this file
            t = np.double(ncfile.variables['Time'][0])
            header["time"].append(epoch_to_datetime(t))

    return header


def parse_line_into_int_and_string(line):
    """Parse a data string into integer and string components.

    Parameters
    ----------
    line : str
        Data line with a known format

    Returns
    -------
    line_num : int
        Integer corresponding to the first number in `line`
    line_str : str
        Whitespace separated string containing the rest of the `line` data

    Notes
    -----
    Splits data line using a single whitespace.  The first element is
    expected to be an integer and the remaining data is returned as a
    single whitespace separated string.

    """

    split_line = line.split(" ")
    line_num = int(split_line[0])
    line_str = " ".join(split_line[1:])

    return line_num, line_str


def read_aether_ascii_header(filelist):
    """Read information from the ascii header file, which mirrors GITM.

    Parameters
    ----------
    filelist : list
        A list of ascii header file names.  The number of files is recorded,
        but only the last file is used.

    Returns
    -------
    header : dict
        A dictionary containing header information from the netCDF files,
        including:
        nfiles - number of files in list
        nlons - number of longitude grids
        nlats - number of latitude grids
        nalts - number of altitude grids
        nvars - number of data variable names
        time - list of datetimes with data

    Notes
    -----


    """

    header = {"nfiles": len(filelist), "version": 0.1, "nlons": 0, "nlats": 0,
              "nalts": 0, "nblockslons": 0, "nblockslats": 0, "nblocksalts": 0,
              "nvars": 0, "vars": [], "time": [], "filename": []}

    for filename in filelist:
        header["filename"].append(filename)
        headerfile = filename.replace(".bin", ".header")

        with open(headerfile, 'r') as fpin:
            for line in fpin:
                if re.match(r'BLOCKS', line):
                    header["nblockslons"] = int(fpin.readline())
                    header["nblockslats"] = int(fpin.readline())
                    header["nblocksalts"] = int(fpin.readline())

                if re.match(r'NUMERICAL', line):
                    header["nvars"], _ = parse_line_into_int_and_string(
                        fpin.readline())
                    header["nlons"], _ = parse_line_into_int_and_string(
                        fpin.readline())
                    header["nlats"], _ = parse_line_into_int_and_string(
                        fpin.readline())
                    header["nalts"], _ = parse_line_into_int_and_string(
                        fpin.readline())

                if re.match(r'VERSION', line):
                    header["version"] = float(fpin.readline())

                if re.match(r'TIME', line):
                    year = int(fpin.readline())
                    month = int(fpin.readline())
                    day = int(fpin.readline())
                    hour = int(fpin.readline())
                    minute = int(fpin.readline())
                    second = int(fpin.readline())
                    msec = int(fpin.readline())
                    header["time"].append(dt.datetime(year, month, day, hour,
                                                      minute, second, msec))

                mvar = re.match(r'VARIABLE', line)
                if mvar and len(header["vars"]) < 1:
                    for i in range(header["nvars"]):
                        nlin, sline = parse_line_into_int_and_string(
                            fpin.readline())
                        header["vars"].append(sline.strip())

    return header


def read_aether_one_binary_file(header, ifile, vars_to_read):
    """Read in list of variables from a single netCDF file.

    Parameters
    ----------
    header : dict
        Dict of header data returned from read_aether_ascii_header
    ifile : int
        Integer corresponding to filename in the header dict
    vars_to_read : list
        List of desired variable names to read

    Returns
    -------
    data : dict
        Dict with keys 'time', which contains a datetime object specifying the
        time of the file and indices ranging from 0 to `len(vars_to_read) - 1`,
        corresponding to the variable names in `vars_to_read` that holds arrays
        of the specified data

    See Also
    --------
    read_aether_ascii_header

    """

    file_to_read = header["filename"][ifile]
    logger.info("Reading file : ", file_to_read)

    data = {hkey: header[hkey] for hkey in ["version", "nlons", "nlats",
                                            "nalts" "nvars", "vars"]}
    data["time"] = header["time"][ifile]

    with open(file_to_read, 'rb') as fin:
        # Read in the header data
        header_len = 0
        num_tot = data["nlons"] * data["nlats"] * data["nalts"]

        # Data were saved as floats, not doubles
        data_len = num_tot * 4
        endChar = '<'
        for ivar in vars_to_read:
            fin.seek(header_len + ivar * data_len)
            data[ivar] = np.array(unpack(endChar + '%if' % num_tot,
                                         fin.read(data_len)))
            data[ivar] = data[ivar].reshape(
                (data["nlons"], data["nlats"], data["nalts"]), order="F")

    return data


def read_gitm_headers(filelist, finds=-1):
    """Read ancillary information from a GITM file.

    Parameters
    ----------
    filelist : array-like
        Array-like object of names for netCDF files
    finds : int, list, or slice
        Index(es) for file(s) from which the header information will be read.
        (default=-1)

    Returns
    -------
    header : dict
        A dictionary containing header information from the netCDF files,
        including:
        nfiles - number of files in list
        nlons - number of longitude grids
        nlats - number of latitude grids
        nalts - number of altitude grids
        vars - list of data variable names
        time - list of datetimes for the processed file start times
        filename - list of the processed filenames

    Raises
    ------
    IOError
        If any one of the input files encounters an unexpected value or
        dimension.

    Notes
    -----
    This routine obtains the same info as `read_aether_headers`

    """
    # Initialize the output
    header = {"nfiles": len(filelist), "version": 0, "nlons": 0, "nlats": 0,
              "nalts": 0, "vars": [], "time": []}

    # Ensure the filelist is array-like, allowing slicing of input
    file_list = np.asarray(filelist)
    header['filename'] = [file_list[finds]]
    for filename in header['filename']:
        # Read in the header from the binary file
        file_vars = list()

        with open(filename, 'rb') as fin:
            # Test to see if the correct endian is being used
            end_char = '>'
            raw_rec_len = fin.read(4)
            rec_len = (unpack(end_char + 'l', raw_rec_len))[0]
            if rec_len > 10000 or rec_len < 0:
                # Ridiculous record length implies wrong endian, fix here
                end_char = '<'
                rec_len = (unpack(end_char + 'l', raw_rec_len))[0]

            # Read version; read fortran footer+header.
            file_version = unpack(end_char + 'd', fin.read(rec_len))[0]
            _, rec_len = unpack(end_char + '2l', fin.read(8))

            # Test the version number
            if header['version'] == 0:
                header["version"] = file_version
            elif header['version'] != file_version:
                raise IOError(''.join(['unexpected version number in file ',
                                       filename]))

            # Read grid size information.
            nlons, nlats, nalts = unpack(end_char + 'lll', fin.read(rec_len))
            _, rec_len = unpack(end_char + '2l', fin.read(8))

            # Test the dimensions
            if(header['nlons'] == 0 and header['nlats'] == 0
               and header['nalts'] == 0):
                header["nlons"] = nlons
                header["nlats"] = nlats
                header["nalts"] = nalts
            elif(header['nlons'] != nlons or header['nlats'] != nlats
                 or header['nalts'] != nalts):
                raise IOError(''.join(['unexpected dimensions in file ',
                                       filename]))

            # Read number of variables.
            num_vars = unpack(end_char + 'l', fin.read(rec_len))[0]
            _, rec_len = unpack(end_char + '2l', fin.read(8))

            # Collect variable names.
            for ivar in range(num_vars):
                vcode = unpack(end_char + '%is' % (rec_len),
                               fin.read(rec_len))[0]
                var = vcode.decode('utf-8').replace(" ", "")
                file_vars.append(var)
                _, rec_len = unpack(end_char + '2l', fin.read(8))

            # Test the variable names
            if len(header["vars"]) == 0:
                header["vars"] = list(file_vars)
            elif header["vars"] != file_vars:
                raise IOError(''.join(['unexpected number or name of ',
                                       'variables in file ', filename]))

            # Extract time
            out_time = np.array(unpack(end_char + 'lllllll', fin.read(rec_len)))
            out_time[-1] *= 1000  # Convert from millisec to microsec
            header["time"].append(dt.datetime(*out_time))

    return header


def read_aether_file(filename, file_vars=None):
    """Read in list of variables from a netCDF file.

    Parameters
    ----------
    filename : str
        Name of netCDF file to read
    file_vars : list or NoneType
        List of desired variable names to read, or None to read all
        (default=None)

    Returns
    -------
    data : dict
        Dict with keys 'time', which contains a datetime object specifying the
        time of the file and zero-offset indices, corresponding to the
        variable names in `file_vars` that holds arrays of the specified data.
        Also contains a list of the variable names.

    """
    data = dict()

    with Dataset(filename, 'r') as ncfile:
        # Save a list of all variable names
        data['vars'] = [var for var in ncfile.variables.keys()]

        # If file variables are not provided, set them here
        if file_vars is None:
            file_vars = data['vars']

        # Save the data as numpy arrays, using variable index as a key
        for i_var, var in enumerate(file_vars):
            data[i_var] = np.array(ncfile.variables[var])

        # Calculate the date and time for this data
        time = np.array(ncfile.variables['Time'])
        data['time'] = epoch_to_datetime(time[0])

    return data


def read_gitm_file(filename, file_vars=None):
    """Read list of variables from one GITM file.

    Parameters
    ----------
    filename : str
        GITM file to read
    file_vars : list or NoneType
        List of desired variable names to read or None to read all
        (default=None)

    Returns
    -------
    data : dict
        Dict with keys 'time', which contains a datetime object specifying the
        time of the file and zero-offset indices, corresponding to the
        variable names in `file_vars` that holds arrays of the specified data.
        Also contains version number, dimensions, and a list of the variable
        names.

    """

    data = {"version": 0, "nlons": 0, "nlats": 0, "nalts": 0, "time": 0,
            "vars": []}

    with open(filename, 'rb') as fin:
        # Determine the correct endian
        end_char = '>'
        raw_rec_len = fin.read(4)
        rec_len = (unpack(end_char + 'l', raw_rec_len))[0]
        if rec_len > 10000 or rec_len < 0:
            # Ridiculous record length implies wrong endian.
            end_char = '<'
            rec_len = (unpack(end_char + 'l', raw_rec_len))[0]

        # Read version; read fortran footer+data.
        data["version"] = unpack(end_char + 'd', fin.read(rec_len))[0]

        _, rec_len = unpack(end_char + '2l', fin.read(8))

        # Read grid size information.
        data["nlons"], data["nlats"], data["nalts"] = unpack(
            end_char + 'lll', fin.read(rec_len))
        _, rec_len = unpack(end_char + '2l', fin.read(8))

        # Read number of variables.
        num_vars = unpack(end_char + 'l', fin.read(rec_len))[0]
        _, rec_len = unpack(end_char + '2l', fin.read(8))

        if file_vars is None:
            file_vars = np.arange(0, num_vars, 1)

        # Collect variable names in a list
        for ivar in range(num_vars):
            vcode = unpack(end_char + '%is' % (rec_len),
                           fin.read(rec_len))[0]
            var = vcode.decode('utf-8').replace(" ", "")
            data['vars'].append(var)
            dummy, rec_lec = unpack(end_char + '2l', fin.read(8))

        # Extract time
        rec_time = np.array(unpack(end_char + 'lllllll', fin.read(28)))
        rec_time[-1] *= 1000  # convert from millisec to microsec
        data["time"] = dt.datetime(*rec_time)

        # Header is this length:
        # Version + start/stop byte
        # nlons, nlats, nalts + start/stop byte
        # num_vars + start/stop byte
        # variable names + start/stop byte
        # time + start/stop byte

        iheader_length = 84 + num_vars * 48

        ntotal = data["nlons"] * data["nlats"] * data["nalts"]
        idata_length = ntotal * 8 + 8

        # Save the data for the desired variables
        for ivar in file_vars:
            fin.seek(iheader_length + ivar * idata_length)
            sdata = unpack(end_char + 'l', fin.read(4))[0]
            data[ivar] = np.array(
                unpack(end_char + '%id' % (ntotal), fin.read(sdata))).reshape(
                    (data["nlons"], data["nlats"], data["nalts"]), order="F")

    return data
