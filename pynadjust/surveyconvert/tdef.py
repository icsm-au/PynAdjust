#!/usr/bin/env python3

"""
PynAdjust - Survey Convert - Trimble Data Exchange Format (TDEF) Module
Author - Josh Batchelor
Purpose - Enables conversion of TDEF GNSS data to DynaML format

For information about DynaML schema definition, see DynAdjust Users Guide
Appendix B.3

DynAdjust Users Guide available here:
https://github.com/icsm-au/DynAdjust/blob/master/resources/DynAdjust%20Users%20Guide.pdf
"""

import os
import shutil
import argparse
import lxml.etree as ET
from geodepy.convert import dec2hp
from pynadjust.surveyconvert.dynaml import (dnaxmlroot,
                                            addstnrecord,
                                            addgnssbaseline)


def stn2xml(tdeffile):
    # Read tdef Projection Zone
    with open(tdeffile) as myfile:
        for line in myfile.readlines():
            if line.startswith('ProjCoordinateZone'):
                zone = line[-3:-1]
    # Create dynaxml root object
    stnroot = dnaxmlroot('stn')
    # Read tdef text, list stations and separate variables by ':'
    stnlines = []
    with open(tdeffile) as myfile:
        for line in myfile.readlines():
            if line.startswith('Station'):
                stnlines.append(line)
    for num, i in enumerate(stnlines):
        stnlines[num] = i.split(':')
    # Write station list components to xml
    for stn in stnlines:
        name = stn[2]
        constraint = 'FFF'
        coordtype = 'LLH'
        # Convert Lats and Lons based on last character
        # convert to ddd.mmssss notation
        if stn[3][-1:] == 'S':
            xaxis = '-' + str(dec2hp(float(stn[3][:-1])))
            hemi = 'S'
        elif stn[3][-1:] == 'N':
            xaxis = str(dec2hp(float(stn[3][:-1])))
            hemi = 'N'
        if stn[4][-1:] == 'W':
            yaxis = '-' + str(dec2hp(float(stn[4][:-1])))
        elif stn[4][-1:] == 'E':
            yaxis = str(dec2hp(float(stn[4][:-1])))
        height = stn[5]
        desc = stn[2] + ', ?, ?'
        stnroot = addstnrecord(stnroot, name, constraint, coordtype,
                               xaxis, yaxis, height, hemi + zone, desc)

    return stnroot


def msr2xml(tdeffile):
    # Create dynaxml root object
    msrroot = dnaxmlroot('msr')
    # Read tdef text, list vectors and separate variables by ':'
    msrlines = []
    with open(tdeffile) as myfile:
        for line in myfile.readlines():
            if line.startswith('Vector'):
                msrlines.append(line)
    for num, i in enumerate(msrlines):
        msrlines[num] = i.split(':')
    # Write measurement list components to xml
    for vector in msrlines:
        refframe = 'ITRF2014'
        epoch = vector[19].replace(' ', '.')
        scale = '1.0'
        firststn = vector[2]
        secondstn = vector[3]
        x = vector[4]
        y = vector[5]
        z = vector[6]
        sxx = vector[7]
        sxy = vector[8]
        sxz = vector[9]
        syy = vector[10]
        syz = vector[11]
        szz = vector[12]
        msrroot = addgnssbaseline(msrroot, refframe, epoch, scale, scale,
                                  scale, scale, firststn, secondstn,
                                  x, y, z, sxx, sxy, sxz, syy, syz, szz)

    return msrroot


if __name__ == "__main__":
    # Build Command Line Argument Parser
    parser = argparse.ArgumentParser(description='Pynadjust TDEF to DynaML'
                                                 'Converter\nConverts Trimble'
                                                 'Data Exchange Format *.asc'
                                                 'file to DynaML *stn.xml and'
                                                 '*msr.xml files, along with'
                                                 'DynaML.xsd file')
    parser.add_argument('-i', '--input', required=True, help='source TDEF file',
                        type=os.path.abspath)
    args = parser.parse_args()
    # Get name of tdef file to convert
    tdef_file = args.input
    tdef_path = os.path.split(tdef_file)[0]
    tdef_name = os.path.splitext(tdef_file)[0]
    # Get directory of DynaML.xsd file (same as Script)
    script_path = os.path.abspath(os.path.realpath(__file__))
    script_dir, script_name = os.path.split(script_path)
    xsd_file = os.path.join(script_dir, 'DynaML.xsd')
    # Create Station and Measure Root Files
    stnroot = stn2xml(tdef_file)
    msrroot = msr2xml(tdef_file)
    # Convert Root Files to strings
    stndata = ET.tostring(stnroot, pretty_print='True',
                          xml_declaration='True', encoding='utf-8')
    msrdata = ET.tostring(msrroot, pretty_print='True',
                          xml_declaration='True', encoding='utf-8')
    # Output to new files
    with open(tdef_name + 'stn.xml', 'wb') as stnxml:
        stnxml.write(stndata)
    with open(tdef_name + 'msr.xml', 'wb') as msrxml:
        msrxml.write(msrdata)
    # Copy DynaML.xsd file to tdef file directory
    shutil.copy(xsd_file, tdef_path)
