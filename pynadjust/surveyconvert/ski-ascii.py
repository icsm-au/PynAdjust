#!/usr/bin/env python3

"""
PynAdjust - Survey Convert - Leica Infinity SKI-ASCII Format Module
Author - Ross Bugg (adopted from Nuddin Tengku's ski-ascii.py)
Purpose - To provide functionality to enable conversion of SKI-ASCII (Leica Infinity output) GNSS data to DynaML format
Notes - SKI-ASCII needs to be exported as 'WGS-84 Cartesian'

MODIFIED CODE FROM ORIGINAL SKI-ASCII.py by Nuddin Tengku
--> See Line 23 for a necessary Pynadjust module suggested location on your PC <--
Modifications start at Line 104, Line 126 and Line 146

For information about DynaML schema definition, see DynAdjust Users Guide Appendix B.3

DynAdjust Users Guide available here:
https://github.com/icsm-au/DynAdjust/blob/master/resources/DynAdjust%20Users%20Guide.pdf
"""

import os
import lxml.etree as ET
import sys

sys.path.append('C:/path-where-pynadjust-folder-is-located-suggest-/Python/site-packages')

from pynadjust.surveyconvert.dynaml import *
from geodepy.transform import xyz2llh, geo2grid
from geodepy.convert import dec2hp

# Example 'Pynadjust\\tests\\resources\\skiasciisample.asc'

def stn2xml(skiasciifile):
    # Create dynaxml root object
    stnroot = dnaxmlroot('stn')
    # Read ski-ascii text, list stations and separate variables by ':'
    stnlines = []
    with open(skiasciifile) as myfile:
        for line in myfile.readlines():
            stnlines.append(line)
    for num, i in enumerate(stnlines):
        stnlines[num] = i.split(' ')
        stnlines[num] = list(filter(None, stnlines[num]))
    # Write station list components to xml
    linebatchdict = {}
    linebatch = 0
    namelist = []
    # Station constants
    constraint = 'FFF'
    coordtype = 'LLH'
    for line in stnlines:
        linetype = line[0][0:2]
        if line[0] == "@%Coordinate":
            outputformat = line[2].rstrip()
        if linetype == '@#':
            # Start of data & station coordinates
            name = line[0][2:]
            llh = xyz2llh(float(line[1]), float(line[2]), float(line[3]))
            xaxis = str(dec2hp(float(llh[0])))
            yaxis = str(dec2hp(float(llh[1])))
            height = str(round(llh[2], 4))
            gridcoords = geo2grid(llh[0], llh[1])
            zone = gridcoords[0][0] + str(gridcoords[1])
            desc = line[0][2:] + ', ?, ?'
            linebatch += 1
            # Add station coordinates, takes only one approximate coordinate
            # Code will only take if output format is Cartesian
            if linebatch >= 1:
                if outputformat == 'Cartesian' and name not in namelist:
                    stnroot = addstnrecord(stnroot, name, constraint, coordtype, xaxis, yaxis, height, zone, desc)
                    namelist.append(line[0][2:])
        if linebatch in linebatchdict:
            linebatchdict[linebatch].append(line)
        else:
            linebatchdict[linebatch] = [line]

    return stnroot


def msr2xml(skiasciifile):
    # Create dynaxml root object
    msrroot = dnaxmlroot('msr')
    # Read ski-ascii text, list vectors and separate variables by ':'
    msrlines = []
    with open(skiasciifile) as myfile:
        for line in myfile.readlines():
            msrlines.append(line)
    for num, i in enumerate(msrlines):
        msrlines[num] = i.split(' ')
        msrlines[num] = list(filter(None, msrlines[num]))
    # Write measurement list components to xml
    linebatchdict = {}
    linebatch = 0
    # Measurement constants
    refframe = 'ITRF2014'
    vscale = '100.000'
    scale = '1.000'
    for line in msrlines:
        linetype = line[0][0:2]
        if linetype == '@#':
        
       # Start of line batch
            linebatch += 1

       # Sigma a posteriori
        if linetype == '@&' and linebatchdict[linebatch][0][4] == 'MEAS':
            sxx = str(float(line[2])/(1/(float(line[1])**2)))
            sxy = str(float(line[3])/(1/(float(line[1])**2)))
            sxz = str(float(line[4])/(1/(float(line[1])**2)))
            syy = str(float(line[5])/(1/(float(line[1])**2)))
            syz = str(float(line[6])/(1/(float(line[1])**2)))
            szz = str(float(line[7].rstrip())/(1/(float(line[1])**2)))
            
        # Individual baseline information (Reference point of baseline and its coordinates)
        if linetype == '@+':
            firststn = line[0][2:]
            
        # Baseline vector components
        if linetype == "@-":
            secondstn = line[0][2:]
            x = line[1]
            y = line[2]
            z = line[3].rstrip()
            
        # Date and time of first common epoch
        if linetype == '@*':
            epoch = line[0][2:]
    
        #Add baseline when next batch of data is cycled
            msrroot = addgnssbaseline(msrroot, refframe, epoch, vscale, scale, scale, scale, firststn, secondstn,
                                              x, y, z, sxx, sxy, sxz, syy, syz, szz)
        if linebatch in linebatchdict:
            linebatchdict[linebatch].append(line)
        else:
            linebatchdict[linebatch] = [line]
    
    return msrroot


if __name__ == "__main__":
    skiasciifile = os.sys.argv[1]
    skiasciiname = os.path.splitext(skiasciifile)[0]
    stnroot = stn2xml(skiasciifile)
    msrroot = msr2xml(skiasciifile)
    stndata = ET.tostring(stnroot, pretty_print='True', xml_declaration='True', encoding='utf-8')
    msrdata = ET.tostring(msrroot, pretty_print='True', xml_declaration='True', encoding='utf-8')

    with open(skiasciiname + '_stn.xml', 'wb') as stnxml:
        stnxml.write(stndata)
    with open(skiasciiname + '_msr.xml', 'wb') as msrxml:
        msrxml.write(msrdata)
