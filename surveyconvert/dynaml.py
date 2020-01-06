import lxml.etree as ET

"""
For information about DynaML schema definition, see DynAdjust Users Guide Appendix B.3

DynAdjust Users Guide available here:
https://github.com/icsm-au/DynAdjust/blob/master/resources/DynAdjust%20Users%20Guide.pdf
"""


def dnaxmlroot(dnaxmlformattype):
    """
    Create Empty DynaML header for either Station or Measurement File as lxml container object
    :param dnaxmlformattype: Either 'stn' or 'msr'
    :return: lxml.etree Element Object
    """
    if dnaxmlformattype == 'stn':
        formattext = 'Station File'
    elif dnaxmlformattype == 'msr':
        formattext = 'Measurement File'
    else:
        raise ValueError("ValueError: dnaxmlformattype must be either 'stn' or 'msr'")
    NS = 'http://www.w3.org/2001/XMLSchema-instance'
    location_attribute = '{%s}noNameSpaceSchemaLocation' % NS
    dnaxmlroot = ET.Element('DnaXmlFormat', attrib={location_attribute: 'DynaML.xsd'})
    dnaxmlroot.set('type', formattext)
    return dnaxmlroot


def addstnrecord(dnaxmlroot, name, constraint, coordtype, xaxis, yaxis, height, hemi_zone, desc):
    """
    Add DynaML Station Record to existing lxml.etree Element Object (using dnaxmlstnroot above)
    :param dnaxmlroot: lxml.etree Element Object
    :param name: Station Name, text
    :param constraint: Station Constraint ('FFF' or 'CCC'), text
    :param coordtype: Coordinate Type, text
    :param xaxis: X Axis Coordinate, text
    :param yaxis: Y Axis Coordinate, text
    :param height: Height, text
    :param hemi_zone: Hemisphere and Zone, text
    :param desc: Station Description, text
    :return: lxml.etree Element Object with new Station Record added
    """
    # Define DnaStation
    dnastn = ET.SubElement(dnaxmlroot, 'DnaStation')

    # Create DnaStation Subelements name, constraints and type
    stnname = ET.SubElement(dnastn, 'Name')
    stnconst = ET.SubElement(dnastn, 'Constraints')
    stntype = ET.SubElement(dnastn, 'Type')

    # Populate Subelements
    stnname.text = name
    stnconst.text = constraint
    stntype.text = coordtype

    # Create StationCoord subelement and subsubelements
    stncoord = ET.SubElement(dnastn, 'StationCoord')
    coordname = ET.SubElement(stncoord, 'Name')
    coordx = ET.SubElement(stncoord, 'XAxis')
    coordy = ET.SubElement(stncoord, 'YAxis')
    coordh = ET.SubElement(stncoord, 'Height')
    coordhem = ET.SubElement(stncoord, 'HemisphereZone')

    # Populate StationCoord subelements
    coordname.text = name
    coordx.text = xaxis
    coordy.text = yaxis
    coordh.text = height
    coordhem.text = hemi_zone

    # Create New DnaStation Subelement, populate
    stndesc = ET.SubElement(dnastn, 'Description')
    stndesc.text = desc

    return dnaxmlroot


def addgnssbaseline(dnaxmlroot, refframe, epoch, scale, firststn, secondstn,
                    x, y, z, sxx, sxy, sxz, syy, syz, szz):
    # Define DnaMeasurement
    dnamsr = ET.SubElement(dnaxmlroot, 'DnaMeasurement')

    # Create DnaMeasurement Subelements
    msrtype = ET.SubElement(dnamsr, 'Type')
    msrrefframe = ET.SubElement(dnamsr, 'ReferenceFrame')
    msrepoch = ET.SubElement(dnamsr, 'Epoch')
    msrvscale = ET.SubElement(dnamsr, 'Vscale')
    msrpscale = ET.SubElement(dnamsr, 'Pscale')
    msrlscale = ET.SubElement(dnamsr, 'Lscale')
    msrhscale = ET.SubElement(dnamsr, 'Hscale')
    msrfirststn = ET.SubElement(dnamsr, 'First')
    msrsecondstn = ET.SubElement(dnamsr, 'Second')

    # Populate DnaMeasurement Subelements
    msrtype.text = 'G'
    msrrefframe.text = refframe
    msrepoch.text = epoch
    msrvscale.text = scale
    msrpscale.text = scale
    msrlscale.text = scale
    msrhscale.text = scale
    msrfirststn.text = firststn
    msrsecondstn.text = secondstn

    # Create DnaMeasurement GPSBaseline and Subelements
    gpsbaseline = ET.SubElement(dnamsr, 'GPSBaseline')
    msrx = ET.SubElement(gpsbaseline, 'X')
    msry = ET.SubElement(gpsbaseline, 'Y')
    msrz = ET.SubElement(gpsbaseline, 'Z')
    msrsxx = ET.SubElement(gpsbaseline, 'SigmaXX')
    msrsxy = ET.SubElement(gpsbaseline, 'SigmaXY')
    msrsxz = ET.SubElement(gpsbaseline, 'SigmaXZ')
    msrsyy = ET.SubElement(gpsbaseline, 'SigmaYY')
    msrsyz = ET.SubElement(gpsbaseline, 'SigmaYZ')
    msrszz = ET.SubElement(gpsbaseline, 'SigmaZZ')

    # Populate GPSBaseline Subelements
    msrx.text = x
    msry.text = y
    msrz.text = z
    msrsxx.text = sxx
    msrsxy.text = sxy
    msrsxz.text = sxz
    msrsyy.text = syy
    msrsyz.text = syz
    msrszz.text = szz

    return dnaxmlroot


# Test Station File Creation
# from surveyconvert.dynaml import *
# root = dnaxmlroot('stn')
# root = addstnrecord(root,
#                     '317701630', 'FFF', 'LLH',
#                     '-36.24390312756', '145.14382243848', '121.7653', 'S55',
#                     '317701630, ?, ?')
# mydata = ET.tostring(root, pretty_print='True', xml_declaration='True', encoding='utf-8')
# print(mydata)

# Write byte string to file
# with open('testdnastn.xml', 'wb') as myfile:
    # myfile.write(mydata)
