import lxml.etree as ET

"""
For information about DynaML schema definition, see DynAdjust Users Guide Appendix B.3

DynAdjust Users Guide available here:
https://github.com/icsm-au/DynAdjust/blob/master/resources/DynAdjust%20Users%20Guide.pdf
"""


def dnaxmlstnroot():
    """
    Create Empty DynaML stn header and lxml container object
    :return: lxml.etree Element Object
    """
    # Define Root Object, Add Components
    NS = 'http://www.w3.org/2001/XMLSchema-instance'
    location_attribute = '{%s}noNameSpaceSchemaLocation' % NS
    dnaxmlroot = ET.Element('DnaXmlFormat', attrib={location_attribute: 'DynaML.xsd'})
    dnaxmlroot.set('type', 'Station File')
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

# Output XML data to byte string
# mydata = ET.tostring(dnaxmlroot, pretty_print='True', xml_declaration='True', encoding='utf-8')

# Write byte string to file
# with open('testdnastn.xml', 'wb') as myfile:
    # myfile.write(mydata)
