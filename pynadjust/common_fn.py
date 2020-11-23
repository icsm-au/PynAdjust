# common functions used across xyz, adj and apu pynadjust modules

import datetime
import shapefile
import os


def write_prj(file_name, ref_frame=None):
    """
    Function to write shape projection files. GDA94 and GDA2020 reference frames are supported. Nothing is written
    for other reference frames.
    :param file_name: name of prj file to be written
    :param ref_frame: reference frame
    :return:
    """

    if ref_frame == 'GDA2020':
        out_str = 'GEOGCS["GDA2020",DATUM["GDA2020",SPHEROID["GRS_1980",6378137.0,298.257222101]],' \
                  'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433],AUTHORITY["EPSG",7844]]'
    elif ref_frame == 'GDA94':
        out_str = 'GEOGCS["GDA94",DATUM["D_GDA_1994",SPHEROID["GRS_1980",6378137,298.257222101]],' \
                  'PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]'
    else:
        out_str = None

    if out_str:
        prj_fh = open(file_name + '.prj', 'w')
        prj_fh.write(out_str)
        prj_fh.close()


def read_metadata(line, search_text, current_value):
    """
    Function to read simple header items in DynAdjust files
    :param line: DynAdjust header line
    :param search_text: header field required.
    :param current_value: stored value. Updated when search_text is successfully found.
    :return: either current value or string corresponding to search_text
    """
    if line[:35] == search_text.ljust(35, ' '):
        return line[35:].strip()
    else:
        return current_value


def read_metadata_tf(line, search_text, current_value):
    """
    function to read simple DynAdjust header items and return True/False
    :param line: DynAdjust header line
    :param search_text: header field desired.
    :param current_value: stored value. Updated when search_text is successfully found.
    :return: either current value or True/False corresponding to search_text
    """
    if line[:35] == search_text.ljust(35, ' '):
        if line[35:].strip() == 'Yes':
            return True
        else:
            return False
    else:
        return current_value


def read_file_date(line, current_value):
    """
    Function to read file created date from DynAdjust header line and return datetime object
    :param line: DynAdjust header line
    :param current_value: stored value. Updated when search_text is successfully found.
    :return: either current value or datetime object corresponding to search_text
    """
    if line[:35] == 'File created:                      ':
        date_str = line[35:].strip()
        file_date = datetime.datetime.strptime(date_str, '%A, %d %B %Y, %I:%M:%S %p')
        return file_date
    else:
        return current_value


def read_epoch(line, current_value):
    """
    function to read epoch header line and return datetime.date object
    :param line: DynAdjust header line
    :param current_value: stored value. Updated when search_text is successfully found.
    :return: either current value or datetime.date object corresponding to search_text
    """
    if line[:35] == 'Epoch:                             ':
        date_str = line[35:].strip()
        d = int(date_str[:2])
        m = int(date_str[3:5])
        y = int(date_str[6:])
        epoch = datetime.date(y, m, d)
        return epoch
    else:
        return current_value


def retrieve_stns(stns_dict, *stn):
    """
    function to retrieve Station objects for requested stations
    :param stns_dict: dictionary of Station objects from adj/apu/xyz read
    :param stn: name of station to be retrieved
    :return: list of Station objects
    """

    stns = []

    for s in stn:
        try:
            stns.append(stns_dict[s])
        except KeyError:
            raise ValueError(f'{s} not found in .apu file')

    return stns


def write_stns_shapefile(stns, network_name, ref_frame=None):
    """
    function to write shapefiles from DynAdjust adj/apu/xyz results
    :param stns: dictionary of Station objects
    :param network_name: name of network (used in shapefile name output)
    :param ref_frame: optional, reference frame of coordinates
    """
    shp_dir = False
    cwd = os.getcwd()
    if cwd.split('\\')[-1] != 'shp':
        check_enter_dir('shp')
    else:
        shp_dir = True

    shp_name = network_name + '_stn'

    write_prj(shp_name, ref_frame)

    w = shapefile.Writer(shp_name, shapeType=1)

    w.autoBalance = 1

    w.field('Station', 'C', size=20)
    w.field('Description', 'C', size=50)
    w.field('Constraints', 'C', size=3)
    w.field('Easting', 'N', decimal=4)
    w.field('Northing', 'N', decimal=4)
    w.field('Zone', 'N')
    w.field('Latitude', 'N', decimal=10)
    w.field('Longitude', 'N', decimal=10)
    w.field('OHGT', 'N', decimal=4)
    w.field('EHGT', 'N', decimal=4)
    w.field('SD_E', 'N', decimal=4)
    w.field('SD_N', 'N', decimal=4)
    w.field('SD_U', 'N', decimal=4)
    w.field('HPU', 'N', decimal=4)
    w.field('VPU', 'N', decimal=4)
    w.field('SMaj', 'N', decimal=4)
    w.field('SMin', 'N', decimal=4)
    w.field('Brg', 'N', decimal=2)

    for s in stns.values():
        w.point(s.lon.dec(), s.lat.dec())

        grid = s.grid()

        w.record(
            s.name, s.description, s.con,
            grid[2], grid[3], grid[1],
            s.lat.dec(), s.lon.dec(),
            s.ohgt, s.ehgt,
            s.sd_e, s.sd_n, s.sd_u,
            s.hpu, s.vpu,
            s.smaj, s.smin, s.brg
        )

    w.close()

    if not shp_dir:
        os.chdir('..')


def check_enter_dir(dir_name):
    """
    Function to check for the existence of a
    :param dir_name:
    :return:
    """
    try:
        os.chdir(dir_name)
    except FileNotFoundError:
        os.mkdir(dir_name)
        os.chdir(dir_name)
