# pynadjust xyz module - retrieve data from DynAdjust .xyz files

import geodepy.convert as gc
from pynadjust.pynadjust_classes import Station, DynaMetadata, Switches
import pynadjust.common_fn as common_fn


class DynaXYZ(object):
    def __init__(self, filename):
        """
        DynaXYZ object which reads all metadata and coordinate information upon initialisation
        :param filename: name of xyz file to be consumed.
        """

        def read_xyz(xyz_file):
            """
            function/method to read metadata and coordinate information from a DynAdjust station file.
            :param xyz_file: DynAdjust *.xyz file
            :return: dictionary of stn objects, and metadata variables
            """
            stns = {}
            version = None
            reference_frame = None
            file_name = None
            file_date = None
            epoch = None
            geoid_model = None

            with open(xyz_file, 'r') as xyz_fh:
                switches = Switches()
                switches.header = True
                line_count = 0
                stn_line = False
                # stn_switch = False
                # metadata_switch = True
                mandatory_coord_types = 'PLHh'
                desc_index = None

                for line_count, line in enumerate(xyz_fh):
                    # line_count += 1

                    if 'Adjusted Coordinates' in line:
                        stn_line = line_count + 5

                    if stn_line:
                        if line_count == stn_line - 2:
                            desc_index = line.find('Description')

                        if line_count == stn_line:
                            switches.reset()
                            switches.stns = True
                            # stn_switch = True
                            # metadata_switch = False

                            metadata = DynaMetadata(reference_frame=reference_frame, epoch=epoch,
                                                    geoid_model=geoid_model, version=version)

                    if switches.header:
                    # if metadata_switch:
                        version = common_fn.read_metadata(line, 'Version:', version)
                        reference_frame = common_fn.read_metadata(line, 'Reference frame:', reference_frame)
                        file_name = common_fn.read_metadata(line, 'File name:', file_name)
                        file_date = common_fn.read_file_date(line, file_date)
                        epoch = common_fn.read_epoch(line, epoch)
                        geoid_model = common_fn.read_metadata(line, 'Geoid model:', geoid_model)

                        if line[:35] == 'Station coordinate types:          ':
                            coord_types = line[35:].strip()
                            missing_coord_types = set()
                            for l in mandatory_coord_types:
                                if l not in coord_types:
                                    missing_coord_types.add(l)
                            if missing_coord_types:
                                raise ValueError(f'Mandatory coordinate types {missing_coord_types} not present in {xyz_file}')

                    # if stn_switch:
                    if switches.stns:
                        if len(line) < 20:
                            switches.reset()
                            # stn_switch = False
                            continue

                        stn_object = read_coord_elements(line, coord_types, desc_index)

                        stns[stn_object.name] = stn_object

            return stns, file_name, file_date, metadata

        results = read_xyz(filename)
        self.stns = results[0]
        self.file_name = results[1]
        self.file_date = results[2]
        self.metadata = results[3]


def read_coord_elements(line, coord_types, desc_index):
    """
    Function to read coordinate components from a line in an adj/xyz file
    :param line: coordinate line
    :param coord_types: string of coordinate types on line
    :param desc_index: position of station descriptions on line
    :return: Station object
    """

    lat = None
    lon = None
    ohgt = None
    ehgt = None

    stn = line[0:20].strip()
    con = line[20:23]
    results = line[25:].split()
    r_count = 0

    for ct in coord_types:
        if ct == 'P':
            lat = gc.hp2dms(float(results[r_count]))
        if ct == 'L':
            lon = gc.hp2dms(float(results[r_count]))
        if ct == 'H':
            ohgt = float(results[r_count])
        if ct == 'h':
            ehgt = float(results[r_count])

        r_count += 1

    # Don't forget about the qualities
    sd_e = float(results[r_count])
    sd_n = float(results[r_count + 1])
    sd_u = float(results[r_count + 2])

    # Station Description
    if desc_index is None or desc_index == -1:
        desc = ''
    else:
        desc = str(line[desc_index:].strip())

    # store details as a Station object, leaving uncertainty fields null
    stn_object = Station(
        name=stn,
        con=con,
        lat=lat,
        lon=lon,
        ehgt=ehgt,
        ohgt=ohgt,
        sd_e=sd_e,
        sd_n=sd_n,
        sd_u=sd_u,
        description=desc,
        hpu=None,
        vpu=None,
        smaj=None,
        smin=None,
        brg=None,
        vcv=None,
        covariances={}
    )

    return stn_object


# ----------------------------------------------------------------------
# Example of usage
# ----------------------------------------------------------------------

xyz_file = ''

if xyz_file:
    xyz_results = DynaXYZ(xyz_file)

    # iterate through and print coordinates
    out_str = ''
    # print metadata
    out_str += f'File name:       {xyz_results.file_name}\n'
    out_str += f'Date created:    {xyz_results.file_date}\n'
    out_str += f'Version:         {xyz_results.metadata.version}\n'
    out_str += f'Reference Frame: {xyz_results.metadata.reference_frame}\n'
    out_str += f'Epoch:           {xyz_results.metadata.epoch}\n'
    out_str += f'Geoid model:     {xyz_results.metadata.geoid_model}\n'

    # print results
    out_str += '\nStations:\n'
    out_str += '{:20s} {:>3s} {:>14s} {:>14s} {:>10s} {:>10s} {:>8s} {:>8s} {:>8s} ' \
               '{:>14s} {:>14s} {:>5s} {:>14s} {:>14s} {:>14s}\n'.format(
        'Station',
        'F/C',
        'Latitude DD',
        'Longitude DD',
        'Ell_Ht',
        'Phys_Ht',
        'SD_E',
        'SD_N',
        'SD_U',
        'Easting',
        'Northing',
        'Zone',
        'X',
        'Y',
        'Z',
    )
    out_str += '-' * 184 + '\n'

    for stn in xyz_results.stns:
        out_str += '{:20s} {:3s} {:14.8f} {:14.8f} {:10.4f} {:10.4f} {:8.4f} {:8.4f} {:8.4f} '.format(
            xyz_results.stns[stn].name,
            xyz_results.stns[stn].con,
            xyz_results.stns[stn].lat.dec(),
            xyz_results.stns[stn].lon.dec(),
            xyz_results.stns[stn].ehgt,
            xyz_results.stns[stn].ohgt,
            xyz_results.stns[stn].sd_e,
            xyz_results.stns[stn].sd_n,
            xyz_results.stns[stn].sd_u
        )

        grid = xyz_results.stns[stn].grid()
        xyz = xyz_results.stns[stn].xyz()

        out_str += '{:14.4f} {:14.4f} {:5d} {:14.4f} {:14.4f} {:14.4f}\n'.format(
            grid[2],
            grid[3],
            grid[1],
            xyz[0],
            xyz[1],
            xyz[2]
        )

    with open('pynadjust_xyz.txt', 'w') as fh:
        fh.write(out_str)
