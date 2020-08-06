# pynadjust apu module - retrieve data from .apu files

import geodepy.convert as gc
import numpy as np
import geodepy.statistics as gstat
import scipy.spatial as sspat
from pynadjust_classes import Station, DynaMetadata
import common_fn
import xyz

# ----------------------------------------------------------------------
# Classes and functions
# ----------------------------------------------------------------------


class DynaApu(object):
    def __init__(self, filename, stns=None, metadata=None):

        def read_apu(apu_file, stns, read_metadata):
            """
            Function to consume DynAdjust *.apu file and read contents into a DynaApu object
            note: Variance matrix components are stored in XYZ rotation ONLY.
                      - ENU Variance matrix components are rotated back to XYZ.
                      - ENU Covariance block matrix components are rotated back to XYZ, at the position of the 'first' station.
            :param apu_file: input *.apu file
            :param stns: either None, or dictionary of station objects from earlier adj/xyz read
            :param read_metadata: True/False switch to read metadata from file, or carry across from earlier xyz/adj read
            :return: DynaApu object fields (dictionary of stations) apu options, metadata
            """

            with open(apu_file) as apu_fh:
                metadata_switch = True
                rotate_vcv = False
                pu_switch = False
                pu_line = None
                cov_count = 0

                version = None
                file_name = None
                file_date = None
                confidence_interval = None
                vcv_blocks = False
                variance_matrix_units = None
                full_covariance_matrix = None
                metadata = None

                # set xyz_stns switch true to add apu info to existing stns dictionary from earlier xyz/adj read
                if stns:
                    xyz_stns = True
                else:
                    xyz_stns = False
                    stns = {}

                for line_count, line in enumerate(apu_fh):

                    if metadata_switch:
                        version = common_fn.read_metadata(line, 'Version:', version)
                        file_name = common_fn.read_metadata(line, 'File name:', file_name)
                        file_date = common_fn.read_file_date(line, file_date)
                        confidence_interval = common_fn.read_metadata(line, 'PU confidence interval', confidence_interval)
                        variance_matrix_units = common_fn.read_metadata(line, 'Variance matrix units ', variance_matrix_units)
                        vcv_blocks = common_fn.read_metadata_tf(line, 'Stations printed in blocks', vcv_blocks)
                        full_covariance_matrix = common_fn.read_metadata_tf(line, 'Full covariance matrix', full_covariance_matrix)

                        if "Positional uncertainty of adjusted station coordinates" in line:
                            pu_line = line_count + 5
                            metadata_switch = False
                            if variance_matrix_units == 'ENU':
                                rotate_vcv = True
                            if read_metadata:
                                metadata = DynaMetadata(epoch=None, reference_frame=None, geoid_model=None, version=version)
                            apu_options = DynaApu.ApuOptions(
                                confidence_interval=confidence_interval,
                                variance_matrix_units=variance_matrix_units,
                                vcv_blocks=vcv_blocks,
                                full_covariance_matrix=full_covariance_matrix
                            )

                    if pu_line:
                        if line_count == pu_line:
                            pu_switch = True

                    if pu_switch:
                        if line == '\n':
                            continue

                        if 'block' in line.lower():
                            continue

                        if '-'*30 in line:
                            continue

                        cols = line.split()
                        num_cols = len(cols)

                        # account for station names with spaces.
                        temp = line[0:20]
                        if temp != (' ' * 20):
                            num_cols = len(line[21:].split()) + 1

                        # Station variance line 1
                        if num_cols == 11:
                            stn = line[:20].strip()
                            lat = gc.hp2dms(float(line[23:36]))
                            lon = gc.hp2dms(float(line[38:51]))
                            hpu = float(line[51:62].strip())
                            vpu = float(line[62:73].strip())
                            smaj = float(line[73:86].strip())
                            smin = float(line[86:99].strip())
                            brg = float(line[99:112].strip())
                            vcv_11 = float(line[112:131].strip())
                            vcv_12 = float(line[131:150].strip())
                            vcv_13 = float(line[150:].strip())
                            continue

                        # Station variance line 2
                        elif num_cols == 2:
                            vcv_22 = float(line[131:150].strip())
                            vcv_23 = float(line[150:].strip())
                            continue

                        #  Station variance line 3
                        elif num_cols == 1:
                            vcv_33 = float(line[150:].strip())

                            vcv = np.array([
                                [vcv_11, vcv_12, vcv_13],
                                [vcv_12, vcv_22, vcv_23],
                                [vcv_13, vcv_23, vcv_33],
                            ])

                            if rotate_vcv:
                                vcv = gstat.vcv_local2cart(vcv, lat.dec(), lon.dec())

                            # update existing or create new station object
                            if xyz_stns:
                                stns[stn].hpu = hpu
                                stns[stn].vpu = vpu
                                stns[stn].vpu = vpu
                                stns[stn].smaj = smaj
                                stns[stn].smin = smin
                                stns[stn].brg = brg
                                stns[stn].vcv = vcv
                            else:
                                apu_stn = Station(
                                    name=stn,
                                    con=None,
                                    lat=lat,
                                    lon=lon,
                                    ehgt=None,
                                    ohgt=None,
                                    sd_e=None,
                                    sd_n=None,
                                    sd_u=None,
                                    description=None,
                                    hpu=hpu,
                                    vpu=vpu,
                                    smaj=smaj,
                                    smin=smin,
                                    brg=brg,
                                    vcv=vcv,
                                    covariances={}
                                )

                                stns[stn] = apu_stn

                        # covariance block information
                        if full_covariance_matrix:
                            if num_cols == 4:
                                cov_stn = line[:20].strip()
                                cov_11 = float(line[113:131].strip())
                                cov_12 = float(line[132:150].strip())
                                cov_13 = float(line[151:169].strip())
                                cov_count = 1

                            elif num_cols == 3:
                                cov_count += 1
                                if cov_count == 2:
                                    cov_21 = float(line[113:131].strip())
                                    cov_22 = float(line[132:150].strip())
                                    cov_23 = float(line[151:169].strip())
                                elif cov_count == 3:
                                    cov_31 = float(line[113:131].strip())
                                    cov_32 = float(line[132:150].strip())
                                    cov_33 = float(line[151:169].strip())

                                    cov = np.array([
                                        [cov_11, cov_12, cov_13],
                                        [cov_21, cov_22, cov_23],
                                        [cov_31, cov_32, cov_33],
                                    ])

                                    if rotate_vcv:
                                        cov = gstat.vcv_local2cart(cov, stns[stn].lat.dec(), stns[stn].lon.dec())

                                    if cov_stn not in stns[stn].covariances:
                                        stns[stn].covariances[cov_stn] = cov

            return stns, metadata, file_name, file_date, apu_options

        read_metadata = True
        if metadata:
            read_metadata = False
            self.metadata = metadata

        result = read_apu(filename, stns, read_metadata=read_metadata)

        self.stns = result[0]
        if not metadata:
            self.metadata = result[1]
        self.file_name = result[2]
        self.file_date = result[3]
        self.options = result[4]
        self.relative_uncertainties = []
        self.ru_strategy = None
        self.ru_neighbours = None

    class ApuOptions(object):
        def __init__(self, confidence_interval, vcv_blocks, variance_matrix_units, full_covariance_matrix):
            self.confidence_interval = confidence_interval
            self.vcv_blocks = vcv_blocks
            self.variance_matrix_units = variance_matrix_units
            self.full_covariance_matrix = full_covariance_matrix

    class RelativeUncertainty(object):
        def __init__(self):
            self.stn1 = None
            self.stn2 = None
            self.ru_hz = None
            self.ru_vt = None
            self.ree_smaj = None
            self.ree_smin = None
            self.ree_brg = None
            self.covariance = None


def retrieve_covariance(apu_object, stn1, stn2):
    """
    function to find 9-element (3 x 3 arrary) covariance block in a DynaApu object
    :param apu_object: DynaApu object type
    :param stn1: Station object for first station
    :param stn2: Station object for second station
    :return: True/False,
             3 x 3 Covariance block, and
             stn covariance block is rotated according to (for covariances in ENU)
    """
    cov = np.zeros([3, 3])
    cov_found = False
    rot_stn = None

    if apu_object.options.full_covariance_matrix:
        if apu_object.stns[stn1].covariances:
            if stn2 in apu_object.stns[stn1].covariances:
                cov = apu_object.stns[stn1].covariances[stn2]
                cov_found = True
                rot_stn = stn1

        if apu_object.stns[stn2].covariances:
            if stn1 in apu_object.stns[stn2].covariances:
                cov = apu_object.stns[stn2].covariances[stn1]
                cov_found = True
                rot_stn = stn2

    return cov_found, cov, rot_stn


def compute_rel_uncertainty(apu_object, strategy, neighbours=15):
    """
    function to compute relative uncertainties between stations in a DynaApu object
    :param apu_object: DynaApu object from read_apu() call
    :param strategy: "all" computes relative uncertainties between all stations (many-to-many)
                     "nearest" computes relative uncertainties between a subject station and its nearest neighbours
    :param neighbours: Number of surrounding stations for which relative uncertainty is computed (using "nearest"
                       strategy)
    """

    def add_rel_unc(apu_object, stn1, stn2):
        """
        function to find variance matrices and covariance components, compute relative uncertainties and add to DynaApu
        object
        :param apu_object: DynaApu object from read_apu() call
        :param stn1: Station object for first station
        :param stn2: Station object for first station
        :return:
        """
        lat = apu_object.stns[stn1].lat.dec()
        lon = apu_object.stns[stn1].lon.dec()
        vcv1 = apu_object.stns[stn1].vcv
        vcv2 = apu_object.stns[stn2].vcv

        cov_used, cov, rot_stn = retrieve_covariance(apu_object, stn1, stn2)

        r_err = gstat.relative_error(lat, lon, vcv1, vcv2, cov)
        r_hu = gstat.circ_hz_pu(r_err[0], r_err[1])
        r_vu = r_err[3] * 1.96

        ru = DynaApu.RelativeUncertainty()
        ru.stn1 = stn1
        ru.stn2 = stn2
        ru.ru_hz = r_hu
        ru.ru_vt = r_vu
        ru.ru_smaj = r_err[0]
        ru.ru_smin = r_err[1]
        ru.ru_brg = r_err[2]
        ru.covariance = cov_used

        apu_object.relative_uncertainties.append(ru)

    # check valid option given for RU strategy
    strategies = (
        'all',
        'nearest'
    )

    if strategy not in strategies:
        raise ValueError(f'{strategy} strategy option not recognised')

    apu_object.ru_strategy = strategy

    stns_list = list(apu_object.stns.keys())

    if strategy == 'all':
        for stn1 in apu_object.stns:

            for stn2 in stns_list:
                if stn1 == stn2:
                    continue

                add_rel_unc(apu_object, stn1, stn2)

    elif strategy == 'nearest':
        apu_object.ru_neighbours = neighbours

        # create coordinates array
        stns_array = np.zeros((len(stns_list), 2))

        stn_count = 0
        for s in stns_list:
            stns_array[stn_count, 0] = apu_object.stns[s].lat.dec()
            stns_array[stn_count, 1] = apu_object.stns[s].lon.dec()
            stn_count += 1

        # form kd tree
        kd = sspat.KDTree(stns_array)

        if neighbours > len(stns_list):
            raise ValueError(f' *** Network too small for {neighbours} connections.\n')

        for stn1 in stns_list:

            point = np.array((apu_object.stns[stn1].lat.dec(), apu_object.stns[stn1].lon.dec()))
            result = kd.query(point, neighbours + 1)
            for i in range(1, neighbours + 1):
                stn2 = stns_list[result[1][i]]

                add_rel_unc(apu_object, stn1, stn2)


# ----------------------------------------------------------------------
# Example use
# ----------------------------------------------------------------------

apu_file = ''
xyz_file = ''

if xyz_file:
    xyz_results = xyz.DynaXYZ(xyz_file)

if apu_file:
    if xyz_file:
        apu_results = DynaApu(apu_file, xyz_results.stations, xyz_results.metadata)
    else:
        apu_results = DynaApu(apu_file)

    # Note: option strategy='all' will compute many-to-many RU. Use with caution on large networks!
    # compute_rel_uncertainty(apu_results, strategy='all')
    compute_rel_uncertainty(apu_results, strategy='nearest', neighbours=5)

    # Metadata
    out_str = 'PynAdjust apu.py example\n'
    out_str += '-' * 80 + '\n'
    out_str += f'File name:              {apu_results.file_name}\n'
    out_str += f'File date:              {apu_results.file_date}\n'
    if xyz_file:
        out_str += f'Reference frame:        {apu_results.metadata.reference_frame}\n'
        out_str += f'Epoch:                  {apu_results.metadata.epoch}\n'
        out_str += f'Geoid model:            {apu_results.metadata.geoid_model}\n'
    out_str += f'Version:                {apu_results.metadata.version}\n'
    out_str += f'Confidence interval:    {apu_results.options.confidence_interval}\n'
    out_str += f'Variance matrix units:  {apu_results.options.variance_matrix_units}\n'
    out_str += f'Full covariances:       {apu_results.options.full_covariance_matrix}\n'
    out_str += '-' * 80 + '\n'

    # Positional Uncertainty
    out_str += '\n\nPositional Uncertainties\n\n'
    out_str += '{:20s} {:>13s} {:>13s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}\n'.format(
        'Station',
        'Latitude',
        'Longitude',
        'Hz_PU',
        'Vt_PU',
        'SMaj',
        'SMin',
        'Brg'
    )
    out_str += '-' * 120 + '\n'
    for s in apu_results.stns:
        out_str += '{:20s} {:13.8f} {:13.8f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.2f}\n'.format(
            s,
            apu_results.stns[s].lat.dec(),
            apu_results.stns[s].lon.dec(),
            apu_results.stns[s].hpu,
            apu_results.stns[s].vpu,
            apu_results.stns[s].smaj,
            apu_results.stns[s].smin,
            apu_results.stns[s].brg
        )

    # Relative Uncertainty
    if apu_results.ru_strategy == 'nearest':
        out_str += f'\n\nRelative Uncertainties ({apu_results.ru_strategy} {apu_results.ru_neighbours} stations)\n\n'
    else:
        out_str += f'\n\nRelative Uncertainties ({apu_results.ru_strategy} stations)\n\n'

    out_str += '{:20s} {:20s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s} {:>11s}\n'.format(
        'Station 1',
        'Station 2',
        'Rel_Hz',
        'Rel_Vt',
        'REE_SMaj',
        'REE_SMin',
        'REE_Brg',
        'Covariance',
    )
    out_str += '-' * 120 + '\n'

    for r in apu_results.relative_uncertainties:
        out_str += '{:20s} {:20s} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:>11s}\n'.format(
            r.stn1,
            r.stn2,
            r.ru_hz,
            r.ru_vt,
            r.ru_smaj,
            r.ru_smin,
            r.ru_brg,
            str(r.covariance),
        )

    # write to file
    with open('apu_results.txt', 'w') as out_fh:
        out_fh.write(out_str)

    # write shapefiles
    common_fn.write_stns_shapefile(apu_results.stns, 'apu_results', ref_frame=apu_results.metadata.reference_frame)

    # # retrieve specific stations and write to file - caution: requested stations must be present in file
    # r_stns = retrieve_stns(apu_results, 'ARMD', 'TS7379')
    #
    # out_str = 'Station,Latitude,Longitude,hpu,vpu,smaj,smin,brg\n'
    # with open('retrieve_key_stns.txt', 'w') as out_fh:
    #     out_fh.write(out_str)
    #     for r in r_stns:
    #         print(r.name, r.lat.dec(), r.lon.dec(), r.hpu, r.vpu, r.smaj, r.smin, r.brg, sep=',', file=out_fh)
