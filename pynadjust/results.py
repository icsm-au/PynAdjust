# pynadjust results module - retrieve data from DynAdjust output files and perform various functions

import datetime
import geodepy.convert as gc
import geodepy.statistics as gstat
import shapefile
import numpy as np
import os
import re
import scipy.spatial as sspat
import warnings

# ----------------------------------------------------------------------
# sets for various msr types
# ----------------------------------------------------------------------

msr_types = ('A', 'B', 'C', 'D', 'E', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'P', 'Q', 'R', 'S', 'V', 'X', 'Y', 'Z')
three_line_msrs = ('G', 'X', 'Y')
multi_line_msrs = ('D', ' ', )
one_stn_msrs = ('H', 'I', 'J', 'P', 'Q', 'R', 'Y')
two_stn_msrs = ('B', 'C', 'E', 'G', 'K', 'L', 'M', 'S', 'V', 'X', 'Z')
three_stn_msrs = ('A', 'D')
angle_msrs = ('A', 'B', 'D', 'V', 'Z', 'P', 'Q')  # angle msr types usually expressed in "DDD MM SS.SSSS" format (not DEC or HP)
hp_msrs = ('I', 'J')


# ----------------------------------------------------------------------
# Classes
# ----------------------------------------------------------------------

class AdjMetadata(object):
    def __init__(self, epoch=None, reference_frame=None, geoid_model=None, version=None, build=None,
                 variance_matrix_units=None, full_covariance_matrix=None, vcv_blocks=None):
        self.epoch = epoch
        self.reference_frame = reference_frame
        self.geoid_model = geoid_model
        self.version = version
        self.build = build
        self.variance_matrix_units = variance_matrix_units
        self.full_covariance_matrix = full_covariance_matrix
        self.vcv_blocks = vcv_blocks


class FileMetadata(object):
    def __init__(self, adj_filename=None, apu_filename=None, xyz_filename=None, adj_date_created=None,
                 apu_date_created=None, xyz_date_created=None):
        self.adj_filename = adj_filename
        self.apu_filename = apu_filename
        self.xyz_filename = xyz_filename
        self.adj_date_created = adj_date_created
        self.apu_date_created = apu_date_created
        self.xyz_date_created = xyz_date_created


class AdjStats(object):
    def __init__(self, soln_type=None, run_time=None, parameters=None, msr_count=None, degrees_of_freedom=None,
                 sigma_zero=None, chi_squared=None, outlier_count=None, global_pelzer=None, chi_square_lower=None,
                 chi_square_upper=None, chi_square_result=None, confidence_interval=None):
        self.soln_type = soln_type
        self.run_time = run_time
        self.parameters = parameters
        self.msr_count = msr_count
        self.degrees_of_freedom = degrees_of_freedom
        self.sigma_zero = sigma_zero
        self.chi_squared = chi_squared
        self.outlier_count = outlier_count
        self.global_pelzer = global_pelzer
        self.chi_square_lower = chi_square_lower
        self.chi_square_upper = chi_square_upper
        self.chi_square_result = chi_square_result
        self.confidence_interval = confidence_interval


class GroupVF(object):
    def __init__(self, group=None, count=None, vs2=None, gvf=None):
        self.group = group
        self.count = count
        self.vs2 = vs2
        self.gvf = gvf


class RelativeUncertainty(object):
    def __init__(self, stn1=None, stn2=None, ru_hz=None, ru_vt=None, ree_smaj=None, ree_smin=None, ree_brg=None,
                 covariance=None):
        self.stn1 = stn1
        self.stn2 = stn2
        self.ru_hz = ru_hz
        self.ru_vt = ru_vt
        self.ree_smaj = ree_smaj
        self.ree_smin = ree_smin
        self.ree_brg = ree_brg
        self.covariance = covariance


class DynaResults(object):
    def __init__(self, adj_file=None, xyz_file=None, apu_file=None, stns={}, msrs=[], adj_metadata=AdjMetadata(),
                 file_metadata=FileMetadata(), adj_stats=AdjStats(), relative_uncertainties=[]):
        self.adj_file = adj_file
        self.apu_file = apu_file
        self.xyz_file = xyz_file
        self.stns = stns
        self.msrs = msrs
        self.adj_metadata = adj_metadata
        self.file_metadata = file_metadata
        self.adj_stats = adj_stats
        self.relative_uncertainties = relative_uncertainties

        # consume files at initialisation
        self.read_results()

    def read_results(self):
        """
        Method to consume DynAdjust output files [*.adj, *.apu and *.xyz]
        :return:
        """
        if self.xyz_file:
            read_xyz_file(
                xyz_file=self.xyz_file,
                stns=self.stns,
                adj_metadata=self.adj_metadata,
                file_metadata=self.file_metadata
            )

        if self.adj_file:
            # clear msrs before (re)loading to avoid doubling up.
            self.msrs = []

            read_adj_file(
                adj_file=self.adj_file,
                stns=self.stns,
                msrs=self.msrs,
                adj_metadata=self.adj_metadata,
                file_metadata=self.file_metadata,
                adj_stats=self.adj_stats
            )

        if self.apu_file:
            read_apu_file(
                apu_file=self.apu_file,
                stns=self.stns,
                adj_metadata=self.adj_metadata,
                file_metadata=self.file_metadata,
                adj_stats=self.adj_stats
            )

    def link_source_with_msr_id(self, xml_msr_file):
        """
        method to link the database id tags between msr.xml file and *.adj file to assign source and epoch to adjusted
        measurements
        :param xml_msr_file: input measurement xml file
        :return:
        """

        metadata = {}

        # find source/epoch info in msr file
        with open(xml_msr_file, 'r') as f:
            for line in f:
                if '<DnaMeasurement>' in line:
                    epoch = None
                    source = None
                    msr_id = None

                if '<Source>' in line:
                    source = re.findall("<Source>(.*?)</Source>", line)[0]
                if '<Epoch>' in line:
                    epoch = re.findall("<Epoch>(.*?)</Epoch>", line)[0]
                if '<MeasurementID>' in line:
                    msr_id = re.findall("<MeasurementID>(.*?)</MeasurementID>", line)[0]

                    if msr_id in metadata:
                        raise ValueError(f' *** repeat occurance of msr_id: {msr_id}')
                    if msr_id != '':
                        metadata[msr_id] = {
                            'epoch': epoch,
                            'source': source,
                        }

        # match with *.adj results using msr id
        for m in self.msrs:
            if m.msr_type != 'D':
                if m.msr_id in metadata:
                    m.epoch = metadata[m.msr_id]['epoch']
                    m.source = metadata[m.msr_id]['source']
            else:
                if m.msr_id:
                    if m.msr_id[0] in metadata:
                        m.epoch = metadata[m.msr_id[0]]['epoch']
                        m.source = metadata[m.msr_id[0]]['source']

    def write_shp(self, network_name):
        """
        Method to write shapefiles of network.
        :param network_name: name of network - must not contain full stops (truncated by pyshp).
        :return:
        """

        if self.stns:
            check_enter_dir('shp')
            write_stn_shapefile(self.stns, network_name, ref_frame=self.adj_metadata.reference_frame)

            if self.msrs:
                write_msr_shapefile(stns=self.stns, msrs=self.msrs, network_name=network_name,
                                    ref_frame=self.adj_metadata.reference_frame)

            os.chdir('..')

    def compute_group_vf(self, source=None):
        """
        Method to compute estimated variances factors for measurement groups/types.
        Note: this computation assumes measurements are uncorrelated. use results with caution for correlated
        measurement types e.g. X, D, etc.
        :param source: Optional => Computation only performed on measurements in a specified source
        :return: list of GroupVF objects (sorted by group), overall estimate of variance factor
        """

        def increment_vs2_sum(group_vf_dict, msr_type, vs2):
            """
            Function to increment the vs2 sum for group vf count
            :param group_vf_dict: dictionary of GroupVF objects
            :param msr_type: String. VF group/msr type
            :param vs2: Float. (correction to msr divided by msr sd)squared
            :return:
            """

            try:
                group_vf_dict[msr_type].count += 1
                group_vf_dict[msr_type].vs2 += vs2
            except KeyError:
                group_vf_dict[msr_type] = GroupVF(group=msr_type, count=1, vs2=vs2)

        if not self.msrs:
            raise ValueError(f'GroupVF cannot be computed; No measurments present')

        groups = {}

        for m in self.msrs:
            if source:
                if m.source != source:
                    continue

            if m.msr_type in three_line_msrs:
                sum_vs2 = 0
                for i in range(0, 3):
                    msr_type = m.msr_type + m.cardinal[i]

                    vs2 = (m.cor[i] / m.msr_sd[i]) ** 2
                    sum_vs2 += vs2

                    # add individual G/X/Y components
                    increment_vs2_sum(groups, msr_type, vs2)

                # add total G/X/Y
                increment_vs2_sum(groups, m.msr_type, sum_vs2)

            elif m.msr_type in multi_line_msrs:
                for i in range(len(m.msr)):
                    vs2 = (m.cor[i] / m.msr_sd[i]) ** 2
                    increment_vs2_sum(groups, m.msr_type, vs2)

            else:
                vs2 = (m.cor / m.msr_sd) ** 2

                increment_vs2_sum(groups, m.msr_type, vs2)

        # sum the overall vs2 and msr counts
        total_vs2 = 0
        group_msr_count = 0
        dof = float(self.adj_stats.degrees_of_freedom)
        total_msrs = int(self.adj_stats.msr_count)
        for msr_type in groups:
            # exclude whole G/X/Y msrs. Their components are included in this summation.
            if msr_type not in three_line_msrs:
                total_vs2 += groups[msr_type].vs2
                group_msr_count += groups[msr_type].count
            else:
                groups[msr_type].count *= 3

            groups[msr_type].gvf = (total_msrs * groups[msr_type].vs2) / (groups[msr_type].count * dof)

        # compute overall estimated variance factor
        est_vf = (total_msrs / group_msr_count) * (total_vs2 / dof)

        # convert dict to list, sort by group
        group_list = list(groups.values())
        group_list.sort(key=lambda x: x.group)

        return group_list, est_vf

    def compute_rel_uncertainty(self, ru_pairs=[], compute_all=False, neighbours=0):
        """
        Method to compute relative uncertainties between station pairs. Computed RU pairs are determined by keyword args:
        :param ru_pairs: Optional, list of lists, [[stn1, stn2], [stn1, stn3], etc, ]. Computes RU for specifically requested RU pairs
        :param compute_all: Optional, True/False. Computes all possible RU pairs
        :param neighbours: Optional. Number of surrounding stations for which relative uncertainty is computed (using
                           "nearest" strategy)
        :return: Appends a RelativeUncertainty object for each RU pair to self.relative_uncertainties
        """

        # check covariances present
        if not self.adj_metadata.full_covariance_matrix:
            warnings.warn(f'Covariances not present; all RU is uncorrelated')

        # create a list of stn pairs to compute RU
        ru_list = []
        stns_list = sorted(list(self.stns.keys()))

        if ru_pairs:
            for stn1, stn2 in ru_pairs:
                if stn1 != stn2:
                    ru_list.append([stn1, stn2])

        if compute_all:
            for stn1 in stns_list:
                for stn2 in stns_list:
                    if stn1 == stn2:
                        continue
                    ru_list.append([stn1, stn2])

        if neighbours > 0:
            # create coordinates array
            stns_array = np.zeros((len(stns_list), 2))

            stn_count = 0
            for s in stns_list:
                stns_array[stn_count, 0] = self.stns[s].lat.dec()
                stns_array[stn_count, 1] = self.stns[s].lon.dec()
                stn_count += 1

            # form kd tree
            kd = sspat.KDTree(stns_array)

            if neighbours > len(stns_list):
                raise ValueError(f' *** Network too small for {neighbours} connections.\n')

            for stn1 in stns_list:
                point = np.array((self.stns[stn1].lat.dec(), self.stns[stn1].lon.dec()))
                result = kd.query(point, neighbours + 1)
                for i in range(1, neighbours + 1):
                    stn2 = stns_list[result[1][i]]
                    ru_list.append([stn1, stn2])

        # compute relative uncertainties
        for stn1, stn2 in ru_list:
            vcv1 = self.stns[stn1].vcv
            vcv2 = self.stns[stn2].vcv

            # find covariance block and determine coords for xyz/enu rotation
            # empty covariance block is used where covariances cannot be found
            cov = np.zeros([3, 3])
            cov_found = False
            rot_stn = stn1

            if self.adj_metadata.full_covariance_matrix:
                if stn2 in self.stns[stn1].covariances:
                    cov = self.stns[stn1].covariances[stn2]
                    cov_found = True
                elif stn1 in self.stns[stn2].covariances:
                    cov = np.transpose(self.stns[stn2].covariances[stn1])
                    cov_found = True
                    rot_stn = stn2

            # RU computation
            r_err = gstat.relative_error(
                lat=self.stns[rot_stn].lat.dec(),
                lon=self.stns[rot_stn].lon.dec(),
                var1=self.stns[stn1].vcv,
                var2=self.stns[stn2].vcv,
                cov12=cov
            )
            r_hz = gstat.circ_hz_pu(r_err[0], r_err[1])
            r_vt = r_err[3] * 1.96

            # store results
            self.relative_uncertainties.append(
                RelativeUncertainty(
                    stn1=stn1,
                    stn2=stn2,
                    ru_hz=r_hz,
                    ru_vt=r_vt,
                    ree_smaj=r_err[0],
                    ree_smin=r_err[1],
                    ree_brg=r_err[2],
                    covariance=cov_found
                )
            )

    def add_type_b_sd(self, stns=[], sd_e=0.006, sd_n=0.006, sd_u=0.012):
        """
        Method to apply one-sigma type b uncertainties to the station variance matrices and re-compute uncertainties.
        Can be applied to specific stations or all stations.
        :param stns: list of stations to apply uncertainties at. Defaults to all stations.
        :param sd_e: float. Type b uncertainty for east (one-sigma) in metres
        :param sd_n: float. Type b uncertainty for north (one-sigma) in metres
        :param sd_u: float. Type b uncertainty for up (one-sigma) in metres
        :return:
        """
        if not stns:
            stns = list(self.stns.keys())

        for stn in stns:
            # rotate xyz to enu
            vcv_enu = gstat.vcv_cart2local(self.stns[stn].vcv, self.stns[stn].lat.dec(), self.stns[stn].lon.dec())

            # add type b variances to enu matrix
            vcv_enu[0, 0] += sd_e**2
            vcv_enu[1, 1] += sd_n**2
            vcv_enu[2, 2] += sd_u**2

            # re-calc uncertainties
            smaj, smin, brg = gstat.error_ellipse(vcv_enu)
            hpu = gstat.circ_hz_pu(smaj, smin)

            # update attributes with new uncertainties
            self.stns[stn].vcv = gstat.vcv_local2cart(vcv_enu, self.stns[stn].lat.dec(), self.stns[stn].lon.dec())
            self.stns[stn].sd_e = vcv_enu[0, 0]**0.5
            self.stns[stn].sd_n = vcv_enu[1, 1]**0.5
            self.stns[stn].sd_u = vcv_enu[2, 2]**0.5
            self.stns[stn].hpu = hpu
            self.stns[stn].vpu = self.stns[stn].sd_u * 1.96
            self.stns[stn].smaj = smaj
            self.stns[stn].smin = smin
            self.stns[stn].brg = brg


class Measurement(object):
    def __init__(self, msr_type=None, stn1=None, stn2=None, stn3=None, ignore=None, cardinal=None, msr=None, adj=None,
                 cor=None, msr_sd=None, adj_sd=None, cor_sd=None, nstat=None, tstat=None, pelzer=None, pre_adj_cor=None,
                 outlier=None, msr_id=None, cluster_id=None, epoch=None, source=None):
        self.msr_type = msr_type
        self.stn1 = stn1
        self.stn2 = stn2
        self.stn3 = stn3
        self.ignore = ignore
        self.cardinal = cardinal
        self.msr = msr
        self.adj = adj
        self.cor = cor
        self.msr_sd = msr_sd
        self.adj_sd = adj_sd
        self.cor_sd = cor_sd
        self.nstat = nstat
        self.tstat = tstat
        self.pelzer = pelzer
        self.pre_adj_cor = pre_adj_cor
        self.outlier = outlier
        self.msr_id = msr_id
        self.cluster_id = cluster_id
        self.epoch = epoch
        self.source = source

    def __repr__(self):
        out_str = f'msr_type {self.msr_type}; stn1 {self.stn1}; stn2 {self.stn2}; stn3 {self.stn3}; ' \
                  f'ignore {self.ignore}; cardinal {self.cardinal}; msr {self.msr}; adj {self.adj}; cor {self.cor}; ' \
                  f'msr_sd {self.msr_sd}; adj_sd {self.adj_sd}; cor_sd {self.cor_sd}; nstat {self.nstat}; ' \
                  f'self.tstat {self.tstat}; pelzer {self.pelzer}; pre_adj_cor {self.pre_adj_cor}; ' \
                  f'outlier {self.outlier}; msr_id; {self.msr_id}; cluster_id; {self.cluster_id}; ' \
                  f'epoch {self.epoch}; source {self.source}'

        return out_str


class Station(object):
    # consider expanding/improving in future to utilise GeodePy Coord Classes.
    def __init__(self, name=None, description=None, con=None, lat=None, lon=None, ohgt=None, ehgt=None,
                 sd_e=None, sd_n=None, sd_u=None, hpu=None, vpu=None, smaj=None, smin=None, brg=None, vcv=None,
                 covariances=None):
        self.name = name
        self.description = description
        self.con = con
        self.lat = lat
        self.lon = lon
        self.ohgt = ohgt
        self.ehgt = ehgt
        self.sd_e = sd_e
        self.sd_n = sd_n
        self.sd_u = sd_u
        self.hpu = hpu
        self.vpu = vpu
        self.smaj = smaj
        self.smin = smin
        self.brg = brg
        self.vcv = vcv
        self.covariances = covariances

    def xyz(self):
        return gc.llh2xyz(self.lat, self.lon, self.ehgt)

    def grid(self):
        return gc.geo2grid(self.lat, self.lon)


class Switches(object):
    def __init__(self, stns=False, msrs=False, header=False, stn_msr_count=False):
        self.stns = stns
        self.msrs = msrs
        self.header = header
        self.stn_msr_count = stn_msr_count

    def reset(self):
        self.stns = False
        self.msrs = False
        self.header = False
        self.stn_msr_count = False


# ----------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------

def read_coord_elements(line, coord_types, desc_index):
    """
    Function to read coordinate components from a line in an adj/xyz file
    :param line: coordinate line
    :param coord_types: string of coordinate types in line
    :param desc_index: position of station descriptions in line
    :return: Station object
    """

    lat, lon, ohgt, ehgt, east, north, zone, x, y, z = None, None, None, None, None, None, None, None, None, None
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
        if ct == 'E':
            east = float(results[r_count])
        if ct == 'N':
            north = float(results[r_count])
        if ct == 'z':
            zone = int(results[r_count])
        if ct == 'X':
            x = float(results[r_count])
        if ct == 'Y':
            y = float(results[r_count])
        if ct == 'Z':
            z = float(results[r_count])

        r_count += 1

    # covert coordinates if lat/lon/ehgt not reported
    if not lat or not lon:
        if x and y and z:
            lat, lon, ehgt = gc.xyz2llh(x, y, z)
            lat = gc.dec2dms(lat)
            lon = gc.dec2dms(lon)
        elif east and north and zone:
            lat, lon, *others = gc.grid2geo(zone, east, north)
            lat = gc.dec2dms(lat)
            lon = gc.dec2dms(lon)
        else:
            raise ValueError(f'Not enough coordinate types present in {coord_types} for successful coordinate read')

    # Qualities
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


def read_xyz_file(xyz_file, stns={}, adj_metadata=AdjMetadata(), file_metadata=FileMetadata()):
    """
    function to read metadata and coordinate information from a DynAdjust station file.
    Can be called standalone or as a DynaResults object method.
    :param xyz_file: DynAdjust *.xyz file
    :param stns: Dictionary of Station objects
    :param adj_metadata: Defaults to new AdjMetadata object
    :param file_metadata: Defaults to new FileMetadata object
    :return: dictionary of Station objects, AdjMetadata and FileMetada objects
    """

    with open(xyz_file, 'r') as f:
        switches = Switches(header=True)
        line_count = 0
        stn_line = False
        desc_index = None

        for line_count, line in enumerate(f):

            if 'Adjusted Coordinates' in line:
                stn_line = line_count + 5

            if stn_line:
                if line_count == stn_line - 2:
                    desc_index = line.find('Description')

                if line_count == stn_line:
                    switches.reset()
                    switches.stns = True

            if switches.header:
                adj_metadata.version = read_metadata(line, 'Version:', adj_metadata.version)
                adj_metadata.build = read_metadata(line, 'Build:', adj_metadata.build)
                adj_metadata.reference_frame = read_metadata(line, 'Reference frame:', adj_metadata.reference_frame)
                file_metadata.xyz_filename = read_metadata(line, 'File name:', file_metadata.xyz_filename)
                file_metadata.xyz_date_created = read_file_date(line, file_metadata.xyz_date_created)
                adj_metadata.epoch = read_epoch(line, adj_metadata.epoch)
                adj_metadata.geoid_model = read_metadata(line, 'Geoid model:', adj_metadata.geoid_model)

                if line[:35] == 'Station coordinate types:          ':
                    coord_types = line[35:].strip()

            if switches.stns:
                if len(line) < 20:
                    switches.reset()
                    continue

                stn_object = read_coord_elements(line, coord_types, desc_index)

                stns[stn_object.name] = stn_object

    return stns, adj_metadata, file_metadata


def read_adj_file(adj_file, stns={}, msrs=[], adj_metadata=AdjMetadata(), file_metadata=FileMetadata(),
                  adj_stats=AdjStats()):
    """
    Function to read all metadata, measurement and coordinate information from an adj file.
    can be called as a standalone function or a DynaResults method
    :param adj_file: DynAdjust *.adj filename
    :param stns: Dictionary of PynAdjust Station objects
    :param msrs: list of PynAdjust Measurement objects
    :param adj_metadata: Defaults to new AdjMetadata object
    :param file_metadata: Defaults to new FileMetadata object
    :param adj_stats: Defaults to new AdjStats object
    :return: tuple of (Dictionary of Station objects, list of Measurement objects, AdjMetadata object, FileMetadata object,
             AdjStats object)
    """

    with open(adj_file, 'r') as f:
        switches = Switches(header=True)
        stn_line = None
        msr_line = None
        msr_count = 0
        dir_set = 0
        msr_header_line = None

        for line_count, line in enumerate(f):

            if 'Adjusted Measurements' in line:
                msr_line = line_count + 5
                msr_header_line = line_count + 3
                three_line_count = 0

            if msr_line:
                if line_count == msr_line:
                    switches.reset()
                    switches.msrs = True

            if msr_header_line:
                if line_count == msr_header_line:
                    tstat_switch, msr_id_switch = read_msr_fields(line)

            # don't bother reading stations from file if already present
            if not stns:
                if 'Adjusted Coordinates' in line:
                    stn_line = line_count + 5
                    print('skipped stn read')

                if stn_line:
                    if line_count == stn_line - 2:
                        desc_index = line.find('Description')

                    if line_count == stn_line:
                        switches.reset()
                        switches.stns = True

            if switches.header:
                # perhaps these function calls could be improved?
                adj_metadata.version = read_metadata(line, 'Version:', adj_metadata.version)
                adj_metadata.reference_frame = read_metadata(line, 'Reference frame:', adj_metadata.reference_frame)
                adj_metadata.epoch = read_epoch(line, adj_metadata.epoch)
                adj_metadata.geoid_model = read_metadata(line, 'Geoid model:', adj_metadata.geoid_model)

                file_metadata.adj_filename = read_metadata(line, 'File name:', file_metadata.adj_filename)
                file_metadata.adj_date_created = read_file_date(line, file_metadata.adj_date_created)

                adj_stats.soln_type = read_metadata(line, 'SOLUTION', adj_stats.soln_type)
                adj_stats.parameters = read_metadata(line, 'Number of unknown parameters', adj_stats.parameters)
                adj_stats.degrees_of_freedom = read_metadata(line, 'Degrees of freedom', adj_stats.degrees_of_freedom)
                adj_stats.chi_squared = read_metadata(line, 'Chi squared', adj_stats.chi_squared)
                adj_stats.sigma_zero = read_metadata(line, 'Rigorous Sigma Zero', adj_stats.sigma_zero)
                adj_stats.confidence_interval = read_metadata(line, 'Test confidence interval:', adj_stats.confidence_interval)

                # more complex reads have their own if statements.
                if line[:35] == 'Total time                         ':
                    items = line[35:].strip().split(':')
                    hour = int(items[0])
                    minute = int(items[1])
                    second = float(items[2])
                    adj_stats.run_time = datetime.timedelta(hours=hour, minutes=minute, seconds=second)

                if line[:35] == 'Number of measurements             ':
                    items = line[35:].split()
                    if len(items) > 1:
                        if '(' in items[1]:
                            outliers = int(items[1].replace('(', ''))
                    else:
                        outliers = 0
                    adj_stats.msr_count = int(items[0])
                    adj_stats.outlier_count = outliers

                if line[:35] == 'Global (Pelzer) Reliability        ':
                    items = line[35:].strip().split()
                    try:
                        adj_stats.global_pelzer = float(items[0])
                    except ValueError:
                        adj_stats.global_pelzer = float('NaN')

                if line[:35] == 'Station coordinate types:          ':
                    coord_types = line[35:].strip()

                # last header category
                if 'Chi-Square test' in line[:35]:
                    items = line[35:].strip().split()
                    adj_stats.chi_square_lower = float(items[0])
                    adj_stats.chi_square_upper = float(items[4])
                    adj_stats.chi_square_result = items[6].strip()

                    switches.reset()

            # read station coords (if not completed earlier)
            if switches.stns:
                if len(line) < 20:
                    switches.reset()
                    continue

                stn_object = read_coord_elements(line, coord_types, desc_index)
                stns[stn_object.name] = stn_object

            # read measurements
            if switches.msrs:
                if line.strip() == '':
                    switches.reset()
                    continue

                line_msr = read_msr_line(line, tstat_switch, msr_id_switch)

                if line_msr.msr_type == 'D':
                    msr_object = Measurement(
                        msr_type=line_msr.msr_type,
                        stn1=line_msr.stn1,
                        stn2=line_msr.stn2,
                        stn3=[],
                        ignore=[line_msr.ignore],
                        cardinal=[line_msr.cardinal],
                        msr=[],
                        adj=[],
                        cor=[],
                        msr_sd=[],
                        adj_sd=[],
                        cor_sd=[],
                        nstat=[],
                        tstat=[],
                        pelzer=[],
                        pre_adj_cor=[],
                        outlier=[],
                        msr_id=[line_msr.msr_id],
                        cluster_id=[line_msr.cluster_id]
                    )

                    msrs.append(msr_object)

                # Direction Set pointings - append attributes to most recent msr
                elif line_msr.msr_type == ' ':
                    msrs[-1].stn3.append(line_msr.stn3)
                    msrs[-1].ignore.append(line_msr.ignore)
                    msrs[-1].msr.append(str_to_dms_angle(line_msr.msr))
                    msrs[-1].adj.append(str_to_dms_angle(line_msr.adj))
                    msrs[-1].cor.append(float(line_msr.cor))
                    msrs[-1].msr_sd.append(float(line_msr.msr_sd))
                    msrs[-1].adj_sd.append(float(line_msr.adj_sd))
                    msrs[-1].cor_sd.append(float(line_msr.cor_sd))
                    msrs[-1].nstat.append(float(line_msr.nstat))
                    msrs[-1].tstat.append(float(line_msr.tstat) if line_msr.tstat else None)
                    msrs[-1].pelzer.append(float(line_msr.pelzer))
                    msrs[-1].pre_adj_cor.append(float(line_msr.pre_adj_cor))
                    msrs[-1].outlier.append(line_msr.outlier)
                    msrs[-1].msr_id.append(line_msr.msr_id)
                    msrs[-1].cluster_id.append(line_msr.cluster_id)

                # Type G/X/Y measurements (split over 3 lines)
                elif line_msr.msr_type in three_line_msrs:

                    if line_msr.cardinal == 'X' or line_msr.cardinal == 'e':
                        # create new msr object for each new triplet, otherwise append to most recent msr
                        msrs.append(
                            Measurement(
                                msr_type=line_msr.msr_type,
                                stn1 = line_msr.stn1,
                                stn2 = line_msr.stn2,
                                ignore = line_msr.ignore,
                                cardinal=[line_msr.cardinal],
                                msr=[float(line_msr.msr)],
                                adj=[float(line_msr.adj)],
                                cor=[float(line_msr.cor)],
                                msr_sd=[float(line_msr.msr_sd)],
                                adj_sd=[float(line_msr.adj_sd)],
                                cor_sd=[float(line_msr.cor_sd)],
                                nstat=[float(line_msr.nstat)],
                                tstat=[float(line_msr.tstat) if line_msr.tstat else None],
                                pelzer=[float(line_msr.pelzer)],
                                pre_adj_cor=[float(line_msr.pre_adj_cor)],
                                outlier=[line_msr.outlier],
                                msr_id=line_msr.msr_id,
                                cluster_id=line_msr.cluster_id
                            )
                        )

                    else:
                        msrs[-1].cardinal.append(line_msr.cardinal)
                        msrs[-1].msr.append(float(line_msr.msr))
                        msrs[-1].adj.append(float(line_msr.adj))
                        msrs[-1].cor.append(float(line_msr.cor))
                        msrs[-1].msr_sd.append(float(line_msr.msr_sd))
                        msrs[-1].adj_sd.append(float(line_msr.adj_sd))
                        msrs[-1].cor_sd.append(float(line_msr.cor_sd))
                        msrs[-1].nstat.append(float(line_msr.nstat))
                        msrs[-1].tstat.append(float(line_msr.tstat) if line_msr.tstat else None)
                        msrs[-1].pelzer.append(float(line_msr.pelzer))
                        msrs[-1].pelzer.append(float(line_msr.pelzer))
                        msrs[-1].pre_adj_cor.append(float(line_msr.pre_adj_cor))
                        msrs[-1].outlier.append(line_msr.outlier)

                # all other msr types
                else:
                    if line_msr.msr_type in angle_msrs:
                        msr = str_to_dms_angle(line_msr.msr)
                        adj = str_to_dms_angle(line_msr.adj)
                    elif line_msr.msr_type in hp_msrs:
                        msr = gc.hp2dec(line_msr.msr)
                        adj = gc.hp2dec(line_msr.adj)
                    else:
                        msr = float(line_msr.msr)
                        adj = float(line_msr.adj)

                    msrs.append(
                        Measurement(
                            msr_type=line_msr.msr_type,
                            stn1=line_msr.stn1,
                            stn2=line_msr.stn2,
                            stn3=line_msr.stn3,
                            ignore=line_msr.ignore,
                            cardinal=line_msr.cardinal,
                            msr=msr,
                            adj=adj,
                            cor=float(line_msr.cor),
                            msr_sd=float(line_msr.msr_sd),
                            adj_sd=float(line_msr.adj_sd),
                            cor_sd=float(line_msr.cor_sd),
                            nstat=float(line_msr.nstat),
                            tstat=float(line_msr.nstat) if line_msr.tstat else None,
                            pelzer=float(line_msr.pelzer),
                            pre_adj_cor=float(line_msr.pre_adj_cor),
                            outlier=line_msr.outlier,
                            msr_id=line_msr.msr_id,
                            cluster_id=line_msr.cluster_id
                        )
                    )

    return stns, msrs, adj_stats, adj_metadata, file_metadata


def read_apu_file(apu_file, stns={}, adj_metadata=AdjMetadata(), file_metadata=FileMetadata(), adj_stats=AdjStats()):
    """
    Function to consume DynAdjust *.apu file
    note: Variance matrix components are stored in XYZ rotation ONLY.
              - ENU Variance matrix components are rotated back to XYZ.
              - ENU Covariance block matrix components are rotated back to XYZ, at the position of the 'first' station.
    Can be called standalone, or from DynaResults object method
    :param apu_file: input *.apu file
    :param adj_metadata: defaults to new AdjMetadata object
    :param file_metadata: defaults to new FileMetadata object
    :param adj_stats: defaults to new AdjStats object
    :return: Updated dictionary of stations, AdjMetadata, FileMetadata and AdjStats objects
    """

    with open(apu_file) as f:
        switches = Switches(header=True)
        rotate_vcv = False
        pu_line = None
        cov_count = 0

        for line_count, line in enumerate(f):

            if switches.header:
                adj_metadata.version = read_metadata(line, 'Version:', adj_metadata.version)
                adj_metadata.vcv_blocks = read_metadata_tf(line, 'Stations printed in blocks', adj_metadata.vcv_blocks)
                adj_metadata.full_covariance_matrix = read_metadata_tf(line, 'Full covariance matrix',
                                                                  adj_metadata.full_covariance_matrix)
                adj_metadata.variance_matrix_units = read_metadata(line, 'Variance matrix units ', adj_metadata.variance_matrix_units)

                adj_stats.confidence_interval = read_metadata(line, 'PU confidence interval',
                                                              adj_stats.confidence_interval)
                file_metadata.apu_filename = read_metadata(line, 'File name:', file_metadata.apu_filename)
                file_metadata.apu_date_created = read_file_date(line, file_metadata.apu_date_created)

                if "Positional uncertainty of adjusted station coordinates" in line:
                    pu_line = line_count + 5
                    switches.reset()
                    rotate_vcv = True if adj_metadata.variance_matrix_units == 'ENU' else False

            if pu_line:
                if line_count == pu_line:
                    switches.stns = True

            if switches.stns:
                if line == '\n':
                    continue

                if 'Block ' in line:
                    if len(line) < 30:
                        continue

                if '-' * 30 in line:
                    continue

                # count number of columns, accounting for station names with spaces present
                num_cols = len(line[20:].split())
                if line[:20].strip() != '':
                    num_cols += 1

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
                    try:
                        vcv_22 = float(line[131:150].strip())
                        vcv_23 = float(line[150:].strip())
                    except ValueError:
                        print(line)
                        exit()
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
                    try:
                        stns[stn].hpu = hpu
                        stns[stn].vpu = vpu
                        stns[stn].vpu = vpu
                        stns[stn].smaj = smaj
                        stns[stn].smin = smin
                        stns[stn].brg = brg
                        stns[stn].vcv = vcv
                    except KeyError:
                        stns[stn] = Station(
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

                # covariance block information
                if adj_metadata.full_covariance_matrix:
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

    return stns, adj_metadata, file_metadata, adj_stats


def write_msr_shapefile(stns, msrs, network_name, ref_frame):
    """
    Function to write shapefiles using the adj results. Results are written to a 'Shp' sub-directory.
    :param stns: Dictionary of Station objects
    :param msrs: List of Measurement objects
    :param network_name: String. Name of network - must not contain full stops
    :param ref_frame: String. reference frame (for prj file creation)
    """

    def find_cardinal(msrs_list):
        """
        function to determine 'cardinal' directions of DynAdjust GNSS measurements
        :param msrs_list: list of Measurement objects (of common msr_type)
        :return: string of cardinal directions: enu/XYZ
        """

        return 'enu' if msrs_list[0].cardinal == ['e', 'n', 'u'] else 'XYZ'

    def max_stat(stat_list):
        """
        function to determine the maximum absolute value in a list
        :param stat_list: list of floats
        :return: maximum value
        """
        max_stat = None
        for s in stat_list:
            if not max_stat:
                max_stat = s
            elif abs(s) > abs(max_stat):
                max_stat = s

        return max_stat

    # create shapefiles for each msr type
    for l in msr_types:
        result = [m for m in msrs if m.msr_type == l]

        shp_name = network_name + '_' + l

        if result:
            # Direction sets
            if l in multi_line_msrs:
                write_prj(shp_name, ref_frame)

                w = shapefile.Writer(shp_name, shapeType=3)
                w.autoBalance = 1

                w.field('Stn1', 'C', size=20)
                w.field('Stn2', 'C', size=20)
                w.field('msr', 'C')
                w.field('adj', 'C')
                w.field('cor', 'N', decimal=4)
                w.field('msr_sd', 'N', decimal=4)
                w.field('adj_sd', 'N', decimal=4)
                w.field('cor_sd', 'N', decimal=4)
                w.field('nstat', 'N', decimal=4)
                w.field('msr_id', 'N')
                w.field('epoch', 'C', size=10)
                w.field('source', 'C')

                for r in result:
                    w.line([
                        [[stns[r.stn1].lon.dec(), stns[r.stn1].lat.dec()],
                         [stns[r.stn2].lon.dec(), stns[r.stn2].lat.dec()]],
                    ])

                    w.record(r.stn1, r.stn2, 'set zero', None, None, None, None, None, None,
                             r.msr_id[0], r.epoch, r.source)

                    for s in range(len(r.stn3)):
                        msr = '{:3d} {:2d} {:7.4f}'.format(r.msr[s].degree, r.msr[s].minute, r.msr[s].second)
                        adj = '{:3d} {:2d} {:7.4f}'.format(r.adj[s].degree, r.adj[s].minute, r.adj[s].second)

                        w.line([
                            [[stns[r.stn1].lon.dec(), stns[r.stn1].lat.dec()],
                             [stns[r.stn3[s]].lon.dec(), stns[r.stn3[s]].lat.dec()]],
                        ])

                        # note: msr_ids are indexed 1 after other fields (1st is for ref direction)
                        w.record(r.stn1, r.stn3[s], msr, adj, r.cor[s], r.msr_sd[s], r.adj_sd[s], r.cor_sd[s],
                                 r.nstat[s], r.msr_id[s + 1], r.epoch, r.source)

                w.close()

            # H/I/J/P/Q/R/Y msrs
            elif l in one_stn_msrs:
                write_prj(shp_name, ref_frame)

                w = shapefile.Writer(shp_name, shapeType=1)
                w.autoBalance = 1

                if l in 'IJPQ':
                    msr_dec = 10
                else:
                    msr_dec = 4

                if l not in three_line_msrs:
                    w.field('Stn1', 'C', size=20)
                    w.field('msr', 'N', decimal=msr_dec)
                    w.field('adj', 'N', decimal=msr_dec)
                    w.field('cor', 'N', decimal=4)
                    w.field('msr_sd', 'N', decimal=4)
                    w.field('adj_sd', 'N', decimal=4)
                    w.field('cor_sd', 'N', decimal=4)
                    w.field('nstat', 'N', decimal=4)
                    w.field('msr_id', 'N')
                    w.field('epoch', 'C', size=10)
                    w.field('source', 'C')

                    for r in result:

                        if l in 'PQ':
                            w.point(stns[r.stn1].lon.dec(), stns[r.stn1].lat.dec())
                            w.record(r.stn1, r.msr.dec(), r.adj.dec(), r.cor, r.msr_sd, r.adj_sd, r.cor_sd, r.nstat,
                                     r.msr_id, r.epoch, r.source)
                        else:
                            w.point(stns[r.stn1].lon.dec(), stns[r.stn1].lat.dec())
                            w.record(r.stn1, r.msr, r.adj, r.cor, r.msr_sd, r.adj_sd, r.cor_sd, r.nstat,
                                     r.msr_id, r.epoch, r.source)

                else:
                    cardinal = find_cardinal(result)

                    w.field('Stn1', 'C', size=20)
                    w.field('msr_' + cardinal[0], 'N', decimal=4)
                    w.field('msr_' + cardinal[1], 'N', decimal=4)
                    w.field('msr_' + cardinal[2], 'N', decimal=4)
                    w.field('adj_' + cardinal[0], 'N', decimal=4)
                    w.field('adj_' + cardinal[1], 'N', decimal=4)
                    w.field('adj_' + cardinal[2], 'N', decimal=4)
                    w.field('cor_' + cardinal[0], 'N', decimal=4)
                    w.field('cor_' + cardinal[1], 'N', decimal=4)
                    w.field('cor_' + cardinal[2], 'N', decimal=4)
                    w.field('msr_sd_' + cardinal[0], 'N', decimal=4)
                    w.field('msr_sd_' + cardinal[1], 'N', decimal=4)
                    w.field('msr_sd_' + cardinal[2], 'N', decimal=4)
                    w.field('adj_sd_' + cardinal[0], 'N', decimal=4)
                    w.field('adj_sd_' + cardinal[1], 'N', decimal=4)
                    w.field('adj_sd_' + cardinal[2], 'N', decimal=4)
                    w.field('cor_sd_' + cardinal[0], 'N', decimal=4)
                    w.field('cor_sd_' + cardinal[1], 'N', decimal=4)
                    w.field('cor_sd_' + cardinal[2], 'N', decimal=4)
                    w.field('nstat_' + cardinal[0], 'N', decimal=4)
                    w.field('nstat_' + cardinal[1], 'N', decimal=4)
                    w.field('nstat_' + cardinal[2], 'N', decimal=4)
                    w.field('max_cor', 'N', decimal=4)
                    w.field('max_nstat', 'N', decimal=4)
                    w.field('msr_id', 'N')
                    w.field('epoch', 'C', size=10)
                    w.field('source', 'C')

                    for r in result:
                        max_nstat = max_stat(r.nstat)
                        max_cor = max_stat(r.cor)

                        w.point(stns[r.stn1].lon.dec(), stns[r.stn1].lat.dec())

                        w.record(r.stn1,
                                 r.msr[0], r.msr[1], r.msr[2],
                                 r.adj[0], r.adj[1], r.adj[2],
                                 r.cor[0], r.cor[1], r.cor[2],
                                 r.msr_sd[0], r.msr_sd[1], r.msr_sd[2],
                                 r.adj_sd[0], r.adj_sd[1], r.adj_sd[2],
                                 r.cor_sd[0], r.cor_sd[1], r.cor_sd[2],
                                 r.nstat[0], r.nstat[1], r.nstat[2],
                                 max_cor, max_nstat,
                                 r.msr_id, r.epoch, r.source
                                 )

                w.close()

            # B/C/E/G/K/L/M/S/V/X/Z
            elif l in two_stn_msrs:
                write_prj(shp_name, ref_frame)

                w = shapefile.Writer(shp_name, shapeType=3)
                w.autoBalance = 1

                # B/C/E/K/L/M/S/V/Z
                if l not in three_line_msrs:
                    if l in angle_msrs:
                        w.field('Stn1', 'C', size=20)
                        w.field('Stn2', 'C', size=20)
                        w.field('msr', 'C')
                        w.field('adj', 'C')
                        w.field('cor', 'N', decimal=4)
                        w.field('msr_sd', 'N', decimal=4)
                        w.field('adj_sd', 'N', decimal=4)
                        w.field('cor_sd', 'N', decimal=4)
                        w.field('nstat', 'N', decimal=4)
                        w.field('msr_id', 'N')
                        w.field('epoch', 'C', size=10)
                        w.field('source', 'C')

                        for r in result:
                            w.line([
                                [[stns[r.stn1].lon.dec(), stns[r.stn1].lat.dec()],
                                 [stns[r.stn2].lon.dec(), stns[r.stn2].lat.dec()]],
                            ])

                            msr = '{:3d} {:2d} {:7.4f}'.format(r.msr.degree, r.msr.minute, r.msr.second)
                            adj = '{:3d} {:2d} {:7.4f}'.format(r.adj.degree, r.adj.minute, r.adj.second)

                            w.record(r.stn1, r.stn2, msr, adj, r.cor, r.msr_sd, r.adj_sd, r.cor_sd, r.nstat,
                                     r.msr_id, r.epoch, r.source)

                    else:
                        w.field('Stn1', 'C', size=20)
                        w.field('Stn2', 'C', size=20)
                        w.field('msr', 'N', decimal=4)
                        w.field('adj', 'N', decimal=4)
                        w.field('cor', 'N', decimal=4)
                        w.field('msr_sd', 'N', decimal=4)
                        w.field('adj_sd', 'N', decimal=4)
                        w.field('cor_sd', 'N', decimal=4)
                        w.field('nstat', 'N', decimal=4)
                        w.field('msr_id', 'N')
                        w.field('epoch', 'C', size=10)
                        w.field('source', 'C')

                        for r in result:
                            w.line([
                                [[stns[r.stn1].lon.dec(), stns[r.stn1].lat.dec()],
                                 [stns[r.stn2].lon.dec(), stns[r.stn2].lat.dec()]],
                            ])

                            w.record(r.stn1, r.stn2, r.msr, r.adj, r.cor, r.msr_sd, r.adj_sd, r.cor_sd, r.nstat,
                                     r.msr_id, r.epoch, r.source)

                # G/X
                else:
                    cardinal = find_cardinal(result)

                    w.field('Stn1', 'C', size=20)
                    w.field('Stn2', 'C', size=20)
                    w.field('msr_' + cardinal[0], 'N', decimal=4)
                    w.field('msr_' + cardinal[1], 'N', decimal=4)
                    w.field('msr_' + cardinal[2], 'N', decimal=4)
                    w.field('adj_' + cardinal[0], 'N', decimal=4)
                    w.field('adj_' + cardinal[1], 'N', decimal=4)
                    w.field('adj_' + cardinal[2], 'N', decimal=4)
                    w.field('cor_' + cardinal[0], 'N', decimal=4)
                    w.field('cor_' + cardinal[1], 'N', decimal=4)
                    w.field('cor_' + cardinal[2], 'N', decimal=4)
                    w.field('msr_sd_' + cardinal[0], 'N', decimal=4)
                    w.field('msr_sd_' + cardinal[1], 'N', decimal=4)
                    w.field('msr_sd_' + cardinal[2], 'N', decimal=4)
                    w.field('adj_sd_' + cardinal[0], 'N', decimal=4)
                    w.field('adj_sd_' + cardinal[1], 'N', decimal=4)
                    w.field('adj_sd_' + cardinal[2], 'N', decimal=4)
                    w.field('cor_sd_' + cardinal[0], 'N', decimal=4)
                    w.field('cor_sd_' + cardinal[1], 'N', decimal=4)
                    w.field('cor_sd_' + cardinal[2], 'N', decimal=4)
                    w.field('nstat_' + cardinal[0], 'N', decimal=4)
                    w.field('nstat_' + cardinal[1], 'N', decimal=4)
                    w.field('nstat_' + cardinal[2], 'N', decimal=4)
                    w.field('max_cor', 'N', decimal=4)
                    w.field('max_nstat', 'N', decimal=4)
                    w.field('msr_id', 'N')
                    w.field('epoch', 'C', size=10)
                    w.field('source', 'C')

                    for r in result:
                        max_nstat = max_stat(r.nstat)
                        max_cor = max_stat(r.cor)

                        w.line([
                            [[stns[r.stn1].lon.dec(), stns[r.stn1].lat.dec()],
                             [stns[r.stn2].lon.dec(), stns[r.stn2].lat.dec()]],
                        ])

                        w.record(r.stn1, r.stn2,
                                 r.msr[0], r.msr[1], r.msr[2],
                                 r.adj[0], r.adj[1], r.adj[2],
                                 r.cor[0], r.cor[1], r.cor[2],
                                 r.msr_sd[0], r.msr_sd[1], r.msr_sd[2],
                                 r.adj_sd[0], r.adj_sd[1], r.adj_sd[2],
                                 r.cor_sd[0], r.cor_sd[1], r.cor_sd[2],
                                 r.nstat[0], r.nstat[1], r.nstat[2],
                                 max_cor, max_nstat,
                                 r.msr_id, r.epoch, r.source
                                 )

                w.close()

            # A (not D)
            elif l in three_stn_msrs:
                write_prj(shp_name, ref_frame)

                w = shapefile.Writer(shp_name, shapeType=3)
                w.autoBalance = 1

                w.field('Stn1', 'C', size=20)
                w.field('Stn2', 'C', size=20)
                w.field('Stn3', 'C', size=20)
                w.field('msr', 'C')
                w.field('adj', 'C')
                w.field('cor', 'N', decimal=4)
                w.field('msr_sd', 'N', decimal=4)
                w.field('adj_sd', 'N', decimal=4)
                w.field('cor_sd', 'N', decimal=4)
                w.field('nstat', 'N', decimal=4)
                w.field('msr_id', 'N')
                w.field('epoch', 'C', size=10)
                w.field('source', 'C')

                for r in result:
                    msr = '{:3d} {:2d} {:7.4f}'.format(r.msr.degree, r.msr.minute, r.msr.second)
                    adj = '{:3d} {:2d} {:7.4f}'.format(r.adj.degree, r.adj.minute, r.adj.second)

                    w.line([
                        [[stns[r.stn1].lon.dec(), stns[r.stn1].lat.dec()],
                         [stns[r.stn2].lon.dec(), stns[r.stn2].lat.dec()]],
                        [[stns[r.stn2].lon.dec(), stns[r.stn2].lat.dec()],
                         [stns[r.stn3].lon.dec(), stns[r.stn3].lat.dec()]],
                    ])

                    w.record(r.stn1, r.stn2, r.stn3, r.msr, r.adj, r.cor, r.msr_sd, r.adj_sd, r.cor_sd, r.nstat)

                w.close()

    # return to parent directory
    os.chdir('..')


def str_to_dms_angle(angle):
    """
    function to create a DMSAngle object from a string
    :param angle: string of text in format "DDD MM SS.SSSS"
    :return: DMSAngle object
    """
    items = angle.split()
    deg = int(items[0])
    min = int(items[1])
    sec = float(items[2])
    dms = gc.DMSAngle(degree=deg, minute=min, second=sec)

    return dms


def read_msr_fields(line):
    """
    Function to detect presence of tstat and msr_id fields in a DynaAdjust adj file
    :param line: Measurement header line
    :return: True/False switches for presence of tstat and msr_id fields
    """
    t_stat_field = False
    msr_id_field = False

    if 'T-stat' in line:
        t_stat_field = True

    if 'Meas. ID' in line:
        msr_id_field = True

    return t_stat_field, msr_id_field


def read_msr_line(line, tstat_switch, msr_id_switch):
    """
    Function to read msr components from a single line in adj file
    :param line: Msr line in adj file
    :param tstat_switch: True/False switch for presence of tstats
    :param msr_id_switch: True/False switch for presence of msr_ids
    :return: Dictionary of components
    """
    # initialise empty measurement, populate attributes
    line_msr = Measurement()

    line_msr.msr_type = line[:1]
    line_msr.stn1 = line[2:22].strip() if line[2:22].strip() else None
    line_msr.stn2 = line[22:42].strip() if line[22:42].strip() else None
    line_msr.stn3 = line[42:62].strip() if line[42:62].strip() else None
    if line[62:65].strip() == '*':
        line_msr.ignore == True

    if line_msr.msr_type == 'D':
        if msr_id_switch:
            items = line.split()
            line_msr.cluster_id = items[-1].strip()
            line_msr.msr_id = items[-2]

        return line_msr
    else:
        line_msr.cardinal = line[65:67].strip()
        line_msr.msr = line[67:86].strip()
        line_msr.adj = line[86:105].strip()
        line_msr.cor = line[105:117].strip()
        line_msr.msr_sd = line[117:130].strip()
        line_msr.adj_sd = line[130:143].strip()
        line_msr.cor_sd = line[143:156].strip()
        line_msr.nstat = line[156:167].strip()

        if tstat_switch:
            line_msr.tstat = line[167:178].strip()
            line_msr.pelzer = line[178:190].strip()
            line_msr.pre_adj_cor = line[190:204].strip()
            line_msr.outlier = line[204:216].strip()
            if msr_id_switch:
                line_msr.msr_id = line[216:226].strip()
                line_msr.cluster_id = line[226:236].strip()

        else:
            line_msr.pelzer = line[167:178].strip()
            line_msr.pre_adj_cor = line[180:193].strip()
            line_msr.outlier = line[204:205]
            if msr_id_switch:
                line_msr.msr_id = line[205:216].strip()
                line_msr.cluster_id = line[216:226].strip()

        return line_msr


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

        if date_str[-1] == 'M':
            file_date = datetime.datetime.strptime(date_str, '%A, %d %B %Y, %I:%M:%S %p')
        else:
            file_date = datetime.datetime.strptime(date_str, '%A, %d %B %Y, %H:%M:%S')

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


def write_stn_shapefile(stns, network_name, ref_frame=None):
    """
    function to write shapefiles from DynAdjust adj/apu/xyz results
    :param stns: dictionary of Station objects
    :param network_name: name of network (used in shapefile name output)
    :param ref_frame: optional, reference frame of coordinates
    """

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


# ----------------------------------------------------------------------
# Example use
# ----------------------------------------------------------------------
adj_file = ''
apu_file = ''
xyz_file = ''
xml_msr_file = ''
network_name = 'pynadjust_example'

if adj_file or apu_file or xyz_file:
    # consume results
    dyna_res = DynaResults(adj_file=adj_file, apu_file=apu_file, xyz_file=xyz_file)

    # allocate epoch and source fields (requires unique msr ids in xml file)
    if xml_msr_file:
        dyna_res.link_source_with_msr_id(xml_msr_file=xml_msr_file)

    # add standard type B StdDevs to all stns:
    dyna_res.add_type_b_sd(sd_e=0.006, sd_n=0.006, sd_u=0.012)

    # compute relative uncertainties between each station and its five nearest neighbours
    dyna_res.compute_rel_uncertainty(neighbours=5)

    # create shapefile
    dyna_res.write_shp(network_name=network_name)

    # print network summary + stn listing
    with open(network_name + '.txt', 'w') as f:
        line_break = '-' * 50 + '\n'
        header_str = 'PynAdjust example report\n' + line_break
        header_str += f'adj file:            {dyna_res.file_metadata.adj_filename}\n'
        header_str += f'adj date:            {dyna_res.file_metadata.adj_date_created}\n'
        header_str += f'apu file:            {dyna_res.file_metadata.apu_filename}\n'
        header_str += f'apu date:            {dyna_res.file_metadata.apu_date_created}\n'
        header_str += f'xyz file:            {dyna_res.file_metadata.xyz_filename}\n'
        header_str += f'xyz date:            {dyna_res.file_metadata.xyz_date_created}\n'
        header_str += line_break
        header_str += f'Stations:            {len(dyna_res.stns)}\n'
        header_str += f'Measurements:        {dyna_res.adj_stats.msr_count}\n'
        header_str += f'Reference Frame:     {dyna_res.adj_metadata.reference_frame}\n'
        header_str += f'Epoch:               {dyna_res.adj_metadata.epoch}\n'
        header_str += line_break
        header_str += f'Sigma Zero:          {dyna_res.adj_stats.sigma_zero}\n'
        header_str += f'Chi Square Test:     {dyna_res.adj_stats.chi_square_lower} < {dyna_res.adj_stats.chi_square_result} < {dyna_res.adj_stats.chi_square_upper}\n'
        header_str += f'Potential Outliers:  {dyna_res.adj_stats.outlier_count}\n'
        header_str += f'Confidence Interval: {dyna_res.adj_stats.confidence_interval}\n'
        header_str += f'Full Covariances:    {dyna_res.adj_metadata.full_covariance_matrix}\n'
        header_str += line_break +'\n'
        header_str += 'Station Coordinates\n' + line_break + '\n'
        header_str += f'{"Station":20s} {"Latitude":>12s} {"Longitude":>12s} {"Ehgt":>10} {"sd_e":>10s} {"sd_n":>10s} {"sd_u":>10s}\n'
        header_str += '-' * 91 + '\n'
        f.write(header_str)

        for stn in dyna_res.stns.values():
            stn_str = f'{stn.name:20s} {stn.lat.dec():12.8f} {stn.lon.dec():12.8f} {stn.ehgt:10.4f} {stn.sd_e:10.4f} {stn.sd_n:10.4f} {stn.sd_u:10.4f}\n'
            f.write(stn_str)

        if dyna_res.relative_uncertainties:
            header_str = '\n\nRelative Uncertainties\n' + line_break + '\n'
            header_str += f'{"Stn1":20s} {"Stn2":20s} {"Rel_Hz":>10s} {"Rel_Vt":>10s} {"REE_SMaj":>10s} {"REE_SMin":>10s} {"REE_Brg":>10s} {"Covariance":>11s}\n'
            header_str += '-' * 109 + '\n'
            f.write(header_str)

            for ru in dyna_res.relative_uncertainties:
                if ru.covariance:
                    cov_str = 'Yes' if ru.covariance else 'No'
                ru_str = f'{ru.stn1:20s} {ru.stn2:20s} {ru.ru_hz:10.4f} {ru.ru_vt:10.4f} {ru.ree_smaj:10.4f} {ru.ree_smin:10.4f} {ru.ree_brg:10.1f} {cov_str:>11s}\n'
                f.write(ru_str)
