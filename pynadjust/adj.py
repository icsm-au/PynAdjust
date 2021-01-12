# pynadjust adj module - retrieve data from DynAdjust.adj files

# For review/consideration:
#     - Confirm all measurement types and components read and stored correctly (test on adj will all msr types)
#     - Store Latitude/Longitude as DMSAngle objects (done here) or as Decimal degrees float
#     - Consider storing measurement objects in a dictionary (done here), or set.

import datetime
import geodepy.convert as gc
import shapefile
import os
import pynadjust.common_fn as common_fn
from pynadjust.pynadjust_classes import Station, DynaMetadata, Switches
from pynadjust.xyz import read_coord_elements, DynaXYZ
import re

# ----------------------------------------------------------------------
# sets for various msr types
# ----------------------------------------------------------------------

msr_types = ('A', 'B', 'C', 'D', 'E', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'P', 'Q', 'R', 'S', 'V', 'X', 'Y', 'Z')
three_line_msrs = ('G', 'X', 'Y')
multi_line_msrs = ('D', ' ', )
one_stn_msrs = ('H', 'I', 'J', 'P', 'Q', 'R', 'Y')
two_stn_msrs = ('B', 'C', 'E', 'G', 'K', 'L', 'M', 'S', 'V', 'X', 'Z')
three_stn_msrs = ('A', 'D')
angle_msrs = ('A', 'B', 'D', 'V', 'Z')  # angle msr types usually expressed in "DDD MM SS.SSSS" format (not DEC or HP)
hp_msrs = ('I', 'J', 'P', 'Q')


# ----------------------------------------------------------------------
# Classes and functions
# ----------------------------------------------------------------------


class DynaAdj(object):
    def __init__(self, filename, stns=None, metadata=None):
        """
        DynaAdj object which reads all metadata, measurement and coordinate information upon initialisation
        :param filename: name of adj file to be consumed.
        :param stns: optional argument to use existing dictionary of station objects (e.g. following xyz read)
        :param metadata: optional argument to use metadata object (e.g. following xyz read)
        """

        def read_adj(adj_file, read_stns, read_metadata):
            """
            Function/method to read all metadata, measurement and coordinate information from an adj file.
            :param adj_file: name of adj file to be consumed.
            :param read_stns: True/False switch to read station results in file
            :param read_metadata: True/False switch to read metadata in file
            :return: components of DynaAdj object for initialisation
            """

            stns = {}
            msrs = {}

            version = None
            reference_frame = None
            file_name = None
            file_date = None
            epoch = None
            geoid_model = None
            soln_type = None
            run_time = None
            parameters = None
            msr_count = None
            outlier_count = None
            degrees_of_freedom = None
            chi_squared = None
            sigma_zero = None
            global_pelzer = None
            chi_square_lower = None
            chi_square_upper = None
            chi_square_result = None

            with open(adj_file, 'r') as adj_fh:
                switches = Switches()
                switches.header = True
                stn_line = None
                mandatory_coord_types = 'PLHh'

                msr_line = None
                msr_count = 0
                dir_set = 0
                msr_header_line = None

                for line_count, line in enumerate(adj_fh):

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

                    if read_stns:
                        if 'Adjusted Coordinates' in line:
                            stn_line = line_count + 5

                        if stn_line:
                            if line_count == stn_line - 2:
                                desc_index = line.find('Description')

                            if line_count == stn_line:
                                switches.reset()
                                switches.stns = True

                    if switches.header:
                        if read_metadata:
                            version = common_fn.read_metadata(line, 'Version:', version)
                            reference_frame = common_fn.read_metadata(line, 'Reference frame:', reference_frame)
                            epoch = common_fn.read_epoch(line, epoch)
                            geoid_model = common_fn.read_metadata(line, 'Geoid model:', geoid_model)

                        # simple text reads (and more complex, but common, types) are done by function
                        file_name = common_fn.read_metadata(line, 'File name:', file_name)
                        file_date = common_fn.read_file_date(line, file_date)
                        soln_type = common_fn.read_metadata(line, 'SOLUTION', soln_type)
                        parameters = common_fn.read_metadata(line, 'Number of unknown parameters', parameters)
                        degrees_of_freedom = common_fn.read_metadata(line, 'Degrees of freedom', degrees_of_freedom)
                        chi_squared = common_fn.read_metadata(line, 'Chi squared', chi_squared)
                        sigma_zero = common_fn.read_metadata(line, 'Rigorous Sigma Zero', sigma_zero)

                        # more complex reads have their own if statements.
                        if line[:35] == 'Total time                         ':
                            items = line[35:].strip().split(':')
                            hour = int(items[0])
                            minute = int(items[1])
                            second = float(items[2])
                            run_time = datetime.timedelta(hours=hour, minutes=minute, seconds=second)

                        if line[:35] == 'Number of measurements             ':
                            items = line[35:].split()
                            if len(items) > 1:
                                if '(' in items[1]:
                                    outliers = int(items[1].replace('(', ''))
                            else:
                                outliers = 0
                            msr_count = int(items[0])
                            outlier_count = outliers

                        if line[:35] == 'Global (Pelzer) Reliability        ':
                            items = line[35:].strip().split()
                            try:
                                global_pelzer = float(items[0])
                            except ValueError:
                                global_pelzer = float('NaN')

                        if line[:35] == 'Station coordinate types:          ':
                            coord_types = line[35:].strip()
                            missing_coord_types = set()
                            for l in mandatory_coord_types:
                                if l not in coord_types:
                                    missing_coord_types.add(l)
                            if missing_coord_types:
                                raise ValueError(f'Mandatory coordinate types {missing_coord_types} '
                                                 f'not present in {adj_file}')

                        # last header category, write to metadata object
                        if line[:35] == 'Chi-Square test (95.0%)            ':
                            items = line[35:].strip().split()
                            chi_square_lower = float(items[0])
                            chi_square_upper = float(items[4])
                            chi_square_result = items[6].strip()

                            metadata = DynaMetadata(epoch=epoch, reference_frame=reference_frame,
                                                    geoid_model=geoid_model, version=version)
                            adj_stats = DynaAdj.AdjStats(
                                soln_type=soln_type,
                                run_time=run_time,
                                parameters=parameters,
                                msr_count=msr_count,
                                degrees_of_freedom=degrees_of_freedom,
                                sigma_zero=sigma_zero,
                                chi_squared=chi_squared,
                                outlier_count=outlier_count,
                                global_pelzer=global_pelzer,
                                chi_square_lower=chi_square_lower,
                                chi_square_upper=chi_square_upper,
                                chi_square_result=chi_square_result
                            )

                            switches.reset()

                    if switches.stns:
                        if len(line) < 20:
                            switches.reset()
                            continue

                        # store details as a Station object, leaving uncertainty fields null
                        stn_object = read_coord_elements(line, coord_types, desc_index)
                        stns[stn_object.name] = stn_object

                    if switches.msrs:
                        if line.strip() == '':
                            switches.reset()
                            continue

                        msr_components = read_msr_line(line, tstat_switch, msr_id_switch)
                        msr_type = msr_components['msr_type']

                        # Direction sets - initialise new msr object
                        if msr_type == 'D':
                            msr_count += 1
                            dir_set += 1
                            d_stn1 = msr_components['stn1']
                            d_stn2 = msr_components['stn2']
                            d_msr_id = [msr_components['msr_id']]

                            msr_object = DynaAdj.Measurement(
                                msr_type=msr_type,
                                stn1=d_stn1,
                                stn2=d_stn2,
                                stn3=[],
                                ignore=[],
                                cardinal=[],
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
                                msr_id=d_msr_id,
                                cluster_id=[]
                            )

                            msrs[msr_count] = msr_object

                        # Direction Set pointings - append msrs to lists
                        elif msr_type == ' ':
                            msrs[msr_count].stn3.append(msr_components['stn3'])
                            msrs[msr_count].ignore.append(msr_components['ignore'])
                            msrs[msr_count].msr.append(str_to_dms_angle(msr_components['msr']))
                            msrs[msr_count].adj.append(str_to_dms_angle(msr_components['adj']))
                            msrs[msr_count].cor.append(float(msr_components['cor']))
                            msrs[msr_count].msr_sd.append(float(msr_components['msr_sd']))
                            msrs[msr_count].adj_sd.append(float(msr_components['adj_sd']))
                            msrs[msr_count].cor_sd.append(float(msr_components['cor_sd']))
                            msrs[msr_count].nstat.append(float(msr_components['nstat']))
                            if tstat_switch:
                                msrs[msr_count].tstat.append(float(msr_components['tstat']))
                            else:
                                msrs[msr_count].tstat.append(None)
                            msrs[msr_count].pelzer.append(float(msr_components['pelzer']))
                            msrs[msr_count].pre_adj_cor.append(float(msr_components['pre_adj_cor']))
                            msrs[msr_count].outlier.append(msr_components['outlier'])
                            msrs[msr_count].msr_id.append(msr_components['msr_id'])
                            msrs[msr_count].cluster_id.append(msr_components['cluster_id'])

                        # Type G/X/Y measurements (split over 3 lines)
                        elif msr_type in three_line_msrs:
                            three_line_count += 1

                            if three_line_count % 3 == 1:
                                # re-initialise lists at new msr
                                cardinal_list = []
                                msr_list = []
                                adj_list = []
                                cor_list = []
                                msr_sd_list = []
                                adj_sd_list = []
                                cor_sd_list = []
                                nstat_list = []
                                tstat_list = []
                                pelzer_list = []
                                pre_adj_cor_list = []
                                outlier_list = []
                                msr_id = msr_components['msr_id']
                                cluster_id = msr_components['cluster_id']

                            stn1 = msr_components['stn1']

                            if msr_components['stn2'] != '':
                                stn2 = msr_components['stn2']
                            else:
                                stn2 = None

                            ignore = msr_components['ignore']

                            cardinal_list.append(msr_components['cardinal'])
                            msr_list.append(float(msr_components['msr']))
                            adj_list.append(float(msr_components['adj']))
                            cor_list.append(float(msr_components['cor']))
                            msr_sd_list.append(float(msr_components['msr_sd']))
                            adj_sd_list.append(float(msr_components['adj_sd']))
                            cor_sd_list.append(float(msr_components['cor_sd']))
                            nstat_list.append(float(msr_components['nstat']))
                            if tstat_switch:
                                tstat_list.append(float(msr_components['tstat']))
                            else:
                                tstat_list.append(None)
                            pelzer_list.append(float(msr_components['pelzer']))
                            pre_adj_cor_list.append(float(msr_components['pre_adj_cor']))
                            outlier_list.append(msr_components['outlier'])

                            if (three_line_count % 3) == 0:
                                msr_count += 1

                                msr_object = DynaAdj.Measurement(
                                    msr_type=msr_type,
                                    stn1=stn1,
                                    stn2=stn2,
                                    stn3=None,
                                    ignore=ignore,
                                    cardinal=cardinal_list,
                                    msr=msr_list,
                                    adj=adj_list,
                                    cor=cor_list,
                                    msr_sd=msr_sd_list,
                                    adj_sd=adj_sd_list,
                                    cor_sd=cor_sd_list,
                                    nstat=nstat_list,
                                    tstat=tstat_list,
                                    pelzer=pelzer_list,
                                    pre_adj_cor=pre_adj_cor_list,
                                    outlier=outlier_list,
                                    msr_id=msr_id,
                                    cluster_id=cluster_id
                                )

                                msrs[msr_count] = msr_object
                            else:
                                pass

                        # all other msr types
                        else:
                            msr_count += 1
                            if msr_type in angle_msrs:
                                msr = str_to_dms_angle(msr_components['msr'])
                                adj = str_to_dms_angle(msr_components['adj'])
                            elif msr_type in hp_msrs:
                                msr = gc.hp2dec(msr_components['msr'])
                                adj = gc.hp2dec(msr_components['adj'])
                            else:
                                msr = float(msr_components['msr'])
                                adj = float(msr_components['adj'])

                            if tstat_switch:
                                tstat = float(msr_components['tstat'])
                            else:
                                tstat = None

                            msr_object = DynaAdj.Measurement(
                                msr_type=msr_type,
                                stn1=msr_components['stn1'],
                                stn2=msr_components['stn2'],
                                stn3=msr_components['stn3'],
                                ignore=msr_components['ignore'],
                                cardinal=msr_components['cardinal'],
                                msr=msr,
                                adj=adj,
                                cor=float(msr_components['cor']),
                                msr_sd=float(msr_components['msr_sd']),
                                adj_sd=float(msr_components['adj_sd']),
                                cor_sd=float(msr_components['cor_sd']),
                                nstat=float(msr_components['nstat']),
                                tstat=tstat,
                                pelzer=float(msr_components['pelzer']),
                                pre_adj_cor=float(msr_components['pre_adj_cor']),
                                outlier=msr_components['outlier'],
                                msr_id=msr_components['msr_id'],
                                cluster_id=msr_components['cluster_id']
                            )

                            msrs[msr_count] = msr_object

            return stns, msrs, metadata, file_name, file_date, adj_stats

        # initialise read stns/metadata switches as True
        read_stns = True
        read_metadata = True

        # turn off if stn/metadata info is present, and take current results
        if stns:
            read_stns = False
            self.stns = stns
        if metadata:
            read_metadata = False
            self.metadata = metadata

        results = read_adj(filename, read_stns=read_stns, read_metadata=read_metadata)

        # assign if not read already
        if read_stns:
            self.stns = results[0]
        if read_metadata:
            self.metadata = results[2]

        # assign rest of information to DynaAdj object
        self.msrs = results[1]
        self.file_name = results[3]
        self.file_date = results[4]
        self.adj_stats = results[5]

    class Measurement(object):
        def __init__(self, msr_type, stn1, stn2, stn3, ignore, cardinal, msr, adj, cor, msr_sd, adj_sd, cor_sd, nstat,
                     tstat, pelzer, pre_adj_cor, outlier, msr_id, cluster_id, epoch=None, source=None):
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

    def link_source_with_msr_id(self, xml_msr_file):
        """
        method to link the database id tags between msr.xml file and *.adj file to assign source and epoch to adjusted
        measurements
        :param xml_msr_file: measurement xml file
        :return:
        """

        metadata = {}

        with open(xml_msr_file, 'r') as msr_fh:
            for line in msr_fh:
                if '<DnaMeasurement>' in line:
                    epoch = None
                    source = None

                if '<Source>' in line:
                    source = re.findall("<Source>(.*?)</Source>", line)[0]
                if '<Epoch>' in line:
                    epoch = re.findall("<Epoch>(.*?)</Epoch>", line)[0]
                if '<MeasurementID>' in line:
                    msr_id = re.findall("<MeasurementID>(.*?)</MeasurementID>", line)[0]

                    if msr_id in metadata:
                        print(f' *** repeat occurance of msr_id: {msr_id}')
                        exit()
                    if msr_id != '':
                        metadata[msr_id] = {
                            'epoch': epoch,
                            'source': source,
                        }

        for m in self.msrs.values():
            msr_type = m.msr_type

            if msr_type != 'D':
                if m.msr_id in metadata:
                    m.epoch = metadata[m.msr_id]['epoch']
                    m.source = metadata[m.msr_id]['source']
            else:
                if m.msr_id:
                    if m.msr_id[0] in metadata:
                        m.epoch = metadata[m.msr_id[0]]['epoch']
                        m.source = metadata[m.msr_id[0]]['source']

    def compute_group_vf(self, source=None):
        """
        method to compute estimated variances factors for measurement groups.
        Note: this computation assumes measurements are uncorrelated. use results with caution for correlated
        measurement types e.g. X, D, etc.
        :param source: Optional => Computation only performed on measurements in a specified source
        :return:
        """

        def increment_vs2_sum(group_vf_dict, msr_type, vs2):
            """
            Function to increment the vs2 sum for group vf count
            :param group_vf_dict:
            :param msr_type:
            :param vs2:
            :return:
            """
            try:
                group_vf_dict[msr_type]['count'] += 1
                group_vf_dict[msr_type]['vs2'] += vs2
            except KeyError:
                group_vf_dict[msr_type] = {
                    'count': 1,
                    'vs2': vs2,
                }

        groups = {}

        for msr in self.msrs.values():
            if source:
                if msr.source != source:
                    continue

            msr_type = msr.msr_type
            if msr_type in three_line_msrs:
                sum_vs2 = 0
                # print(msr.cardinal)
                for i in range(0, 3):
                    msr_type = msr.msr_type + msr.cardinal[i]

                    vs2 = (msr.cor[i] / msr.msr_sd[i]) ** 2
                    sum_vs2 += vs2

                    # add individual G/X/Y components
                    increment_vs2_sum(groups, msr_type, vs2)

                # add total G/X/Y
                increment_vs2_sum(groups, msr.msr_type, sum_vs2)

            elif msr_type in multi_line_msrs:
                for i in range(len(msr.msr)):
                    vs2 = (msr.cor[i] / msr.msr_sd[i]) ** 2
                    increment_vs2_sum(groups, msr_type, vs2)

            else:
                vs2 = (msr.cor / msr.msr_sd) ** 2

                increment_vs2_sum(groups, msr_type, vs2)

        total_vs2 = 0
        group_msr_count = 0
        group_vs2 = 0
        dof = float(self.adj_stats.degrees_of_freedom)
        total_msrs = float(self.adj_stats.msr_count)
        for msr_type in groups:
            # exclude whole G/X/Y msrs. Their components are included in this summation.
            if msr_type not in three_line_msrs:
                total_vs2 += groups[msr_type]['vs2']
                group_msr_count += groups[msr_type]['count']
            else:
                groups[msr_type]['count'] *= 3

            gvf = (total_msrs * groups[msr_type]['vs2']) / (groups[msr_type]['count'] * dof)
            groups[msr_type]['gvf'] = gvf

            group_vs2 += groups[msr_type]['vs2']

        est_vf = (total_msrs / group_msr_count) * (total_vs2 / dof)

        return groups, est_vf

    class AdjStats(object):
        def __init__(self, soln_type, run_time, parameters, msr_count, degrees_of_freedom, sigma_zero, chi_squared,
                     outlier_count, global_pelzer, chi_square_lower, chi_square_upper, chi_square_result):
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


def write_adj_shapefile(adj_object, network_name):
    """
    Function to write shapefiles using the adj results. Results are written to a 'Shp' sub-directory.
    Existing results are overwritten.
    :param adj_object: DynaAdj object
    :param network_name: Name of network
    """

    def find_cardinal(msrs_list):
        """
        function to determine 'cardinal' directions of DynAdjust GNSS measurements
        :param msrs_object: list of Measurement objects
        :return: string of cardinal directions: enu/XYZ
        """
        for m in msrs_list:
            if m.cardinal == ['e', 'n', 'u']:
                cardinal = 'enu'
            else:
                cardinal = 'XYZ'
            break

        return cardinal

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

    common_fn.check_enter_dir('shp')

    if adj_object.stns:
        common_fn.write_stns_shapefile(adj_object.stns, network_name, ref_frame=adj_object.metadata.reference_frame)

    for l in msr_types:
        result = [adj_object.msrs[m] for m in adj_object.msrs if adj_object.msrs[m].msr_type == l]

        shp_name = network_name + '_' + l

        if result:
            # Direction sets
            if l in multi_line_msrs:
                common_fn.write_prj(shp_name, adj_object.metadata.reference_frame)

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
                        [[adj_object.stns[r.stn1].lon.dec(), adj_object.stns[r.stn1].lat.dec()],
                        [adj_object.stns[r.stn2].lon.dec(), adj_object.stns[r.stn2].lat.dec()]],
                    ])

                    w.record(r.stn1, r.stn2, 'set zero', None, None, None, None, None, None,
                             r.msr_id[0], r.epoch, r.source)

                    for s in range(len(r.stn3)):
                        msr = '{:3d} {:2d} {:7.4f}'.format(r.msr[s].degree, r.msr[s].minute, r.msr[s].second)
                        adj = '{:3d} {:2d} {:7.4f}'.format(r.adj[s].degree, r.adj[s].minute, r.adj[s].second)

                        w.line([
                            [[adj_object.stns[r.stn1].lon.dec(), adj_object.stns[r.stn1].lat.dec()],
                            [adj_object.stns[r.stn3[s]].lon.dec(), adj_object.stns[r.stn3[s]].lat.dec()]],
                        ])

                        # note: msr_ids are indexed 1 after other fields (1st is for ref direction)
                        w.record(r.stn1, r.stn3[s], msr, adj, r.cor[s], r.msr_sd[s], r.adj_sd[s], r.cor_sd[s],
                                 r.nstat[s], r.msr_id[s + 1], r.epoch, r.source)

                w.close()

            # H/I/J/P/Q/R/Y msrs
            elif l in one_stn_msrs:
                common_fn.write_prj(shp_name, adj_object.metadata.reference_frame)

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
                        w.point(adj_object.stns[r.stn1].lon.dec(), adj_object.stns[r.stn1].lat.dec())

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

                        w.point(adj_object.stns[r.stn1].lon.dec(), adj_object.stns[r.stn1].lat.dec())

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

                common_fn.write_prj(shp_name, adj_object.metadata.reference_frame)

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
                                [[adj_object.stns[r.stn1].lon.dec(), adj_object.stns[r.stn1].lat.dec()],
                                 [adj_object.stns[r.stn2].lon.dec(), adj_object.stns[r.stn2].lat.dec()]],
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
                                [[adj_object.stns[r.stn1].lon.dec(), adj_object.stns[r.stn1].lat.dec()],
                                 [adj_object.stns[r.stn2].lon.dec(), adj_object.stns[r.stn2].lat.dec()]],
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
                            [[adj_object.stns[r.stn1].lon.dec(), adj_object.stns[r.stn1].lat.dec()],
                             [adj_object.stns[r.stn2].lon.dec(), adj_object.stns[r.stn2].lat.dec()]],
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
                common_fn.write_prj(shp_name, adj_object.metadata.reference_frame)

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
                        [[adj_object.stns[r.stn1].lon.dec(), adj_object.stns[r.stn1].lat.dec()],
                         [adj_object.stns[r.stn2].lon.dec(), adj_object.stns[r.stn2].lat.dec()]],
                        [[adj_object.stns[r.stn2].lon.dec(), adj_object.stns[r.stn2].lat.dec()],
                         [adj_object.stns[r.stn3].lon.dec(), adj_object.stns[r.stn3].lat.dec()]],
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

    msr_type = line[:1]
    stn1 = line[2:22].strip()
    stn2 = line[22:42].strip()
    if stn2 == '':
        stn2 = None
    stn3 = line[42:62].strip()
    if stn3 == '':
        stn3 = None

    if msr_type == 'D':
        if msr_id_switch:
            items = line.split()
            cluster_id = items[-1].strip()
            msr_id = items[-2]
        else:
            msr_id = None
            cluster_id = None
        msr_components = {
            'msr_type': msr_type,
            'stn1': stn1,
            'stn2': stn2,
            'stn3': stn3,
            'msr_id': msr_id,
            'cluster_id': cluster_id
        }

        return msr_components
    else:
        if line[62:65].strip() == '*':
            ignore = True
        else:
            ignore = False

        cardinal = line[65:67].strip()
        msr = line[67:86].strip()
        adj = line[86:105].strip()
        cor = line[105:117].strip()
        msr_sd = line[117:130].strip()
        adj_sd = line[130:143].strip()
        cor_sd = line[143:156].strip()
        nstat = line[156:167].strip()

        if tstat_switch:
            tstat = line[167:178].strip()
            pelzer = line[178:190].strip()
            pre_adj_cor = line[190:204].strip()
            outlier = line[204:216]
            if msr_id_switch:
                msr_id = line[216:226].strip()
                cluster_id = line[226:236].strip()
            else:
                msr_id = None
                cluster_id = None
        else:
            tstat = None
            pelzer = line[167:178].strip()
            pre_adj_cor = line[180:193].strip()
            outlier = line[204:205]
            if msr_id_switch:
                msr_id = line[205:216].strip()
                cluster_id = line[216:226].strip()
            else:
                msr_id = None
                cluster_id = None

        msr_components = {
            'msr_type': msr_type,
            'stn1': stn1,
            'stn2': stn2,
            'stn3': stn3,
            'ignore': ignore,
            'cardinal': cardinal,
            'msr': msr,
            'adj': adj,
            'cor': cor,
            'msr_sd': msr_sd,
            'adj_sd': adj_sd,
            'cor_sd': cor_sd,
            'nstat': nstat,
            'tstat': tstat,
            'pelzer': pelzer,
            'pre_adj_cor': pre_adj_cor,
            'outlier': outlier,
            'msr_id': msr_id,
            'cluster_id': cluster_id
        }

        return msr_components


# ----------------------------------------------------------------------
# Sample usage
# ----------------------------------------------------------------------

adj_file = ''
xyz_file = ''

if xyz_file:
    dyna_stations = DynaXYZ(xyz_file)

if adj_file:
    # consume adj file
    if xyz_file:
        dyna_results = DynaAdj(adj_file, dyna_stations.stns, dyna_stations.metadata)
    else:
        dyna_results = DynaAdj(adj_file)

    # print file metadata
    out_str = 'DynAdjust Network Metadata\n'
    out_str += '-'*80 + '\n'
    out_str += f'File name:           {dyna_results.file_name}\n'
    out_str += f'File date:           {dyna_results.file_date}\n'
    out_str += f'Version:             {dyna_results.metadata.version}\n'
    out_str += f'Reference Frame:     {dyna_results.metadata.reference_frame}\n'
    out_str += f'Epoch:               {dyna_results.metadata.epoch}\n'
    out_str += f'Geoid Model:         {dyna_results.metadata.geoid_model}\n'
    out_str += f'Solution Type:       {dyna_results.adj_stats.soln_type}\n'
    out_str += f'Run time:            {dyna_results.adj_stats.run_time}\n'
    out_str += f'Measurements:        {dyna_results.adj_stats.msr_count}\n'
    out_str += f'Parameters:          {dyna_results.adj_stats.parameters}\n'
    out_str += f'Degrees of Freedom:  {dyna_results.adj_stats.degrees_of_freedom}\n'
    out_str += f'Global Pelzer:       {dyna_results.adj_stats.global_pelzer}\n'
    out_str += f'Potential Outliers:  {dyna_results.adj_stats.outlier_count}\n'
    out_str += f'Chi Squared (95%):   {dyna_results.adj_stats.chi_squared}\n'
    out_str += f'Chi Square window:   {dyna_results.adj_stats.chi_square_lower} - {dyna_results.adj_stats.chi_square_upper}\n'
    out_str += f'Sigma Zero:          {dyna_results.adj_stats.sigma_zero}\n'
    out_str += f'Chi Square Result:   {dyna_results.adj_stats.chi_square_result}\n'
    out_str += '-' * 80 + '\n'

    # print measurements
    if dyna_results.msrs:
        out_str += '\n\nMeasurements\n'
        out_str += 'M {:20s} {:20s} {:20s} {:>2s} {:>2s} {:>16s} {:>16s}\n'.format(
            'Station 1',
            'Station 2',
            'Station 3',
            '*',
            'C',
            'Measured',
            'Adjusted',
        )
        out_str += '-' * 120 + '\n'

        for m in dyna_results.msrs.values():
            if m.ignore:
                ignore = '*'
            else:
                ignore = ''

            # prepare msr/adj str
            if m.msr_type in angle_msrs:
                if m.msr_type != 'D':
                    msr_str = '  {:3d} {:2d} {:7.4f}'.format(m.msr.degree, m.msr.minute, m.msr.second)
                    adj_str = '  {:3d} {:2d} {:7.4f}'.format(m.adj.degree, m.adj.minute, m.adj.second)
            else:
                msr_str = str(m.msr)
                adj_str = str(m.adj)

            # direction sets
            if m.msr_type == 'D':
                out_str += '{:s} {:20s} {:20s}\n'.format(m.msr_type, m.stn1, m.stn2)
                for d in range(len(m.stn3)):
                    if m.ignore[d]:
                        ignore = '*'
                    else:
                        ignore = ''

                    out_str += '{:43s} {:20s} {:>2s} {:>2s}   {:3d} {:2d} {:7.4f}   {:3d} {:2d} {:7.4f}\n'.format(
                        '',
                        m.stn3[d],
                        ignore,
                        '',
                        m.msr[d].degree,
                        m.msr[d].minute,
                        m.msr[d].second,
                        m.adj[d].degree,
                        m.adj[d].minute,
                        m.adj[d].second
                    )

            # three line msrs (G/X/Y)
            elif m.msr_type in three_line_msrs:
                # G/X
                if m.stn2:
                    out_str += '{:s} {:20s} {:20s} {:20s} {:>2s} {:>2s} {:16.4f} {:16.4f}\n'.format(
                        m.msr_type,
                        m.stn1,
                        m.stn2,
                        '',
                        ignore,
                        m.cardinal[0],
                        m.msr[0],
                        m.adj[0]
                    )
                    for i in range(1, 3):
                        out_str += '{:64s} {:>2s} {:>2s} {:16.4f} {:16.4f}\n'.format(
                            '',
                            ignore,
                            m.cardinal[i],
                            m.msr[i],
                            m.adj[i]
                        )
                # Y
                else:
                    out_str += '{:s} {:20s} {:20s} {:20s} {:>2s} {:>2s} {:16.4f} {:16.4f}\n'.format(
                        m.msr_type,
                        m.stn1,
                        '',
                        '',
                        ignore,
                        m.cardinal[0],
                        m.msr[0],
                        m.adj[0]
                    )
                    for i in range(1, 3):
                        out_str += '{:64s} {:>2s} {:>2s} {:16.4f} {:16.4f}\n'.format(
                            '',
                            ignore,
                            m.cardinal[i],
                            m.msr[i],
                            m.adj[i]
                        )

            # single-line three-stn msrs
            elif m.stn3:
                out_str += '{:s} {:20s} {:20s} {:20s} {:>2s} {:>2s} {:16s} {:16s}\n'.format(
                    m.msr_type,
                    m.stn1,
                    m.stn2,
                    m.stn3,
                    ignore,
                    m.cardinal,
                    msr_str,
                    adj_str
                )

            # single-line two-stn msrs
            elif m.stn2:
                out_str += '{:s} {:20s} {:20s} {:20s} {:>2s} {:>2s} {:>16s} {:>16s}\n'.format(
                    m.msr_type,
                    m.stn1,
                    m.stn2,
                    '',
                    ignore,
                    m.cardinal,
                    msr_str,
                    adj_str
                )

            # single-line one-stn msrs
            else:
                out_str += '{:s} {:20s} {:20s} {:20s} {:>2s} {:>2s} {:>16s} {:>16s}\n'.format(
                    m.msr_type,
                    m.stn1,
                    '',
                    '',
                    ignore,
                    m.cardinal,
                    msr_str,
                    adj_str
                )

    # print adjusted coordinates
    if dyna_results.stns:
        out_str += '\n\nAdjusted Coordinates\n'
        out_str += '{:20s} {:>15s} {:>15s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s}\n'.format(
            'Station',
            'Latitude',
            'Longitude',
            'Ell_Ht',
            'Phys_Ht',
            'SD_East',
            'SD_North',
            'SD_Up',
        )

        out_str += '-'*120 + '\n'

        for m in dyna_results.stns.values():
            out_str += '{:20s} {:15.8f} {:15.8f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f}\n'.format(
                m.name,
                m.lat.dec(),
                m.lon.dec(),
                m.ehgt,
                m.ohgt,
                m.sd_e,
                m.sd_n,
                m.sd_u,
            )

        out_str += '\n'

    with open('pynadjust_adj.txt', 'w') as out_fh:
        out_fh.write(out_str)

    write_adj_shapefile(dyna_results, 'test_network')
