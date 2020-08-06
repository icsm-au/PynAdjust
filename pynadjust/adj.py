# pynadjust adj module - retrieve data from .adj files

# For review/consideration:
#     - Confirm all measurement types and components read and stored correctly (test on adj will all msr types)
#     - Store Latitude/Longitude as DMSAngle objects (or as Decimal degrees float)
#     - Consider storing measurement objects in a dictionary (done here), or set.
#     - expand shapefile function to include all fields
#     - expand shapefile function to include option for uncertainty figures using apu object (once apu.py complete)


import datetime
import geodepy.convert as gc
import shapefile
import os


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
    def __init__(self):
        self.stns = {}
        self.msrs = {}
        self.epoch = None
        self.reference_frame = None
        self.geoid_model = None
        self.version = None
        self.file_name = None
        self.file_date = None
        self.soln_type = None
        self.run_time = None
        self.parameters = None
        self.msr_count = None
        self.degrees_of_freedom = None
        self.sigma_zero = None
        self.chi_squared = None
        self.outlier_count = None
        self.global_pelzer = None
        self.chi_square_lower = None
        self.chi_square_upper = None
        self.chi_square_result = None

    class Station(object):
        def __init__(self):
            self.name = None
            self.con = None
            self.lat = None
            self.lon = None
            self.ohgt = None
            self.ehgt = None
            self.sd_e = None
            self.sd_n = None
            self.sd_u = None

        def xyz(self):
            return gc.llh2xyz(self.lat.dec(), self.lon.dec(), self.ehgt)

        def grid(self):
            return gc.geo2grid(self.lat.dec(), self.lon.dec())

    class Measurement(object):
        def __init__(self):
            self.msr_type = None
            self.stn1 = None
            self.stn2 = None
            self.stn3 = None
            self.ignore = None
            self.cardinal = None
            self.msr = None
            self.adj = None
            self.cor = None
            self.msr_sd = None
            self.adj_sd = None
            self.cor_sd = None
            self.nstat = None
            self.tstat = None
            self.pelzer = None
            self.pre_adj_cor = None
            self.outlier = None
            self.msr_id = None
            self.cluster_id = None


def write_shapefile(dyna_adj, network_name):

    def write_prj(file_name, ref_frame):

        if ref_frame == 'GDA2020':
            out_str = 'GEOGCS["GDA2020",DATUM["GDA2020",SPHEROID["GRS_1980",6378137.0,298.257222101]],' \
                      'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433],AUTHORITY["EPSG",7844]]'
        elif ref_frame == 'GDA94':
            out_str = 'GEOGCS["GDA94",DATUM["D_GDA_1994",SPHEROID["GRS_1980",6378137,298.257222101]],' \
                      'PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]'
        else:
            out_str = None

        if out_str != None:
            prj_fh = open(file_name + '.prj', 'w')
            prj_fh.write(out_str)
            prj_fh.close()

    def find_cardinal(msrs_object):
        for m in msrs_object:
            if m.cardinal == ['e', 'n', 'u']:
                cardinal = 'enu'
            else:
                cardinal = 'XYZ'
            break

        return cardinal

    def max_stat(stat_list):
        max_stat = None
        for s in stat_list:
            if not max_stat:
                max_stat = s
            elif abs(s) > abs(max_stat):
                max_stat = s

        return max_stat

    if os.path.exists('shp'):
        os.chdir('shp')
    else:
        os.mkdir('shp')
        os.chdir('shp')

    if dyna_adj.stns:
        shp_name = network_name + '_stn'

        write_prj(shp_name, dyna_adj.reference_frame)

        w = shapefile.Writer(shp_name, shapeType=1)

        w.autoBalance = 1

        w.field('Station', 'C', size=20)
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

        for s in dyna_adj.stns.values():
            w.point(s.lon.dec(), s.lat.dec())

            grid = s.grid()

            w.record(
                s.name, s.con,
                grid[2], grid[3], grid[1],
                s.lat.dec(), s.lon.dec(),
                s.ohgt, s.ehgt,
                s.sd_e, s.sd_n, s.sd_u
            )

        w.close()

    for l in msr_types:
        result = [dyna_adj.msrs[m] for m in dyna_adj.msrs if dyna_adj.msrs[m].msr_type == l]

        shp_name = network_name + '_' + l

        if result:
            # Direction sets
            if l in multi_line_msrs:
                write_prj(shp_name, dyna_adj.reference_frame)

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

                for r in result:
                    w.line([
                        [[dyna_adj.stns[r.stn1].lon.dec(), dyna_adj.stns[r.stn1].lat.dec()],
                         [dyna_adj.stns[r.stn2].lon.dec(), dyna_adj.stns[r.stn2].lat.dec()]],
                    ])

                    w.record(r.stn1, r.stn2, 'set zero', None, None, None, None, None, None)

                    for s in range(len(r.stn3)):
                        msr = '{:3d} {:2d} {:7.4f}'.format(r.msr[s].degree, r.msr[s].minute, r.msr[s].second)
                        adj = '{:3d} {:2d} {:7.4f}'.format(r.adj[s].degree, r.adj[s].minute, r.adj[s].second)

                        w.line([
                            [[dyna_adj.stns[r.stn1].lon.dec(), dyna_adj.stns[r.stn1].lat.dec()],
                             [dyna_adj.stns[r.stn3[s]].lon.dec(), dyna_adj.stns[r.stn3[s]].lat.dec()]],
                        ])

                        w.record(r.stn1, r.stn3[s], msr, adj, r.cor[s], r.msr_sd[s], r.adj_sd[s], r.cor_sd[s], r.nstat[s])

                w.close()

            # H/I/J/P/Q/R/Y msrs
            elif l in one_stn_msrs:
                write_prj(shp_name, dyna_adj.reference_frame)

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

                    for r in result:
                        w.point(dyna_adj.stns[r.stn1].lon.dec(), dyna_adj.stns[r.stn1].lat.dec())

                        w.record(r.stn1, r.msr, r.adj, r.cor, r.msr_sd, r.adj_sd, r.cor_sd, r.nstat)

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

                    for r in result:
                        max_nstat = max_stat(r.nstat)
                        max_cor = max_stat(r.cor)

                        w.point(dyna_adj.stns[r.stn1].lon.dec(), dyna_adj.stns[r.stn1].lat.dec())

                        w.record(r.stn1,
                                 r.msr[0], r.msr[1], r.msr[2],
                                 r.adj[0], r.adj[1], r.adj[2],
                                 r.cor[0], r.cor[1], r.cor[2],
                                 r.msr_sd[0], r.msr_sd[1], r.msr_sd[2],
                                 r.adj_sd[0], r.adj_sd[1], r.adj_sd[2],
                                 r.cor_sd[0], r.cor_sd[1], r.cor_sd[2],
                                 r.nstat[0], r.nstat[1], r.nstat[2],
                                 max_cor, max_nstat
                                 )

                w.close()

            # B/C/E/G/K/L/M/S/V/X/Z
            elif l in two_stn_msrs:

                write_prj(shp_name, dyna_adj.reference_frame)

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

                        for r in result:
                            w.line([
                                [[dyna_adj.stns[r.stn1].lon.dec(), dyna_adj.stns[r.stn1].lat.dec()],
                                 [dyna_adj.stns[r.stn2].lon.dec(), dyna_adj.stns[r.stn2].lat.dec()]],
                            ])

                            msr = '{:3d} {:2d} {:7.4f}'.format(r.msr.degree, r.msr.minute, r.msr.second)
                            adj = '{:3d} {:2d} {:7.4f}'.format(r.adj.degree, r.adj.minute, r.adj.second)

                            w.record(r.stn1, r.stn2, msr, adj, r.cor, r.msr_sd, r.adj_sd, r.cor_sd, r.nstat)

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

                        for r in result:
                            w.line([
                                [[dyna_adj.stns[r.stn1].lon.dec(), dyna_adj.stns[r.stn1].lat.dec()],
                                [dyna_adj.stns[r.stn2].lon.dec(), dyna_adj.stns[r.stn2].lat.dec()]],
                            ])

                            w.record(r.stn1, r.stn2, r.msr, r.adj, r.cor, r.msr_sd, r.adj_sd, r.cor_sd, r.nstat)

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

                    for r in result:
                        max_nstat = max_stat(r.nstat)
                        max_cor = max_stat(r.cor)

                        w.line([
                            [[dyna_adj.stns[r.stn1].lon.dec(), dyna_adj.stns[r.stn1].lat.dec()],
                            [dyna_adj.stns[r.stn2].lon.dec(), dyna_adj.stns[r.stn2].lat.dec()]],
                        ])

                        w.record(r.stn1, r.stn2,
                                 r.msr[0], r.msr[1], r.msr[2],
                                 r.adj[0], r.adj[1], r.adj[2],
                                 r.cor[0], r.cor[1], r.cor[2],
                                 r.msr_sd[0], r.msr_sd[1], r.msr_sd[2],
                                 r.adj_sd[0], r.adj_sd[1], r.adj_sd[2],
                                 r.cor_sd[0], r.cor_sd[1], r.cor_sd[2],
                                 r.nstat[0], r.nstat[1], r.nstat[2],
                                 max_cor, max_nstat
                                 )

                w.close()

            # A (not D)
            elif l in three_stn_msrs:
                write_prj(shp_name, dyna_adj.reference_frame)

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

                for r in result:
                    msr = '{:3d} {:2d} {:7.4f}'.format(r.msr.degree, r.msr.minute, r.msr.second)
                    adj = '{:3d} {:2d} {:7.4f}'.format(r.adj.degree, r.adj.minute, r.adj.second)

                    w.line([
                        [[dyna_adj.stns[r.stn1].lon.dec(), dyna_adj.stns[r.stn1].lat.dec()],
                         [dyna_adj.stns[r.stn2].lon.dec(), dyna_adj.stns[r.stn2].lat.dec()]],
                        [[dyna_adj.stns[r.stn2].lon.dec(), dyna_adj.stns[r.stn2].lat.dec()],
                         [dyna_adj.stns[r.stn3].lon.dec(), dyna_adj.stns[r.stn3].lat.dec()]],
                    ])

                    w.record(r.stn1, r.stn2, r.stn3, r.msr, r.adj, r.cor, r.msr_sd, r.adj_sd, r.cor_sd, r.nstat)

                w.close()


def str_to_dms_angle(angle):
    items = angle.split()
    deg = int(items[0])
    min = int(items[1])
    sec = float(items[2])
    dms = gc.DMSAngle(degree=deg, minute=min, second=sec)

    return dms


def read_msr_fields(line):
    t_stat_field = False
    msr_id_field = False

    if 'T-stat' in line:
        t_stat_field = True

    if 'Msr ID' in line:
        msr_id_field = False

    return t_stat_field, msr_id_field


def read_msr_line(line, tstat_switch, msr_id_switch):
    msr_type = line[:1]
    stn1 = line[2:22].strip()
    stn2 = line[22:42].strip()
    if stn2 == '':
        stn2 = None
    stn3 = line[42:62].strip()
    if stn3 == '':
        stn3 = None
    if msr_type == 'D':
        msr_components = {
            'msr_type': msr_type,
            'stn1': stn1,
            'stn2': stn2,
            'stn3': stn3
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


def read_adj(adj_file):

    dyna_adj = DynaAdj()

    with open(adj_file, 'r') as adj_fh:
        line_count = 0
        stn_line = None
        stn_switch = False
        mandatory_coord_types = 'PLHh'

        msr_line = None
        msr_switch = False
        msr_count = 0
        dir_set = 0
        msr_header_line = None

        metadata_switch = True

        for line in adj_fh:
            line_count += 1

            if 'Adjusted Measurements' in line:
                msr_line = line_count + 5
                msr_header_line = line_count + 3
                three_line_count = 0

            if msr_line:
                if line_count == msr_line:
                    metadata_switch = False
                    msr_switch = True

            if msr_header_line:
                if line_count == msr_header_line:
                    tstat_switch, msr_id_switch = read_msr_fields(line)

            if 'Adjusted Coordinates' in line:
                stn_line = line_count + 5

            if stn_line:
                if line_count == stn_line:
                    stn_switch = True

            if metadata_switch:
                if line[:35] == 'Version:                           ':
                    dyna_adj.version = line[35:].strip()

                if line[:35] == 'Reference frame:                   ':
                    dyna_adj.reference_frame = line[35:].strip()

                if line[:35] == 'File name:                         ':
                    dyna_adj.file_name = line[35:].strip()

                if line[:35] == 'File created:                      ':
                    date_str = line[35:].strip()
                    dyna_adj.file_date = datetime.datetime.strptime(date_str, '%A, %d %B %Y, %I:%M:%S %p')
                    dyna_adj.file_name = line[35:].strip()

                if line[:35] == 'Epoch:                             ':
                    date_str = line[35:].strip()
                    d = int(date_str[:2])
                    m = int(date_str[3:5])
                    y = int(date_str[6:])
                    dyna_adj.epoch = datetime.date(y, m, d)

                if line[:35] == 'Geoid model:                       ':
                    dyna_adj.geoid_model = line[35:].strip()

                if line[:35] == 'SOLUTION                           ':
                    dyna_adj.soln_type = line[35:].strip()

                if line[:35] == 'Total time                         ':
                    items = line[35:].strip().split(':')
                    hour = int(items[0])
                    minute = int(items[1])
                    second = float(items[2])
                    dyna_adj.run_time = datetime.timedelta(hours=hour, minutes=minute, seconds=second)

                if line[:35] == 'Number of unknown parameters       ':
                    dyna_adj.parameters = int(line[35:])

                if line[:35] == 'Number of measurements             ':
                    items = line[35:].split()
                    if len(items) > 1:
                        if '(' in items[1]:
                            outliers = int(items[1].replace('(', ''))
                    dyna_adj.msr_count = int(items[0])
                    dyna_adj.outlier_count = outliers

                if line[:35] == 'Degrees of freedom                 ':
                    dyna_adj.degrees_of_freedom = int(line[35:])

                if line[:35] == 'Chi squared                        ':
                    dyna_adj.chi_squared = float(line[35:])

                if line[:35] == 'Rigorous Sigma Zero                ':
                    dyna_adj.sigma_zero = float(line[35:])

                if line[:35] == 'Global (Pelzer) Reliability        ':
                    items = line[35:].strip().split()
                    dyna_adj.global_pelzer = float(items[0])

                if line[:35] == 'Chi-Square test (95.0%)            ':
                    items = line[35:].strip().split()
                    dyna_adj.chi_square_lower = float(items[0])
                    dyna_adj.chi_square_upper = float(items[4])
                    dyna_adj.chi_square_result = items[6].strip()

                if line[:35] == 'Station coordinate types:          ':
                    coord_types = line[35:].strip()
                    missing_coord_types = set()
                    for l in mandatory_coord_types:
                        if l not in coord_types:
                            missing_coord_types.add(l)
                    if missing_coord_types:
                        raise ValueError(f'Mandatory coordinate types {missing_coord_types} not present in {adj_file}')

            if stn_switch:
                if len(line) < 20:
                    stn_switch = False
                    continue

                P = None
                L = None
                H = None
                h = None

                stn = line[0:20].strip()
                con = line[20:23]
                results = line[25:].split()
                r_count = 0

                for ct in coord_types:
                    if ct == 'P':
                        P = gc.dec2dms(gc.hp2dec(float(results[r_count])))
                    if ct == 'L':
                        L = gc.dec2dms(gc.hp2dec(float(results[r_count])))
                    if ct == 'H':
                        H = float(results[r_count])
                    if ct == 'h':
                        h = float(results[r_count])

                    r_count += 1

                # Don't forget about the qualities
                sd_e = float(results[r_count])
                sd_n = float(results[r_count + 1])
                sd_u = float(results[r_count + 2])

                # initialise a Station object in dictionary, then populate fields
                dyna_adj.stns[stn] = dyna_adj.Station()
                dyna_adj.stns[stn].name = stn
                dyna_adj.stns[stn].con = con
                dyna_adj.stns[stn].lat = P
                dyna_adj.stns[stn].lon = L
                dyna_adj.stns[stn].ehgt = h
                dyna_adj.stns[stn].ohgt = H
                dyna_adj.stns[stn].sd_e = sd_e
                dyna_adj.stns[stn].sd_n = sd_n
                dyna_adj.stns[stn].sd_u = sd_u

            if msr_switch:
                if line.strip() == '':
                    msr_switch = False
                    continue

                msr_components = read_msr_line(line, tstat_switch, msr_id_switch)
                msr_type = msr_components['msr_type']

                # Direction sets - initialise new msr object
                if msr_type == 'D':
                    msr_count += 1
                    dir_set += 1
                    d_stn1 = msr_components['stn1']
                    d_stn2 = msr_components['stn2']

                    msr = dyna_adj.Measurement()
                    msr.msr_type = msr_type
                    msr.stn1 = d_stn1
                    msr.stn2 = d_stn2
                    msr.stn3 = []
                    msr.ignore = []
                    msr.msr = []
                    msr.adj = []
                    msr.cor = []
                    msr.msr_sd = []
                    msr.adj_sd = []
                    msr.cor_sd = []
                    msr.nstat = []
                    msr.tstat = []
                    msr.pelzer = []
                    msr.pre_adj_cor = []
                    msr.outlier = []
                    msr.msr_id = []
                    msr.cluster_id = []

                    dyna_adj.msrs[msr_count] = msr

                # Direction Set pointings - append msrs to lists
                elif msr_type == ' ':
                    dyna_adj.msrs[msr_count].stn3.append(msr_components['stn3'])
                    dyna_adj.msrs[msr_count].ignore.append(msr_components['ignore'])
                    dyna_adj.msrs[msr_count].msr.append(str_to_dms_angle(msr_components['msr']))
                    dyna_adj.msrs[msr_count].adj.append(str_to_dms_angle(msr_components['adj']))
                    dyna_adj.msrs[msr_count].cor.append(float(msr_components['cor']))
                    dyna_adj.msrs[msr_count].msr_sd.append(float(msr_components['msr_sd']))
                    dyna_adj.msrs[msr_count].adj_sd.append(float(msr_components['adj_sd']))
                    dyna_adj.msrs[msr_count].cor_sd.append(float(msr_components['cor_sd']))
                    dyna_adj.msrs[msr_count].nstat.append(float(msr_components['nstat']))
                    if tstat_switch:
                        dyna_adj.msrs[msr_count].tstat.append(float(msr_components['tstat']))
                    else:
                        dyna_adj.msrs[msr_count].tstat.append(None)
                    dyna_adj.msrs[msr_count].pelzer.append(float(msr_components['pelzer']))
                    dyna_adj.msrs[msr_count].pre_adj_cor.append(float(msr_components['pre_adj_cor']))
                    dyna_adj.msrs[msr_count].outlier.append(msr_components['outlier'])
                    dyna_adj.msrs[msr_count].msr_id.append(msr_components['msr_id'])
                    dyna_adj.msrs[msr_count].cluster_id.append(msr_components['cluster_id'])

                # Type X/Y/Z measurements (split over 3 lines)
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

                        msr = dyna_adj.Measurement()
                        msr.msr_type = msr_type
                        msr.stn1 = stn1
                        msr.stn2 = stn2
                        msr.stn3 = None
                        msr.ignore = ignore
                        msr.cardinal = cardinal_list
                        msr.msr = msr_list
                        msr.adj = adj_list
                        msr.cor = cor_list
                        msr.msr_sd = msr_sd_list
                        msr.adj_sd = adj_sd_list
                        msr.cor_sd = cor_sd_list
                        msr.nstat = nstat_list
                        msr.tstat = tstat_list
                        msr.pelzer = pelzer_list
                        msr.pre_adj_cor = pre_adj_cor_list
                        msr.outlier = outlier_list
                        msr.msr_id = msr_id
                        msr.cluster_id = cluster_id

                        dyna_adj.msrs[msr_count] = msr
                    else:
                        pass

                # all other msr types
                else:
                    msr_count += 1
                    msr = dyna_adj.Measurement()
                    msr.msr_type = msr_type
                    msr.stn1 = msr_components['stn1']
                    msr.stn2 = msr_components['stn2']
                    msr.stn3 = msr_components['stn3']
                    msr.ignore = msr_components['ignore']
                    msr.cardinal = msr_components['cardinal']
                    if msr_type in angle_msrs:
                        msr.msr = str_to_dms_angle(msr_components['msr'])
                        msr.adj = str_to_dms_angle(msr_components['adj'])
                    elif msr_type in hp_msrs:
                        msr.msr = gc.hp2dec(msr_components['msr'])
                        msr.adj = gc.hp2dec(msr_components['adj'])
                    else:
                        msr.msr = float(msr_components['msr'])
                        msr.adj = float(msr_components['adj'])
                    msr.cor = float(msr_components['cor'])
                    msr.msr_sd = float(msr_components['msr_sd'])
                    msr.adj_sd = float(msr_components['adj_sd'])
                    msr.cor_sd = float(msr_components['cor_sd'])
                    msr.nstat = float(msr_components['nstat'])
                    if tstat_switch:
                        msr.tstat = float(msr_components['tstat'])
                    msr.pelzer = float(msr_components['pelzer'])
                    msr.pre_adj_cor = float(msr_components['pre_adj_cor'])
                    msr.outlier = msr_components['outlier']
                    if msr_id:
                        msr.msr_id = float(msr_components['msr_id'])
                        msr.cluster_id = float(msr_components['cluster_id'])

                    dyna_adj.msrs[msr_count] = msr

    return dyna_adj


# ----------------------------------------------------------------------
# Sample usage
# ----------------------------------------------------------------------

adj_file = ''

if adj_file:
    # consume adj file
    dyna_results = read_adj(adj_file)

    # print file metadata
    out_str = 'DynAdjust Network Metadata\n'
    out_str += '-'*80 + '\n'
    out_str += f'File name:           {dyna_results.file_name}\n'
    out_str += f'File date:           {dyna_results.file_date}\n'
    out_str += f'Version:             {dyna_results.version}\n'
    out_str += f'Reference Frame:     {dyna_results.reference_frame}\n'
    out_str += f'Epoch:               {dyna_results.epoch}\n'
    out_str += f'Geoid Model:         {dyna_results.geoid_model}\n'
    out_str += f'Solution Type:       {dyna_results.soln_type}\n'
    out_str += f'Run time:            {dyna_results.run_time}\n'
    out_str += f'Measurements:        {dyna_results.msr_count}\n'
    out_str += f'Parameters:          {dyna_results.parameters}\n'
    out_str += f'Degrees of Freedom:  {dyna_results.degrees_of_freedom}\n'
    out_str += f'Global Pelzer:       {dyna_results.global_pelzer}\n'
    out_str += f'Potential Outliers:  {dyna_results.outlier_count}\n'
    out_str += f'Chi Squared (95%):   {dyna_results.chi_squared}\n'
    out_str += f'Chi Square window:   {dyna_results.chi_square_lower} - {dyna_results.chi_square_upper}\n'
    out_str += f'Sigma Zero:          {dyna_results.sigma_zero}\n'
    out_str += f'Chi Square Result:   {dyna_results.chi_square_result}\n'
    out_str += '-' * 80 + '\n'

    # print measurements
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

    with open('adj_read.txt', 'w') as out_fh:
        out_fh.write(out_str)

    write_shapefile(dyna_results, 'test_network')
