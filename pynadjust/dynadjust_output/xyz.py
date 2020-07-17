# pynadjust xyz module - retrieve data from .xyz files
import geodepy.convert as gc
import datetime


class StationLibrary(object):
    def __init__(self):
        self.stations = {}
        self.epoch = None
        self.reference_frame = None
        self.geoid_model = None
        self.version = None
        self.file_name = None
        self.file_date = None

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
            return gc.llh2xyz(self.lat, self.lon, self.ehgt)

        def grid(self):
            return gc.geo2grid(self.lat, self.lon)


def read_xyz(xyz_file):

    dyna_xyz = StationLibrary()

    with open(xyz_file, 'r') as xyz_fh:
        line_count = 0
        stn_line = False
        stn_switch = False
        mandatory_coord_types = 'PLHh'

        for line in xyz_fh:
            line_count += 1

            if 'Adjusted Coordinates' in line:
                stn_line = line_count + 5

            if stn_line:
                if line_count == stn_line:
                    stn_switch = True

            if line[:35] == 'Version:                           ':
                dyna_xyz.version = line[35:].strip()

            if line[:35] == 'Reference frame:                   ':
                dyna_xyz.reference_frame = line[35:].strip()

            if line[:35] == 'File name:                         ':
                dyna_xyz.file_name = line[35:].strip()

            if line[:35] == 'File created:                      ':
                date_str = line[35:].strip()
                dyna_xyz.file_date = datetime.datetime.strptime(date_str, '%A, %d %B %Y, %I:%M:%S %p')
                dyna_xyz.file_name = line[35:].strip()

            if line[:35] == 'Epoch:                             ':
                date_str = line[35:].strip()
                d = int(date_str[:2])
                m = int(date_str[3:5])
                y = int(date_str[6:])
                dyna_xyz.epoch = datetime.date(y, m, d)

            if line[:35] == 'Geoid model:                       ':
                dyna_xyz.geoid_model = line[35:].strip()

            if line[:35] == 'Station coordinate types:          ':
                coord_types = line[35:].strip()
                missing_coord_types = set()
                for l in mandatory_coord_types:
                    if l not in coord_types:
                        missing_coord_types.add(l)
                if missing_coord_types:
                    raise ValueError(f'Mandatory coordinate types {missing_coord_types} not present in {xyz_file}')

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
                        P = gc.hp2dec(float(results[r_count]))
                    if ct == 'L':
                        L = gc.hp2dec(float(results[r_count]))
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
                dyna_xyz.stations[stn] = dyna_xyz.Station()
                dyna_xyz.stations[stn].name = stn
                dyna_xyz.stations[stn].con = con
                dyna_xyz.stations[stn].lat = P
                dyna_xyz.stations[stn].lon = L
                dyna_xyz.stations[stn].ehgt = h
                dyna_xyz.stations[stn].ohgt = H
                dyna_xyz.stations[stn].sd_e = sd_e
                dyna_xyz.stations[stn].sd_n = sd_n
                dyna_xyz.stations[stn].sd_u = sd_u

    return dyna_xyz


# ----------------------------------------------------------------------
# Example of usage
# ----------------------------------------------------------------------

xyz_file = ''

if xyz_file:
    xyz_results = read_xyz(xyz_file)

    # iterate through and print coordinates
    out_str = ''
    # print metadata
    out_str += f'File name:       {xyz_results.file_name}\n'
    out_str += f'Version:         {xyz_results.version}\n'
    out_str += f'Date created:    {xyz_results.file_date}\n'
    out_str += f'Reference Frame: {xyz_results.reference_frame}\n'
    out_str += f'Epoch:           {xyz_results.epoch}\n'
    out_str += f'Geoid model:     {xyz_results.geoid_model}\n'

    # print results
    out_str += '\nStations:\n'
    out_str += '-' * 178 + '\n'
    for stn in xyz_results.stations:
        out_str += '{:20s} {:3s} {:14.8f} {:14.8f} {:10.4f} {:10.4f} {:8.4f} {:8.4f} {:8.4f} '.format(
            xyz_results.stations[stn].name,
            xyz_results.stations[stn].con,
            xyz_results.stations[stn].lat,
            xyz_results.stations[stn].lon,
            xyz_results.stations[stn].ehgt,
            xyz_results.stations[stn].ohgt,
            xyz_results.stations[stn].sd_e,
            xyz_results.stations[stn].sd_n,
            xyz_results.stations[stn].sd_u
        )

        grid = xyz_results.stations[stn].grid()
        xyz = xyz_results.stations[stn].xyz()

        out_str += '{:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f}\n'.format(
            grid[2],
            grid[3],
            xyz[0],
            xyz[1],
            xyz[2]
        )

    with open('scratch.txt', 'w') as scratch_fh:
        scratch_fh.write(out_str)
