# ----------------------------------------------------------------------
#                          DynaDiff.py
# ----------------------------------------------------------------------
# Author:   Nicholas Gowans
# Date:     17 May 2019
# Purpose:  To compare difference between common stations in two
#           DynAdjust *xyz files.
#               
# ----------------------------------------------------------------------          
# Usage:    Cmd:\>python DynaDiff.py <xyz_file_1> <xyz_file_2>
# ----------------------------------------------------------------------
# Notes:    Updated 30 SEP 2019:
#               - Reads coordinate info based on DynAdjust coordinate type option
#               - Print error message if .xyz files are missing essential coord types
#               - Write point differences shapefile
#
#           Updated 02 OCT 2019:
#               - Altered read_switch trigger for different coord_type formats
#               - Incorporated read_xyz function from Josh Batchelor
# ----------------------------------------------------------------------

import os
import math
import geodepy.convert as gc
import shapefile

dateUpdated = '20190210'
mandatory_coord_types = 'ENzPLHh'

xyz1_file = os.sys.argv[1]
xyz2_file = os.sys.argv[2]
output_rpt_file = '{:s}-{:s}_DynaDiff.rpt'.format(xyz2_file, xyz1_file)
output_csv_file = '{:s}-{:s}_DynaDiff.csv'.format(xyz2_file, xyz1_file)


def shortname(filename):
    """returns the string truncated at the first period"""
    count = 0
    len_fn = len(filename)

    for l in filename:
        count += 1
        if l == '.':
            short = filename[:count - 1]
            break
        # no period in file name
        if count == len_fn:
            short = filename

    return short


# pyshp truncates shapefile names at the first period.
# Provide a meaningful name, with no periods.
short1 = shortname(xyz1_file)
short2 = shortname(xyz2_file)
shp_file = '{:s}-{:s}_p'.format(short1, short2)

# open output files, and print header information
rep_fh = open(output_rpt_file, 'w')
csv_fh = open(output_csv_file, 'w')

headerStr = '-'*80 + '\n'
headerStr += 'DynaDiff.py Report File\n'
headerStr += '-'*30 + '\n'
headerStr += 'Coordinates differences computed by File 2 minus File 1\n'
headerStr += 'File 1:  ' + xyz1_file + '\n'
headerStr += 'File 2:  ' + xyz2_file + '\n'
headerStr += 'Version: ' + dateUpdated + '\n'
headerStr += '-'*80 + '\n'
print(headerStr, file=rep_fh)

csvHeaderStr = 'Station,Latitude,Longitude,dE,dN,Dist,Brg,dZ,dOHGT,dEHGT,dSD_E,dSD_N,dSD_U'
print(csvHeaderStr, file=csv_fh)

# initialisations
coords1 = {}
coords2 = {}


# ----------------------------------------------------------------------
# Read in First and Second xyz files
# ----------------------------------------------------------------------


def read_xyz(xyz_file, read_switch=False):
    coords = {}

    print()
    print(' Reading {:s}'.format(xyz_file))

    with open(xyz_file, 'r') as xyz_fh:

        for line in xyz_fh:
            # skip empty lines
            if line == '':
                continue
            if line == '\n':
                continue

            if read_switch is False:
                if line[:35] == 'Version:                           ':
                    vers = line[35:].strip()
                if line[:35] == 'File name:                         ':
                    file = line[35:].strip()
                if line[:35] == 'Reference frame:                   ':
                    ref_frame = line[35:].strip()
                if line[:35] == 'Epoch:                             ':
                    epoch = line[35:].strip()
                if line[:35] == 'Geoid model:                       ':
                    geoid = line[35:].strip()
                if line[:35] == 'Station coordinate types:          ':
                    coord_types = line[35:].strip()
                    for l in mandatory_coord_types:
                        if l not in coord_types:
                            print('*' * 20)
                            print(' Warning! Mandatory coordinate types not present in {:s}'.format(xyz_file))
                            print(' .xyz file must contain coord types ENzPLHh')
                            print('')
                            print('Exiting...')
                            exit()

            if line[87:88] == '-':
                read_switch = True
                continue

            if read_switch:
                stn = line[0:20].strip()
                results = line[25:].split()
                r_count = 0

                for ct in coord_types:
                    if ct == 'E':
                        E = float(results[r_count])
                    if ct == 'N':
                        N = float(results[r_count])
                    if ct == 'z':
                        z = int(results[r_count])
                    if ct == 'P':
                        P = gc.hp2dec(float(results[r_count]))
                    if ct == 'L':
                        L = gc.hp2dec(float(results[r_count]))
                    if ct == 'H':
                        H = float(results[r_count])
                    if ct == 'h':
                        h = float(results[r_count])
                    # Cartesian coords disabled for now
                    # if ct == 'X':
                    #     X = float(results[r_count])
                    # if ct == 'Y':
                    #     Y = float(results[r_count])
                    # if ct == 'Z':
                    #     Z = float(results[r_count])

                    r_count += 1

                # Don't forget about the qualities
                SE = float(results[r_count])
                SN = float(results[r_count + 1])
                SU = float(results[r_count + 2])

                coords[stn] = {'E': E,
                               'N': N,
                               'Z': z,
                               'LAT': P,
                               'LON': L,
                               'OHGT': H,
                               'EHGT': h,
                               'SE': SE,
                               'SN': SN,
                               'SU': SU
                               }

    return vers, file, ref_frame, epoch, geoid, coord_types, coords


vers1, file1, ref_frame1, epoch1, geoid1, coord_types1, coords1 = read_xyz(xyz1_file)
vers2, file2, ref_frame2, epoch2, geoid2, coord_types2, coords2 = read_xyz(xyz2_file)

print()
print(' Read:')
print('   - {:8d} stations in {:s}'.format(len(coords1), xyz1_file))
print('   - {:8d} stations in {:s}'.format(len(coords2), xyz2_file))


# ----------------------------------------------------------------------
# Print warnings
# ----------------------------------------------------------------------

warn_str = ''

if vers1 != vers2:
    warn_str = warn_str + ' * Different version of DynAdjust used between adjustments.\n'
if ref_frame1 != ref_frame2:
    warn_str = warn_str + ' * Reference frame differs between adjustments.\n'
if epoch1 != epoch2:
    warn_str = warn_str + ' * Epoch differs between adjustments.\n'
if geoid1 != geoid2:
    warn_str = warn_str + ' * Geoid model differs between adjustments.\n'
if coord_types1 != coord_types2:
    warn_str = warn_str + ' * Different coordinate types printed between adjustments.\n'

print()
print('-'*20, file=rep_fh)
print('WARNINGS:', file=rep_fh)
print('-'*20, file=rep_fh)

if warn_str != '':
    rep_fh.write(warn_str)
else:
    rep_fh.write(' None\n')

print('-'*20, file=rep_fh)
print('', file=rep_fh)
    

# ----------------------------------------------------------------------
# Compare coordinates
# ----------------------------------------------------------------------

diffs = {}
dist_list = []
maxDist = 0
minDist = 0

for stn in coords1:
    if stn in coords2:
        dE = coords2[stn]['E'] - coords1[stn]['E']
        dN = coords2[stn]['N'] - coords1[stn]['N']
        dZ = coords2[stn]['Z'] - coords1[stn]['Z']
        dOHGT = coords2[stn]['OHGT'] - coords1[stn]['OHGT']
        dEHGT = coords2[stn]['EHGT'] - coords1[stn]['EHGT']
        dSE = coords2[stn]['SE'] - coords1[stn]['SE']
        dSN = coords2[stn]['SN'] - coords1[stn]['SN']
        dSU = coords2[stn]['SU'] - coords1[stn]['SU']
        dist = math.sqrt(dE**2 + dN**2)
        brg = (180 / math.pi) * math.atan2(dE, dN)
        brg = brg % 360

        if dist < minDist:
            minDist = dist
        if dist > maxDist:
            maxDist = dist

        dist_list.append(dist)

        diffs[stn] = {'dE': dE,
                      'dN': dN,
                      'DIST': dist,
                      'BRG': brg,
                      'dZ': dZ,
                      'dOHGT': dOHGT,
                      'dEHGT': dEHGT,
                      'dSE': dSE,
                      'dSN': dSN,
                      'dSU': dSU
                      }


# ----------------------------------------------------------------------
# Print coordinate differences
# ----------------------------------------------------------------------
common_stns = len(diffs)

print('Station                    Lat (dd)   Lon (dd)        dE        dN        Dist     Brg     dZ    dOHGT     '
      'dEHGT     dSD_E     dSD_N     dSD_U', file=rep_fh)
print('-'*145,file=rep_fh)
for stn in diffs:
    print('{:20s}{:15.8f}{:15.8f}{:10.4f}{:10.4f}{:10.4f}{:8.1f}{:5d}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}'
        .format(stn,
                coords1[stn]['LAT'],
                coords1[stn]['LON'],
                diffs[stn]['dE'],
                diffs[stn]['dN'],
                diffs[stn]['DIST'],
                diffs[stn]['BRG'],
                diffs[stn]['dZ'],
                diffs[stn]['dOHGT'],
                diffs[stn]['dEHGT'],
                diffs[stn]['dSE'],
                diffs[stn]['dSN'],
                diffs[stn]['dSU'],),
          file=rep_fh)

    print('{:s},{:0.8f},{:0.8f},{:0.4f},{:0.4f},{:0.4f},{:0.1f},{:d},{:0.4f},{:0.4f},{:0.4f},{:0.4f},{:0.4f}'
          .format(stn,
                  coords1[stn]['LAT'],
                  coords1[stn]['LON'],
                  diffs[stn]['dE'],
                  diffs[stn]['dN'],
                  diffs[stn]['DIST'],
                  diffs[stn]['BRG'],
                  diffs[stn]['dZ'],
                  diffs[stn]['dOHGT'],
                  diffs[stn]['dEHGT'],
                  diffs[stn]['dSE'],
                  diffs[stn]['dSN'],
                  diffs[stn]['dSU'],),
          file=csv_fh)

rep_fh.close()
csv_fh.close()


# ----------------------------------------------------------------------
#  write shapefile
# ----------------------------------------------------------------------

if common_stns > 0:
    print(" writing point shapefile ...")

    w = shapefile.Writer(shp_file, shapeType=1)    # type 1 for points.

    w.autoBalance = 1

    w.field('Station', 'C', size=20)
    w.field('Latitude', 'N', decimal=10)
    w.field('Longitude', 'N', decimal=10)
    w.field('dEast', 'N', decimal=4)
    w.field('dNorth', 'N', decimal=4)
    w.field('dist', 'N', decimal=4)
    w.field('brg', 'N', decimal=4)
    w.field('dZone', 'N', decimal=4)
    w.field('dOHGT', 'N', decimal=4)
    w.field('dEHGT', 'N', decimal=4)
    w.field('dSD_E', 'N', decimal=4)
    w.field('dSD_N', 'N', decimal=4)
    w.field('dSD_U', 'N', decimal=4)

    # write out point for each dictionary entry
    for stn in diffs:

        w.point(coords1[stn]['LON'], coords1[stn]['LAT'])

        w.record(stn,
                 float(coords1[stn]['LAT']),
                 float(coords1[stn]['LON']),
                 float(diffs[stn]['dE']),
                 float(diffs[stn]['dN']),
                 float(diffs[stn]['DIST']),
                 float(diffs[stn]['BRG']),
                 float(diffs[stn]['dZ']),
                 float(diffs[stn]['dOHGT']),
                 float(diffs[stn]['dEHGT']),
                 float(diffs[stn]['dSE']),
                 float(diffs[stn]['dSN']),
                 float(diffs[stn]['dSU'])
                 )

    w.close()


# ----------------------------------------------------------------------
# Print summary
# ----------------------------------------------------------------------

print()
print('-'*40)
print(' Compared {:d} common stations'.format(common_stns))
print('   Max hz change: {:10.4f}'.format(maxDist))
print('   Min hz change: {:10.4f}'.format(minDist))
print('')
