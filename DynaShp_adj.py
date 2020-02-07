# ----------------------------------------------------------------------
#                          DynaShp.py
# ----------------------------------------------------------------------
#  Author:  Nicholas Gowans
#    Date:  05 February 2020
# Purpose:  Script to create stn/msr shapefiles from DynAdjust .adj file
# ----------------------------------------------------------------------
#   Usage:  cmd:\> python DynaShp_adj.py <*.adj_file>
# ----------------------------------------------------------------------
#   Notes:  - Requires coordinates to be in form ENzPLHhXYZ
#           - Currently handles msr types BDEGHLMRXY. Other types are ignored.
#           - GDA2020 and GDA94 reference frames are supported. Others are
#             written as datum-less.
#
#           Future work:
#               - Include remaining msr types
#               - Introduce additional epsg codes for other reference frames
#               - Handle variable coordinate formats
#               - Introduce functions for shapefile writes (not one per msr type)
# ----------------------------------------------------------------------

import geodepy.convert as gc
import os
import shapefile

adj_file = os.sys.argv[1]


# ----------------------------------------------------------------------
# prepare file names
# ----------------------------------------------------------------------

# pyshp truncates shapefile names at the first period.
count = 0
for l in adj_file:
    count += 1
    if l == '.':
        short = adj_file[:count-1]
        break
    # no period in file name
    if count == len(adj_file):
        short = adj_file

point_name = short + '_stn'
b_msr_name = short + '_b'
d_msr_name = short + '_d'
e_msr_name = short + '_e'
g_msr_name = short + '_g'
h_msr_name = short + '_h'
l_msr_name = short + '_l'
m_msr_name = short + '_m'
r_msr_name = short + '_r'
x_msr_name = short + '_x'
y_msr_name = short + '_y'


# ----------------------------------------------------------------------
# prepare counters/dictionaries
# ----------------------------------------------------------------------

msr_switch = False
stn_switch = False
l_count = 0
h_count = 0
r_count = 0
g_count = 0
x_count = 0
y_count = 0
d_count = 0
d_set = 0
e_count = 0
m_count = 0
b_count = 0

l_msrs = {}
h_msrs = {}
r_msrs = {}
g_msrs = {}
x_msrs = {}
y_msrs = {}
d_msrs = {}
m_msrs = {}
e_msrs = {}
b_msrs = {}
stns = {}


# ----------------------------------------------------------------------
# read adj results
# ----------------------------------------------------------------------

adj_fh = open(adj_file,'r')

print(" reading adj results ...")

for line in adj_fh:

    if line.strip() == '-'*228:
        msr_switch = True
        continue

    if line.strip() == '-'*248:
        stn_switch = True
        continue

    if line[:35] == 'Reference frame:                   ':
        ref_frame = line[35:].strip()

    if line[:35] == 'Station coordinate types:          ':
        coord_types = line[35:].strip()
        if coord_types != 'ENzPLHhXYZ':
            print('\n ERROR!\n\n Coordinate types must be in form ENzPLHhXYZ\n\n Exiting...')
            exit()

    if(msr_switch):
        msrType = line[0:1]

        if msrType =='L':
            l_count += 1
            station1 = line[2:22].strip()
            station2 = line[22:42].strip()
            flag = line[62:63]
            data = line[67:].split()
            msr = data[0]
            adj = data[1]
            cor = data[2]
            msr_SD = data[3]
            adj_SD = data[4]
            res = data[5]
            nStat = data[6]
            pelzer = data[7]
            PreAdjCor = data[8] 
            Outlier = line[204:205]
            
            l_msrs[l_count] = {
                'Stn1': station1,
                'Stn2': station2,
                'Msr': msr,
                'Msr_SD': msr_SD,
                'Adj': adj,
                'Cor': cor,
                'nStat': nStat,
                'outlier':Outlier,
            }

            continue

        if msrType =='B':
            b_count += 1
            station1 = line[2:22].strip()
            station2 = line[22:42].strip()
            flag = line[62:63]
            data = line[67:].split()
            msr = data[0]
            adj = data[1]
            cor = data[2]
            msr_SD = data[3]
            adj_SD = data[4]
            res = data[5]
            nStat = data[6]
            pelzer = data[7]
            PreAdjCor = data[8] 
            Outlier = line[204:205]
            
            b_msrs[b_count] = {
                'Stn1': station1,
                'Stn2': station2,
                'Msr': msr,
                'Msr_SD': msr_SD,
                'Adj': adj,
                'Cor': cor,
                'nStat': nStat,
                'outlier': Outlier,
            }

            continue

        if msrType =='E':
            e_count += 1
            station1 = line[2:22].strip()
            station2 = line[22:42].strip()
            flag = line[62:63]
            data = line[67:].split()
            msr = data[0]
            adj = data[1]
            cor = data[2]
            msr_SD = data[3]
            adj_SD = data[4]
            res = data[5]
            nStat = data[6]
            pelzer = data[7]
            PreAdjCor = data[8] 
            Outlier = line[204:205]
            
            e_msrs[e_count] = {
                'Stn1': station1,
                'Stn2': station2,
                'Msr': msr,
                'Msr_SD': msr_SD,
                'Adj': adj,
                'Cor': cor,
                'nStat': nStat,
                'outlier': Outlier,
            }

            continue

        if msrType =='M':
            m_count += 1
            station1 = line[2:22].strip()
            station2 = line[22:42].strip()
            flag = line[62:63]
            data = line[67:].split()
            msr = data[0]
            adj = data[1]
            cor = data[2]
            msr_SD = data[3]
            adj_SD = data[4]
            res = data[5]
            nStat = data[6]
            pelzer = data[7]
            PreAdjCor = data[8] 
            Outlier = line[204:205]
            
            m_msrs[m_count] = {
                'Stn1': station1,
                'Stn2': station2,
                'Msr': msr,
                'Msr_SD': msr_SD,
                'Adj': adj,
                'Cor': cor,
                'nStat': nStat,
                'outlier': Outlier,
            }

            continue

        if msrType == 'G':
            g_count += 1
            station1 = line[2:22].strip()
            station2 = line[22:42].strip()
            flag = line[62:63]
            data = line[67:].split()

            # reset lists if first in tuple
            if g_count % 3 == 1:
                g_msr = []
                g_adj = []
                g_cor = []
                g_msr_SD = []
                g_adj_SD = []
                g_res = []
                g_nStat = []
                g_pelzer = []
                g_PreAdjCor = []
                g_Outlier = []

            g_msr.append(data[0])
            g_adj.append(data[1])
            g_cor.append(float(data[2]))
            g_msr_SD.append(data[3])
            g_adj_SD.append(data[4])
            g_res.append(data[5])
            g_nStat.append(float(data[6]))
            try:
                g_pelzer.append(data[7])
            except:
                print(line)
                print(data)
                exit()
            g_PreAdjCor.append(data[8])
            g_Outlier.append(line[204:205])

            #  write to dictionary if 3rd in tuple
            if g_count % 3 == 0:
                # find max nStat
                max_nStat = 0.0
                for n in g_nStat:
                    if abs(n) > max_nStat:
                        max_nStat = n
                
                # find max nStat
                max_Cor = 0.0
                for c in g_cor:
                    if abs(c) > max_Cor:
                        max_Cor = c
                
                # write to dictionary
                g_msrs[g_count/3] = {
                    'Stn1': station1,
                    'Stn2': station2,
                    'Msr_X': g_msr[0],
                    'Msr_Y': g_msr[1],
                    'Msr_Z': g_msr[2],
                    'Msr_SD_X': g_msr_SD[0],
                    'Msr_SD_Y': g_msr_SD[1],
                    'Msr_SD_Z': g_msr_SD[2],
                    'Adj_X': g_adj[0],
                    'Adj_Y': g_adj[1],
                    'Adj_Z': g_adj[2],
                    'Cor_X': g_cor[0],
                    'Cor_Y': g_cor[1],
                    'Cor_Z': g_cor[2],
                    'Max_Cor': max_Cor,
                    'nStat_X': g_nStat[0],
                    'nStat_Y': g_nStat[1],
                    'nStat_Z': g_nStat[2],
                    'Max_nStat': max_nStat,
                    'outlier_X': g_Outlier[0],
                    'outlier_Y': g_Outlier[1],
                    'outlier_Z': g_Outlier[2]
                }

            continue

        if msrType == 'X':
            x_count += 1
            station1 = line[2:22].strip()
            station2 = line[22:42].strip()
            flag = line[62:63]
            data = line[67:].split()

            # reset lists if first in tuple
            if x_count % 3 == 1:
                x_msr = []
                x_adj = []
                x_cor = []
                x_msr_SD = []
                x_adj_SD = []
                x_res = []
                x_nStat = []
                x_pelzer = []
                x_PreAdjCor = []
                x_Outlier = []


            x_msr.append(data[0])
            x_adj.append(data[1])
            x_cor.append(float(data[2]))
            x_msr_SD.append(data[3])
            x_adj_SD.append(data[4])
            x_res.append(data[5])
            x_nStat.append(float(data[6]))
            x_pelzer.append(data[7])
            x_PreAdjCor.append(data[8])
            x_Outlier.append(line[204:205])

            #  write to dictionary if 3rd in tuple
            if x_count % 3 == 0:
                # find max nStat
                max_nStat = 0.0
                for n in x_nStat:
                    if abs(n) > max_nStat:
                        max_nStat = n
                
                # find max nStat
                max_Cor = 0.0
                for c in x_cor:
                    if abs(c) > max_Cor:
                        max_Cor = c
                
                # write to dictionary

                x_msrs[x_count/3] = {
                    'Stn1': station1,
                    'Stn2': station2,
                    'Msr_X': x_msr[0],
                    'Msr_Y': x_msr[1],
                    'Msr_Z': x_msr[2],
                    'Msr_SD_X': x_msr_SD[0],
                    'Msr_SD_Y': x_msr_SD[1],
                    'Msr_SD_Z': x_msr_SD[2],
                    'Adj_X': x_adj[0],
                    'Adj_Y': x_adj[1],
                    'Adj_Z': x_adj[2],
                    'Cor_X': x_cor[0],
                    'Cor_Y': x_cor[1],
                    'Cor_Z': x_cor[2],
                    'Max_Cor': max_Cor,
                    'nStat_X': x_nStat[0],
                    'nStat_Y': x_nStat[1],
                    'nStat_Z': x_nStat[2],
                    'Max_nStat': max_nStat,
                    'outlier_X': x_Outlier[0],
                    'outlier_Y': x_Outlier[1],
                    'outlier_Z': x_Outlier[2]
                }

            continue

        if msrType == 'Y':
            y_count += 1
            station1 = line[2:22].strip()
            # station2 = line[22:42].strip()
            flag = line[62:63]
            data = line[67:].split()

            # reset lists if first in tuple
            if y_count % 3 == 1:
                y_msr = []
                y_adj = []
                y_cor = []
                y_msr_SD = []
                y_adj_SD = []
                y_res = []
                y_nStat = []
                y_pelzer = []
                y_PreAdjCor = []
                y_Outlier = []


            y_msr.append(data[0])
            y_adj.append(data[1])
            y_cor.append(float(data[2]))
            y_msr_SD.append(data[3])
            y_adj_SD.append(data[4])
            y_res.append(data[5])
            y_nStat.append(float(data[6]))
            y_pelzer.append(data[7])
            y_PreAdjCor.append(data[8])
            y_Outlier.append(line[204:205])

            #  write to dictionary if 3rd in tuple
            if y_count % 3 == 0:
                # find max nStat
                max_nStat = 0.0
                for n in y_nStat:
                    if abs(n) > max_nStat:
                        max_nStat = n
                
                # find max nStat
                max_Cor = 0.0
                for c in y_cor:
                    if abs(c) > max_Cor:
                        max_Cor = c
                
                # write to dictionary
                y_msrs[y_count/3] = {
                    'Stn1': station1,
                    'Msr_X': y_msr[0],
                    'Msr_Y': y_msr[1],
                    'Msr_Z': y_msr[2],
                    'Msr_SD_X': y_msr_SD[0],
                    'Msr_SD_Y': y_msr_SD[1],
                    'Msr_SD_Z': y_msr_SD[2],
                    'Adj_X': y_adj[0],
                    'Adj_Y': y_adj[1],
                    'Adj_Z': y_adj[2],
                    'Cor_X': y_cor[0],
                    'Cor_Y': y_cor[1],
                    'Cor_Z': y_cor[2],
                    'Max_Cor': max_Cor,
                    'nStat_X': y_nStat[0],
                    'nStat_Y': y_nStat[1],
                    'nStat_Z': y_nStat[2],
                    'Max_nStat': max_nStat,
                    'outlier_X': y_Outlier[0],
                    'outlier_Y': y_Outlier[1],
                    'outlier_Z': y_Outlier[2]
                }

            continue

        if msrType == 'H':
            h_count += 1
            station1 = line[2:22].strip()
            flag = line[62:63]
            data = line[67:].split()
            msr = data[0]
            adj = data[1]
            cor = data[2]
            msr_SD = data[3]
            adj_SD = data[4]
            res = data[5]
            nStat = data[6]
            pelzer = data[7]
            PreAdjCor = data[8]
            Outlier = line[204:205]

            h_msrs[h_count] = {
                'Stn1': station1,
                'Msr': msr,
                'Msr_SD': msr_SD,
                'Adj': adj,
                'Cor': cor,
                'nStat': nStat,
                'outlier': Outlier,
                'matched': False
            }

        if msrType == 'R':
            r_count += 1
            station1 = line[2:22].strip()
            flag = line[62:63]
            data = line[67:].split()
            msr = data[0]
            adj = data[1]
            cor = data[2]
            msr_SD = data[3]
            adj_SD = data[4]
            res = data[5]
            nStat = data[6]
            pelzer = data[7]
            PreAdjCor = data[8]
            Outlier = line[204:205]

            r_msrs[r_count] = {
                'Stn1': station1,
                'Msr': msr,
                'Msr_SD': msr_SD,
                'Adj': adj,
                'Cor': cor,
                'nStat': nStat,
                'outlier': Outlier,
                'matched': False
            }
        
        if msrType == 'D':
            output = []
            d_set +=1
            pointing = 0
            station1 = line[2:22].strip()
            station2 = line[22:42].strip()
            d_msrs[d_set] = {
                'Stn1': station1,
                'Stn2': station2
            }

        if msrType == ' ':
            d_count += 1
            pointing += 1
            target = line[42:63].strip()
            msr = line[72:86].strip()
            adj = line[86:105].strip()
            cor = line[105:117].strip()
            dmsr_SD = line[117:130].strip()
            dadj_SD = line[130:143].strip()
            dres = line[143:156].strip()
            dnStat = line[156:167].strip()
            dpelzer = line[167:179].strip()
            dPreAdjCor = line[179:193].strip()
            doutlier = line[204:205].strip()

            d_msrs[d_set][pointing] = {
                'target': target,
                'msr': msr,
                'adj': adj,
                'cor': cor,
                'msr_SD': dmsr_SD,
                'adj_SD': dadj_SD,
                'res': dres,
                'nStat': dnStat,
                'pelzer': dpelzer,
                'PreAdjCor': dPreAdjCor,
                'outlier': doutlier
            }


        if line == '\n':
            msr_switch = False
            continue
        
    if(stn_switch):
        if len(line) < 20:
            stn_switch = False
            continue

        stn = line[0:20].strip()
        con = line[20:23]
        east = float(line[28:39])
        north = float(line[40:54])
        zone = int(line[60:62])
        lat = gc.hp2dec(float(line[62:76]))
        lon = gc.hp2dec(float(line[77:91]))
        OHGT = float(line[91:102])
        EHGT = float(line[102:113])
        SD_E = float(line[161:170])
        SD_N = float(line[170:180])
        SD_U = float(line[180:190])

        stns[stn] = {
            'con':con,
            'E': east,
            'N': north,
            'Z': zone,
            'P': lat,
            'L': lon,
            'H': OHGT,
            'h': EHGT,
            'SD_E': SD_E,
            'SD_N': SD_N,
            'SD_U': SD_U
        }


# projection function
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


# ----------------------------------------------------------------------
#  write stations
# ----------------------------------------------------------------------
if stns:
    print(" writing stn shapefile ...")

    write_prj(point_name, ref_frame)

    w = shapefile.Writer(point_name, shapeType=1)    # type 1 for points.

    w.autoBalance = 1

    # w.field('ID', 'N')
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



    for s in stns:


        w.point(stns[s]['L'], stns[s]['P'])

        w.record(s, str(stns[s]['con']), float(stns[s]['E']), float(stns[s]['N']), int(stns[s]['Z']), float(stns[s]['P']),
                 float(stns[s]['L']), float(stns[s]['H']), float(stns[s]['h']), float(stns[s]['SD_E']),
                 float(stns[s]['SD_N']), float(stns[s]['SD_U']))

    w.close()

# ----------------------------------------------------------------------
#  write h_msrs
# ----------------------------------------------------------------------

if h_msrs:
    print(" writing type H msr shapefile ...")

    write_prj(h_msr_name, ref_frame)

    w = shapefile.Writer(h_msr_name,shapeType=1)    # type 1 for points.

    w.autoBalance = 1

    # w.field('ID', 'N')
    w.field('Stn1', 'C', size=20)
    w.field('msr', 'N', decimal=4)
    w.field('StdDev', 'N', decimal=4)
    w.field('adj', 'N', decimal=4)
    w.field('cor', 'N', decimal=4)
    w.field('nStat', 'N', decimal=4)

    for h in h_msrs:


        w.point(stns[h_msrs[h]['Stn1']]['L'], stns[h_msrs[h]['Stn1']]['P'])

        w.record(h_msrs[h]['Stn1'], h_msrs[h]['Msr'], h_msrs[h]['Msr_SD'], h_msrs[h]['Adj'],
                 h_msrs[h]['Cor'], h_msrs[h]['nStat'])

    w.close()


# ----------------------------------------------------------------------
#  write r_msrs
# ----------------------------------------------------------------------

if r_msrs:
    print(" writing type R msr shapefile ...")

    write_prj(r_msr_name, ref_frame)

    w = shapefile.Writer(r_msr_name, shapeType=1)    # type 1 for points.

    w.autoBalance = 1

    # w.field('ID', 'N')
    w.field('Stn1', 'C', size=20)
    w.field('msr', 'N', decimal=4)
    w.field('StdDev', 'N', decimal=4)
    w.field('adj', 'N', decimal=4)
    w.field('cor', 'N', decimal=4)
    w.field('nStat', 'N', decimal=4)

    for r in r_msrs:

        w.point(stns[r_msrs[r]['Stn1']]['L'], stns[r_msrs[r]['Stn1']]['P'])

        w.record(r_msrs[r]['Stn1'], r_msrs[r]['Msr'], r_msrs[r]['Msr_SD'], r_msrs[r]['Adj'],
                 r_msrs[r]['Cor'], r_msrs[r]['nStat'])

    w.close()


# ----------------------------------------------------------------------
#  write b_msrs
# ----------------------------------------------------------------------
if b_msrs:
    print(" writing type B msr shapefile ...")

    write_prj(b_msr_name, ref_frame)

    w = shapefile.Writer(b_msr_name, shapeType=3)    # type 3 for polylines.

    w.autoBalance = 1

    w.field('Stn1', 'C', size=20)
    w.field('Stn2', 'C', size=20)
    w.field('msr', 'N', decimal=4)
    w.field('StdDev', 'N', decimal=4)
    w.field('adj', 'N', decimal=4)
    w.field('cor', 'N', decimal=4)
    w.field('nStat', 'N', decimal=4)

    for b in b_msrs:

        lon1 = stns[b_msrs[b]['Stn1']]['L']
        lat1 = stns[b_msrs[b]['Stn1']]['P']

        lon2 = stns[b_msrs[b]['Stn2']]['L']
        lat2 = stns[b_msrs[b]['Stn2']]['P']

        w.line([
            [[lon1, lat1],[lon2, lat2]]
        ])

        w.record(b_msrs[b]['Stn1'], b_msrs[b]['Stn2'], b_msrs[b]['Msr'], b_msrs[b]['Msr_SD'], b_msrs[b]['Adj'],
                 b_msrs[b]['Cor'], b_msrs[b]['nStat'])

    w.close()


# ----------------------------------------------------------------------
#  write l_msrs
# ----------------------------------------------------------------------
if l_msrs:
    print(" writing type L msr shapefile ...")

    write_prj(l_msr_name, ref_frame)

    w = shapefile.Writer(l_msr_name, shapeType=3)    # type 3 for polylines.

    w.autoBalance = 1

    w.field('Stn1', 'C', size=20)
    w.field('Stn2', 'C', size=20)
    w.field('msr', 'N', decimal=4)
    w.field('StdDev', 'N', decimal=4)
    w.field('adj', 'N', decimal=4)
    w.field('cor', 'N', decimal=4)
    w.field('nStat', 'N', decimal=4)

    for l in l_msrs:

        lon1 = stns[l_msrs[l]['Stn1']]['L']
        lat1 = stns[l_msrs[l]['Stn1']]['P']

        lon2 = stns[l_msrs[l]['Stn2']]['L']
        lat2 = stns[l_msrs[l]['Stn2']]['P']

        w.line([
            [[lon1, lat1],[lon2, lat2]]
        ])

        w.record(l_msrs[l]['Stn1'], l_msrs[l]['Stn2'], l_msrs[l]['Msr'], l_msrs[l]['Msr_SD'], l_msrs[l]['Adj'],
                 l_msrs[l]['Cor'], l_msrs[l]['nStat'])

    w.close()

# ----------------------------------------------------------------------
#  write e_msrs
# ----------------------------------------------------------------------
if e_msrs:
    print(" writing type E msr shapefile ...")

    write_prj(e_msr_name, ref_frame)

    w = shapefile.Writer(e_msr_name, shapeType=3)    # type 3 for polylines.

    w.autoBalance = 1

    w.field('Stn1', 'C', size=20)
    w.field('Stn2', 'C', size=20)
    w.field('msr', 'N', decimal=4)
    w.field('StdDev', 'N', decimal=4)
    w.field('adj', 'N', decimal=4)
    w.field('cor', 'N', decimal=4)
    w.field('nStat', 'N', decimal=4)

    for e in e_msrs:

        lon1 = stns[e_msrs[e]['Stn1']]['L']
        lat1 = stns[e_msrs[e]['Stn1']]['P']

        lon2 = stns[e_msrs[e]['Stn2']]['L']
        lat2 = stns[e_msrs[e]['Stn2']]['P']

        w.line([
            [[lon1, lat1],[lon2, lat2]]
        ])

        w.record(e_msrs[e]['Stn1'], e_msrs[e]['Stn2'], e_msrs[e]['Msr'], e_msrs[e]['Msr_SD'], e_msrs[e]['Adj'],
                 e_msrs[e]['Cor'], e_msrs[e]['nStat'])

    w.close()


# ----------------------------------------------------------------------
#  write m_msrs
# ----------------------------------------------------------------------
if m_msrs:
    print(" writing type M msr shapefile ...")

    write_prj(m_msr_name, ref_frame)

    w = shapefile.Writer(m_msr_name, shapeType=3)    # type 3 for polylines.

    w.autoBalance = 1

    w.field('Stn1', 'C', size=20)
    w.field('Stn2', 'C', size=20)
    w.field('msr', 'N', decimal=4)
    w.field('StdDev', 'N', decimal=4)
    w.field('adj', 'N', decimal=4)
    w.field('cor', 'N', decimal=4)
    w.field('nStat', 'N', decimal=4)

    for m in m_msrs:

        lon1 = stns[m_msrs[m]['Stn1']]['L']
        lat1 = stns[m_msrs[m]['Stn1']]['P']

        lon2 = stns[m_msrs[m]['Stn2']]['L']
        lat2 = stns[m_msrs[m]['Stn2']]['P']

        w.line([
            [[lon1, lat1],[lon2, lat2]]
        ])

        w.record(m_msrs[m]['Stn1'], m_msrs[m]['Stn2'], m_msrs[m]['Msr'], m_msrs[m]['Msr_SD'], m_msrs[m]['Adj'],
                 m_msrs[m]['Cor'], m_msrs[m]['nStat'])

    w.close()


# ----------------------------------------------------------------------
#  write g_msrs
# ----------------------------------------------------------------------

if g_msrs:
    print(" writing type G msr shapefile ...")

    write_prj(g_msr_name, ref_frame)

    w = shapefile.Writer(g_msr_name, shapeType=3)    # type 3 for polylines.

    w.autoBalance = 1

    w.field('Stn1', 'C', size=20)
    w.field('Stn2', 'C', size=20)
    w.field('msr_X', 'N', decimal=4)
    w.field('msr_Y', 'N', decimal=4)
    w.field('msr_Z', 'N', decimal=4)
    w.field('StdDev_X', 'N', decimal=4)
    w.field('StdDev_Y', 'N', decimal=4)
    w.field('StdDev_Z', 'N', decimal=4)
    w.field('adj_X', 'N', decimal=4)
    w.field('adj_Y', 'N', decimal=4)
    w.field('adj_Z', 'N', decimal=4)
    w.field('cor_X', 'N', decimal=4)
    w.field('cor_Y', 'N', decimal=4)
    w.field('cor_Z', 'N', decimal=4)
    w.field('max_cor', 'N', decimal=4)
    w.field('nStat_X', 'N', decimal=4)
    w.field('nStat_Y', 'N', decimal=4)
    w.field('nStat_Z', 'N', decimal=4)
    w.field('max_nStat', 'N', decimal=4)


    for g in g_msrs:

        lon1 = stns[g_msrs[g]['Stn1']]['L']
        lat1 = stns[g_msrs[g]['Stn1']]['P']

        lon2 = stns[g_msrs[g]['Stn2']]['L']
        lat2 = stns[g_msrs[g]['Stn2']]['P']

        w.line([
            [[lon1, lat1],[lon2, lat2]]
        ])

        w.record(g_msrs[g]['Stn1'], g_msrs[g]['Stn2'],
                 g_msrs[g]['Msr_X'], g_msrs[g]['Msr_Y'], g_msrs[g]['Msr_Z'],
                 g_msrs[g]['Msr_SD_X'], g_msrs[g]['Msr_SD_Y'], g_msrs[g]['Msr_SD_Z'],
                 g_msrs[g]['Adj_X'], g_msrs[g]['Adj_Y'], g_msrs[g]['Adj_Z'],
                 g_msrs[g]['Cor_X'], g_msrs[g]['Cor_Y'], g_msrs[g]['Cor_Z'], g_msrs[g]['Max_Cor'],
                 g_msrs[g]['nStat_X'], g_msrs[g]['nStat_Y'], g_msrs[g]['nStat_Z'], g_msrs[g]['Max_nStat']
                 )

    w.close()

# ----------------------------------------------------------------------
#  write x_msrs
# ----------------------------------------------------------------------

if x_msrs:
    print(" writing type X msr shapefile ...")

    write_prj(x_msr_name, ref_frame)

    w = shapefile.Writer(x_msr_name, shapeType=3)  # type 3 for polylines.

    w.autoBalance = 1

    w.field('Stn1', 'C', size=20)
    w.field('Stn2', 'C', size=20)
    w.field('msr_X', 'N', decimal=4)
    w.field('msr_Y', 'N', decimal=4)
    w.field('msr_Z', 'N', decimal=4)
    w.field('StdDev_X', 'N', decimal=4)
    w.field('StdDev_Y', 'N', decimal=4)
    w.field('StdDev_Z', 'N', decimal=4)
    w.field('adj_X', 'N', decimal=4)
    w.field('adj_Y', 'N', decimal=4)
    w.field('adj_Z', 'N', decimal=4)
    w.field('cor_X', 'N', decimal=4)
    w.field('cor_Y', 'N', decimal=4)
    w.field('cor_Z', 'N', decimal=4)
    w.field('max_cor', 'N', decimal=4)
    w.field('nStat_X', 'N', decimal=4)
    w.field('nStat_Y', 'N', decimal=4)
    w.field('nStat_Z', 'N', decimal=4)
    w.field('max_nStat', 'N', decimal=4)


    for x in x_msrs:
        lon1 = stns[x_msrs[x]['Stn1']]['L']
        lat1 = stns[x_msrs[x]['Stn1']]['P']

        lon2 = stns[x_msrs[x]['Stn2']]['L']
        lat2 = stns[x_msrs[x]['Stn2']]['P']

        w.line([
            [[lon1, lat1], [lon2, lat2]]
        ])

        w.record(x_msrs[x]['Stn1'], x_msrs[x]['Stn2'],
                 x_msrs[x]['Msr_X'], x_msrs[x]['Msr_Y'], x_msrs[x]['Msr_Z'],
                 x_msrs[x]['Msr_SD_X'], x_msrs[x]['Msr_SD_Y'], x_msrs[x]['Msr_SD_Z'],
                 x_msrs[x]['Adj_X'], x_msrs[x]['Adj_Y'], x_msrs[x]['Adj_Z'],
                 x_msrs[x]['Cor_X'], x_msrs[x]['Cor_Y'], x_msrs[x]['Cor_Z'], x_msrs[x]['Max_Cor'],
                 x_msrs[x]['nStat_X'], x_msrs[x]['nStat_Y'], x_msrs[x]['nStat_Z'], x_msrs[x]['Max_nStat']
                 )

    w.close()

# ----------------------------------------------------------------------
#  write y_msrs
# ----------------------------------------------------------------------

if y_msrs:
    print(" writing type Y msr shapefile ...")

    write_prj(y_msr_name, ref_frame)

    w = shapefile.Writer(y_msr_name, shapeType=1)  # type 1 for points.

    w.autoBalance = 1

    w.field('Stn1', 'C', size=20)
    # w.field('Stn2', 'C', size=20)
    w.field('msr_X', 'N', decimal=4)
    w.field('msr_Y', 'N', decimal=4)
    w.field('msr_Z', 'N', decimal=4)
    w.field('StdDev_X', 'N', decimal=4)
    w.field('StdDev_Y', 'N', decimal=4)
    w.field('StdDev_Z', 'N', decimal=4)
    w.field('adj_X', 'N', decimal=4)
    w.field('adj_Y', 'N', decimal=4)
    w.field('adj_Z', 'N', decimal=4)
    w.field('cor_X', 'N', decimal=4)
    w.field('cor_Y', 'N', decimal=4)
    w.field('cor_Z', 'N', decimal=4)
    w.field('max_cor', 'N', decimal=4)
    w.field('nStat_X', 'N', decimal=4)
    w.field('nStat_Y', 'N', decimal=4)
    w.field('nStat_Z', 'N', decimal=4)
    w.field('max_nStat', 'N', decimal=4)

    for y in y_msrs:
        lon1 = stns[y_msrs[y]['Stn1']]['L']
        lat1 = stns[y_msrs[y]['Stn1']]['P']

        w.point(lon1, lat1)

        w.record(y_msrs[y]['Stn1'],
                 y_msrs[y]['Msr_X'], y_msrs[y]['Msr_Y'], y_msrs[y]['Msr_Z'],
                 y_msrs[y]['Msr_SD_X'], y_msrs[y]['Msr_SD_Y'], y_msrs[y]['Msr_SD_Z'],
                 y_msrs[y]['Adj_X'], y_msrs[y]['Adj_Y'], y_msrs[y]['Adj_Z'],
                 y_msrs[y]['Cor_X'], y_msrs[y]['Cor_Y'], y_msrs[y]['Cor_Z'],y_msrs[y]['Max_Cor'],
                 y_msrs[y]['nStat_X'], y_msrs[y]['nStat_Y'], y_msrs[y]['nStat_Z'], y_msrs[y]['Max_nStat']
                 )

    w.close()


# ----------------------------------------------------------------------
#  write d_msrs
# ----------------------------------------------------------------------

if d_msrs:
    print(" writing type D msr shapefile ...")

    write_prj(d_msr_name, ref_frame)

    w = shapefile.Writer(d_msr_name, shapeType=3)  # type 3 for polylines.

    w.autoBalance = 1

    w.field('D_set', 'N', decimal=0)
    w.field('Stn1', 'C', size=20)
    w.field('Stn2', 'C', size=20)
    w.field('msr', 'C' )
    w.field('StdDev', 'N', decimal=4)
    w.field('adj', 'C')
    w.field('adj_SD', 'N', decimal=4)
    w.field('cor', 'N', decimal=4)
    w.field('nStat', 'N', decimal=4)
    
    for d_set in d_msrs:
        stn1 = d_msrs[d_set]['Stn1']
        stn2 = d_msrs[d_set]['Stn2']
        lon1 = stns[stn1]['L']
        lat1 = stns[stn1]['P']

        lon2 = stns[stn2]['L']
        lat2 = stns[stn2]['P']

        w.line([
            [[lon1, lat1], [lon2, lat2]]
        ])

        # set zeroes for direction set reference.
        w.record(d_set, stn1, stn2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        pointings = len(d_msrs[d_set])-2

        # note, pointing key starts at 1
        for i in range(1, pointings+1):
            target = d_msrs[d_set][i]['target']
            lonT = stns[target]['L']
            latT = stns[target]['P']
            msr = d_msrs[d_set][i]['msr']
            msr_SD = d_msrs[d_set][i]['msr_SD']
            adj = d_msrs[d_set][i]['adj']
            adj_SD = d_msrs[d_set][i]['adj_SD']
            cor = d_msrs[d_set][i]['cor']
            nStat = d_msrs[d_set][i]['nStat']

            w.line([
                [[lon1, lat1], [lonT, latT]]
            ])

            # set zeroes for direction set reference.
            w.record(d_set, stn1, target, msr, msr_SD, adj, adj_SD, cor, nStat)

    w.close()
