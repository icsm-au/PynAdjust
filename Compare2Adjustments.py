# ----------------------------------------------------------------------
#                          Compare2Adjustments.py
# ----------------------------------------------------------------------
#  Author:  Kent Wheeler
#    Date:  12 March 2020
# Purpose:  Script to create stn/msr shapefiles from DynAdjust .adj file
# ----------------------------------------------------------------------
#   Usage:  cmd:\> python Compare2Adjustments.py <*.adj_file1> <*.adj_file2>
#           Or
#           cmd:\> python Compare2Adjustments.py <*.adj_file1>
# ----------------------------------------------------------------------
#   Notes:  - Currently handles the following msr types and produces these shp files
#                         Azimuths_TYPE_B_K_V_Z
#                         Angles_TYPE_D_A
#                         Distances_TYPE_C_E_M_S
#                         GNSS_TYPE_G_X_Y
#                         Height_Diff_TYPE_L
#                         Heights_TYPE_H_R
#                         Astro_TYPE_P_Q_I_J
#                          Other types are ignored.
#           - GDA2020 and GDA94 reference frames are supported.
#           - angular_msr_types are required to be in the same format in each file
#

import os, sys,sqlite3
import shapefile as shp
import math

# WGS 84
a = 6378137  # meters
f = 1 / 298.257222101
b = (1 - f)*a

MAX_ITERATIONS = 200
CONVERGENCE_THRESHOLD = 1e-12  # .000,000,000,001

def vincenty_inverse(point1, point2, miles=False):
    # short-circuit coincident points
    if point1[0] == point2[0] and point1[1] == point2[1]:
        return [0.0,0.0]

    U1 = math.atan((1 - f) * math.tan(math.radians(point1[0])))
    U2 = math.atan((1 - f) * math.tan(math.radians(point2[0])))
    L = math.radians(point2[1] - point1[1])
    Lambda = L

    sinU1 = math.sin(U1)
    cosU1 = math.cos(U1)
    sinU2 = math.sin(U2)
    cosU2 = math.cos(U2)
    s12 = 0.0
    Alpha12 = 0.0
    for iteration in range(MAX_ITERATIONS):
        sinLambda = math.sin(Lambda)
        cosLambda = math.cos(Lambda)
        sinSigma = math.sqrt((cosU2 * sinLambda) ** 2 +
                             (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ** 2)
        if sinSigma == 0:
            return 0.0,0.0  # coincident points
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
        sigma = math.atan2(sinSigma, cosSigma)
        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
        cosSqAlpha = 1 - sinAlpha ** 2
        try:
            cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha
        except ZeroDivisionError:
            cos2SigmaM = 0
        C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
        LambdaPrev = Lambda
        Lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma *
                                               (cos2SigmaM + C * cosSigma *
                                                (-1 + 2 * cos2SigmaM ** 2)))
        if abs(Lambda - LambdaPrev) < CONVERGENCE_THRESHOLD:
            break  # successful convergence
    else:
        return None  # failure to converge

    uSq = cosSqAlpha * (a ** 2 - b ** 2) / (b ** 2)
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma *
                 (-1 + 2 * cos2SigmaM ** 2) - B / 6 * cos2SigmaM *
                 (-3 + 4 * sinSigma ** 2) * (-3 + 4 * cos2SigmaM ** 2)))
    s12 = b * A * (sigma - deltaSigma)
    
    Alpha12 = math.degrees(math.atan2((cosU2*sinLambda),(cosU1*sinU2-sinU1*cosU2*cosLambda)))
    
    if Alpha12 < 0:
        Alpha12 = Alpha12+360.0
    return [round(s12, 4), round(Alpha12,5)]
    
def create_connection(db_file):
    """ create a database connection to a SQLite database """
    try:
        conn = sqlite3.connect(db_file)
        print(sqlite3.version)
    except:
        print ('Cannot create a database') 
    finally:
        conn.close()
        
def dec2hp(dec):
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    hp = degree + (minute / 100) + (second / 10000)
    return hp if dec >= 0 else -hp

def hp2dec(hp):
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    dec = degree + (minute / 60) + (second / 360)
    return dec if hp >= 0 else -dec

def hms2dd(HMS_Ang):
    #Input: HH MM SS.ssss used by Geolab
    #Output: HH.MMSSSsssss used by DynAdjust
    sign=1
    if HMS_Ang.find(' ')==-1:
        return float(HMS_Ang)
    else:
        if HMS_Ang.find('S')!=-1 or HMS_Ang.find('-')!=-1:
            sign=-1
        while HMS_Ang.find('  ')!=-1:
            HMS_Ang=HMS_Ang.replace('  ',' ')
        HMS_Ang=HMS_Ang.replace('S','')
        HMS_Ang=HMS_Ang.replace('E','')
        aAng=HMS_Ang.split()
        return sign*abs(int(float(aAng[0]))) +  float(aAng[1])/60 + float(aAng[2])/3600

###############################################################################
######################## Compare two adjustments ##############################
###############################################################################

adj_file1 ='20200320_(03)_WA_GDA2020.phased-stage.adj'
adj_file2 ='20200330_(03)_WA_GDA2020.phased-stage.adj'
if len(sys.argv)==2: adj_file1=sys.argv[1]; adj_file2 =sys.argv[1]
if len(sys.argv)==3: adj_file1=sys.argv[1]; adj_file2 =sys.argv[2]
adj_files=[adj_file1,adj_file2]

apu_file1 = adj_file1.replace('.adj','_typeB.apu') if os.path.isfile(adj_file1.replace('.adj','_typeB.apu')) else adj_file1.replace('.adj','.apu')
apu_file2 = adj_file2.replace('.adj','_typeB.apu') if os.path.isfile(adj_file2.replace('.adj','_typeB.apu')) else adj_file2.replace('.adj','.apu')
apu_files=[apu_file1,apu_file2]

###############################################################################
############################# CREATE A DATABASE ###############################
###############################################################################
print ('Creating a database of adjustments ......')

dbname = "Compare_Adj.db"
if not os.path.exists(dbname):
    create_connection(dbname)
    conn = sqlite3.connect(dbname) # or use :memory: to put it in RAM
else:
    conn = sqlite3.connect(dbname) # or use :memory: to put it in RAM     
cursor = conn.cursor()
###############################################################################
######################### CREATE TABLES FOR ADJUSMENT RECORDS #################
###############################################################################
print ('Create a few tables for Marks and Obs ... ')
for i in range(1, 3):
    cursor.execute('DROP TABLE IF EXISTS ADJ'+str(i)+'_POINTS')
    cursor.execute("""CREATE TABLE IF NOT EXISTS ADJ"""+str(i)+"""_POINTS (
                          ID INTEGER,
                          STATION TEXT PRIMARY KEY NOT NULL,
                          CONST TEXT,
                          EASTING REAL,
                          NORTHING REAL,
                          ZONE REAL,
                          LATITUDE REAL,
                          LONGITUDE REAL,
                          H_ORTHO REAL,
                          h_ELLIPSE REAL,
                          X REAL,
                          Y REAL,
                          Z REAL,
                          X_SDEV REAL,
                          Y_SDEV REAL,
                          Z_SDEV REAL,
                          DESC TEXT);""")
    conn.commit()
    cursor.execute('DROP TABLE IF EXISTS APU'+str(i)+'_POINTS')
    cursor.execute("""CREATE TABLE IF NOT EXISTS APU"""+str(i)+"""_POINTS (
                          ID INTEGER,
                          STATION TEXT PRIMARY KEY NOT NULL,
                          LATITUDE REAL,
                          LONGITUDE REAL,
                          Hz_PosU REAL,
                          Vt_PosU REAL,
                          Semi_major REAL,
                          Semi_minor REAL,
                          Orientation REAL);""")
    conn.commit()
    cursor.execute('DROP TABLE IF EXISTS ADJ'+str(i)+'_MEASUREMENTS')
    cursor.execute("""CREATE TABLE IF NOT EXISTS ADJ"""+str(i)+"""_MEASUREMENTS (
                          ID INTEGER,
                          MEAS_TYPE TEXT,
                          STN1 TEXT,
                          STN2 TEXT,
                          STN3 TEXT,
                          C TEXT,
                          VALUE REAL,
                          ADJUSTED_VAL REAL,
                          CORRECTION REAL,
                          MEAS_SDEV REAL,
                          ADJ_SDEV REAL,
                          RESIDUAL REAL,
                          N_STAT REAL,
                          PELZER_REL REAL,
                          PRE_ADJ_CORR REAL,
                          OUTLIER TEXT,
                          UNIQUE_ID TEXT PRIMARY KEY NOT NULL);""")
    conn.commit()

#######################################################################
### strip the adjustment results from the adj files and into Tables ###
#######################################################################
file_count=0
rfs=[]
angular_msr_type='hms'
print ('Reading the adj files ... ')
for fn in adj_files:
    print (fn)
    count_lines=0; file_count=file_count +1
    read_coords = 'false' ; read_meas = 'false'
    prevline='' ; easting = '' ; northing = '' ; zone=''; angular_msr_type=''
    with open(fn, 'r') as file:
        dnadata = file.readlines()
        for line in dnadata:
            if line.find('Reference frame:')!=-1:rfs.append(line[35:])
            if line.find('Command line arguments:')!=-1 and line.find('--angular-msr-type 1')!=-1:angular_msr_type='dd'
            if line.find('Adjusted Coordinates')!=-1:
                read_coords = 'true'
                count_lines=0
                conn.commit()
            if len(line) <=1 and count_lines>=5: 
                read_coords = 'false'
                read_meas = 'false'
            if read_coords =='true' and count_lines==3:
                grid_clm=line.find('       Easting')
                geog_clm=line.find('      Latitude')
                cart_clm=line.find('             X              Y              Z')
                desc_clm=line.find('Description')
            if read_coords =='true' and count_lines>=5:
                pointid = line[0:20]
                const = line[20:25]
                if grid_clm!=-1:
                    easting = float(line[grid_clm:grid_clm+14])
                    northing = float(line[grid_clm+14:grid_clm+29])
                    zone = int(line[grid_clm+29:grid_clm+37])
                if geog_clm!=-1:
                    lat = float(line[geog_clm:geog_clm+14])
                    long = float(line[geog_clm+14:geog_clm+29])
                    ortho_ht = float(line[geog_clm+29:geog_clm+40])
                    ell_ht = float(line[geog_clm+40:geog_clm+51])
                if cart_clm!=-1:
                    x = float(line[cart_clm:cart_clm+15])
                    y = float(line[cart_clm+16:cart_clm+30])
                    z = float(line[cart_clm+31:cart_clm+45])
                    x_sd = float(line[cart_clm+45:cart_clm+57])
                    y_sd = float(line[cart_clm+57:cart_clm+67])
                    z_sd = float(line[cart_clm+67:cart_clm+77])
                desc = line[desc_clm:-1]
                cursor.execute("INSERT OR REPLACE INTO ADJ" + str(file_count) + "_POINTS (ID, STATION, CONST, EASTING, NORTHING, ZONE, LATITUDE, LONGITUDE, H_ORTHO, h_ELLIPSE, X, Y, Z, X_SDEV, Y_SDEV, Z_SDEV, DESC) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", 
                                [count_lines, pointid.strip(), const.strip(), easting, northing, zone, hp2dec(lat),hp2dec(long), ortho_ht, ell_ht, x, y, z, x_sd, y_sd, z_sd, desc.strip()])
            if line.find('Adjusted Measurements')!=-1:
                read_meas = 'true'
                count_lines=0
            if read_meas =='true' and count_lines>=5:
                line=line.replace('-nan(ind)','      0  ')
                if line[0:32]!='                                ':
                    typ = line[0:2].strip()
                    stn1 = line[2:22].strip()
                    stn2 = line[22:42].strip()
                if line[0:2].strip()!='D':
                    stn3 = line[42:63].strip()
                    C = line[65:67].strip()
                    if typ =='D' or typ =='B' or typ =='K':
                        if angular_msr_type=='dd':
                            value = float(line[68:87])
                            Adjusted = float(line[87:105])
                        else:
                            value = hms2dd(line[68:87])
                            Adjusted = hms2dd(line[87:105])
                    else:
                        value = float(line[68:87])
                        Adjusted = float(line[87:105])
                    Correction = float(line[106:118])
                    Meas_SD = float(line[118:131])
                    Adj_SD = float(line[131:144])
                    Residual = float(line[144:157])
                    N_stat = float(line[157:168])
                    Pelzer_Rel = float(line[168:180])
                    Pre_Adj_Corr = float(line[180:194])
                    Outlier = line[194:206].strip()
                    uniqueid=typ + stn1 + stn2 + stn3 + str(value) + str(Meas_SD)
                    cursor.execute("INSERT OR REPLACE INTO ADJ" + str(file_count) + "_MEASUREMENTS (ID, MEAS_TYPE, STN1, STN2, STN3, C, VALUE, ADJUSTED_VAL, CORRECTION, MEAS_SDEV, ADJ_SDEV, RESIDUAL, N_STAT, PELZER_REL, PRE_ADJ_CORR, OUTLIER, UNIQUE_ID) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", 
                                   [count_lines, typ, stn1,stn2,stn3,C,value,Adjusted,Correction,Meas_SD,Adj_SD,Residual,N_stat,Pelzer_Rel,Pre_Adj_Corr,Outlier,uniqueid])
            count_lines = count_lines + 1
            prevline=line
conn.commit()

###############################################################################
##################### Open and read the APU files #############################
###############################################################################
print ('Reading the apu files ... ')
file_count=0
for fn in apu_files:
    print (fn)
    count_lines=0; file_count=file_count +1
    read_coords = 'false'
    if os.path.exists(fn):
        with open(fn, 'r') as file:
            dnadata = file.readlines()
            for line in dnadata:
                if line.find('Positional uncertainty of adjusted station coordinates')!=-1:
                    read_coords = 'true'
                    count_lines=0
                if read_coords =='true' and count_lines>=5 and len(line)>1:
                    if line[0:20]!='                    ':
                        pointid = line[0:20]
                        lat = float(line[20:37])
                        long = float(line[37:52])
                        Hz_PU= float(line[52:63])
                        Vt_PU= float(line[63:74])
                        S_Maj= float(line[74:87])
                        S_Min= float(line[87:100])
                        Orien= float(line[100:113])
                    cursor.execute("INSERT OR REPLACE INTO APU" + str(file_count) + "_POINTS (ID, STATION, LATITUDE, LONGITUDE, Hz_PosU, Vt_PosU, Semi_major, Semi_minor, Orientation) VALUES (?,?,?,?,?,?,?,?,?)", 
                                    [count_lines, pointid.strip(), hp2dec(lat),hp2dec(long), Hz_PU, Vt_PU, S_Maj, S_Min, Orien])
                count_lines = count_lines + 1
conn.commit()
            
print ('Now creating an index on tables .....')
for i in range(1, 3):
    cursor.execute("CREATE UNIQUE INDEX IF NOT EXISTS idxPoints ON ADJ"+str(i)+"_POINTS (ID ASC)")
    conn.commit()
    cursor.execute("CREATE UNIQUE INDEX IF NOT EXISTS idxPoints ON ADJ"+str(i)+"_MEASUREMENTS (ID ASC)")
    conn.commit()
    cursor.execute("CREATE UNIQUE INDEX IF NOT EXISTS idxPoints ON APU"+str(i)+"_POINTS (ID ASC)")
    conn.commit()
shp_sql=[]
shp_fn=[]
shp_fn.append("Removed_Measurements")
shp_sql.append("SELECT ADJ1_MEASUREMENTS.*, ADJ1_POINTS.LATITUDE AS LATITUDE1, ADJ1_POINTS.LONGITUDE AS LONGITUDE1, ADJ1_POINTS_1.LATITUDE AS LATITUDE2, ADJ1_POINTS_1.LONGITUDE AS LONGITUDE2, ADJ1_POINTS_2.LATITUDE AS LATITUDE3, ADJ1_POINTS_2.LONGITUDE AS LONGITUDE3 \
                    FROM (((ADJ1_MEASUREMENTS LEFT JOIN ADJ2_MEASUREMENTS ON ADJ1_MEASUREMENTS.UNIQUE_ID = ADJ2_MEASUREMENTS.UNIQUE_ID) LEFT JOIN ADJ1_POINTS ON ADJ1_MEASUREMENTS.STN1 = ADJ1_POINTS.STATION) LEFT JOIN ADJ1_POINTS AS ADJ1_POINTS_1 ON ADJ1_MEASUREMENTS.STN2 = ADJ1_POINTS_1.STATION) LEFT JOIN ADJ1_POINTS AS ADJ1_POINTS_2 ON ADJ1_MEASUREMENTS.STN3 = ADJ1_POINTS_2.STATION \
                    WHERE (((ADJ2_MEASUREMENTS.UNIQUE_ID) Is Null));")
shp_fn.append("added_Measurements")
shp_sql.append("SELECT ADJ2_MEASUREMENTS.*, ADJ2_POINTS.LATITUDE AS LATITUDE1, ADJ2_POINTS.LONGITUDE AS LONGITUDE1, ADJ2_POINTS_1.LATITUDE AS LATITUDE2, ADJ2_POINTS_1.LONGITUDE AS LONGITUDE2, ADJ2_POINTS_2.LATITUDE AS LATITUDE3, ADJ2_POINTS_2.LONGITUDE AS LONGITUDE3 \
                    FROM (((ADJ2_MEASUREMENTS LEFT JOIN ADJ1_MEASUREMENTS ON ADJ2_MEASUREMENTS.UNIQUE_ID = ADJ1_MEASUREMENTS.UNIQUE_ID) LEFT JOIN ADJ2_POINTS ON ADJ2_MEASUREMENTS.STN1 = ADJ2_POINTS.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_1 ON ADJ2_MEASUREMENTS.STN2 = ADJ2_POINTS_1.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_2 ON ADJ2_MEASUREMENTS.STN3 = ADJ2_POINTS_2.STATION \
                    WHERE (((ADJ1_MEASUREMENTS.UNIQUE_ID) Is Null));")
shp_fn.append("Azimuths_TYPE_B_K_V_Z")
shp_sql.append("SELECT ADJ2_MEASUREMENTS.*, ADJ2_POINTS.LATITUDE, ADJ2_POINTS.LONGITUDE, ADJ2_POINTS_1.LATITUDE, ADJ2_POINTS_1.LONGITUDE, ADJ2_POINTS_2.LATITUDE, ADJ2_POINTS_2.LONGITUDE \
                    FROM ((ADJ2_MEASUREMENTS LEFT JOIN ADJ2_POINTS ON ADJ2_MEASUREMENTS.STN1 = ADJ2_POINTS.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_1 ON ADJ2_MEASUREMENTS.STN2 = ADJ2_POINTS_1.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_2 ON ADJ2_MEASUREMENTS.STN3 = ADJ2_POINTS_2.STATION \
                    WHERE (((ADJ2_MEASUREMENTS.MEAS_TYPE)='B' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='K' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='V' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='Z'));")
shp_fn.append("Angles_TYPE_D_A")
shp_sql.append("SELECT ADJ2_MEASUREMENTS.*, ADJ2_POINTS.LATITUDE, ADJ2_POINTS.LONGITUDE, ADJ2_POINTS_1.LATITUDE, ADJ2_POINTS_1.LONGITUDE, ADJ2_POINTS_2.LATITUDE, ADJ2_POINTS_2.LONGITUDE \
                    FROM ((ADJ2_MEASUREMENTS LEFT JOIN ADJ2_POINTS ON ADJ2_MEASUREMENTS.STN1 = ADJ2_POINTS.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_1 ON ADJ2_MEASUREMENTS.STN2 = ADJ2_POINTS_1.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_2 ON ADJ2_MEASUREMENTS.STN3 = ADJ2_POINTS_2.STATION \
                    WHERE (((ADJ2_MEASUREMENTS.MEAS_TYPE)='D' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='A'));")
shp_fn.append("Distances_TYPE_C_E_M_S")
shp_sql.append("SELECT ADJ2_MEASUREMENTS.*, ADJ2_POINTS.LATITUDE, ADJ2_POINTS.LONGITUDE, ADJ2_POINTS_1.LATITUDE, ADJ2_POINTS_1.LONGITUDE, ADJ2_POINTS_2.LATITUDE, ADJ2_POINTS_2.LONGITUDE \
                    FROM ((ADJ2_MEASUREMENTS LEFT JOIN ADJ2_POINTS ON ADJ2_MEASUREMENTS.STN1 = ADJ2_POINTS.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_1 ON ADJ2_MEASUREMENTS.STN2 = ADJ2_POINTS_1.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_2 ON ADJ2_MEASUREMENTS.STN3 = ADJ2_POINTS_2.STATION \
                    WHERE (((ADJ2_MEASUREMENTS.MEAS_TYPE)='C' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='E' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='M' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='S'));")
shp_fn.append("GNSS_TYPE_G_X_Y")
shp_sql.append("SELECT ADJ2_MEASUREMENTS.*, ADJ2_POINTS.LATITUDE, ADJ2_POINTS.LONGITUDE, ADJ2_POINTS_1.LATITUDE, ADJ2_POINTS_1.LONGITUDE, ADJ2_POINTS_2.LATITUDE, ADJ2_POINTS_2.LONGITUDE \
                    FROM ((ADJ2_MEASUREMENTS LEFT JOIN ADJ2_POINTS ON ADJ2_MEASUREMENTS.STN1 = ADJ2_POINTS.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_1 ON ADJ2_MEASUREMENTS.STN2 = ADJ2_POINTS_1.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_2 ON ADJ2_MEASUREMENTS.STN3 = ADJ2_POINTS_2.STATION \
                    WHERE (((ADJ2_MEASUREMENTS.MEAS_TYPE)='G' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='X' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='Y'));")
shp_fn.append("Height_Diff_TYPE_L")
shp_sql.append("SELECT ADJ2_MEASUREMENTS.*, ADJ2_POINTS.LATITUDE, ADJ2_POINTS.LONGITUDE, ADJ2_POINTS_1.LATITUDE, ADJ2_POINTS_1.LONGITUDE, ADJ2_POINTS_2.LATITUDE, ADJ2_POINTS_2.LONGITUDE \
                    FROM ((ADJ2_MEASUREMENTS LEFT JOIN ADJ2_POINTS ON ADJ2_MEASUREMENTS.STN1 = ADJ2_POINTS.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_1 ON ADJ2_MEASUREMENTS.STN2 = ADJ2_POINTS_1.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_2 ON ADJ2_MEASUREMENTS.STN3 = ADJ2_POINTS_2.STATION \
                    WHERE (((ADJ2_MEASUREMENTS.MEAS_TYPE)='L'));")
shp_fn.append("Heights_TYPE_H_R")
shp_sql.append("SELECT ADJ2_MEASUREMENTS.*, ADJ2_POINTS.LATITUDE, ADJ2_POINTS.LONGITUDE, ADJ2_POINTS_1.LATITUDE, ADJ2_POINTS_1.LONGITUDE, ADJ2_POINTS_2.LATITUDE, ADJ2_POINTS_2.LONGITUDE \
                    FROM ((ADJ2_MEASUREMENTS LEFT JOIN ADJ2_POINTS ON ADJ2_MEASUREMENTS.STN1 = ADJ2_POINTS.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_1 ON ADJ2_MEASUREMENTS.STN2 = ADJ2_POINTS_1.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_2 ON ADJ2_MEASUREMENTS.STN3 = ADJ2_POINTS_2.STATION \
                    WHERE (((ADJ2_MEASUREMENTS.MEAS_TYPE)='H' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='R'));")
shp_fn.append("Astro_TYPE_P_Q_I_J")
shp_sql.append("SELECT ADJ2_MEASUREMENTS.*, ADJ2_POINTS.LATITUDE, ADJ2_POINTS.LONGITUDE, ADJ2_POINTS_1.LATITUDE, ADJ2_POINTS_1.LONGITUDE, ADJ2_POINTS_2.LATITUDE, ADJ2_POINTS_2.LONGITUDE \
                    FROM ((ADJ2_MEASUREMENTS LEFT JOIN ADJ2_POINTS ON ADJ2_MEASUREMENTS.STN1 = ADJ2_POINTS.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_1 ON ADJ2_MEASUREMENTS.STN2 = ADJ2_POINTS_1.STATION) LEFT JOIN ADJ2_POINTS AS ADJ2_POINTS_2 ON ADJ2_MEASUREMENTS.STN3 = ADJ2_POINTS_2.STATION \
                    WHERE (((ADJ2_MEASUREMENTS.MEAS_TYPE)='P' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='Q' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='I' Or (ADJ2_MEASUREMENTS.MEAS_TYPE)='J'));")

##########################################################################
############ Print the queries to some shape files #######################
##########################################################################
print ('Printing the shape files ......')
if not os.path.exists(r'shp Layers'):
    os.mkdir(r'shp Layers')
i=0
for fn in shp_fn:
    f1 = r'./shp Layers/'+fn+'_line'; w1 = shp.Writer(f1)
    w1.field('ID','C'); w1.field('MEAS_TYPE', 'C'); w1.field('STN1','C'); w1.field('STN2','C'); w1.field('STN3','C'); w1.field('C','C');
    w1.field('VALUE', 'N', decimal=10); w1.field('ADJUSTED_VAL', 'N', decimal = 10); w1.field('CORRECTION', 'N', decimal = 10);
    w1.field('MEAS_SDEV', 'N', decimal=10); w1.field('ADJ_SDEV', 'N', decimal = 10); w1.field('RESIDUAL', 'N', decimal = 10);
    w1.field('N_STAT', 'N', decimal=10); w1.field('PELZER_REL', 'N', decimal = 10); w1.field('PRE_ADJ_CORR', 'N', decimal = 10); w1.field('OUTLIER', 'C');
    qry=cursor.execute(shp_sql[i]).fetchall()
    row_cnt=0
    for row in qry:    
        if not row[20]is None and row[22] is None:
            w1.line([
                    [[row[18], row[17]], [row[20], row[19]]]
                    ])
        if not row[22] is None:
            w1.line([
                    [[row[20], row[19]], [row[18], row[17]], [row[22], row[21]]]
                    ])
        if not row[20] is None:
            w1.record(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13], row[14], row[15])
            row_cnt=row_cnt+1
    w1.close()
    if row_cnt==0:os.remove(f1+'.dbf');os.remove(f1+'.shp');os.remove(f1+'.shx')
    else:
        with open(f1 +'.prj', 'w') as prj:
            if shp_sql[i].find('ADJ1_POINTS')!=-1:rf=rfs[0]
            else:rf=rfs[1]
            if rf=='GDA2020\n':prj.write("GEOGCS['GDA2020',DATUM['D_GDA_2020',SPHEROID['GRS_1980',6378137,298.257222101]],PRIMEM['Greenwich',0],UNIT['Degree',0.017453292519943295]]")
            else: prj.write("GEOGCS['GCS_GDA_1994',DATUM['D_GDA_1994',SPHEROID['GRS_1980',6378137,298.257222101]],PRIMEM['Greenwich',0],UNIT['Degree',0.017453292519943295]]")
    
    f1 = r'./shp Layers/'+fn+'_point'; w1 = shp.Writer(f1)
    w1.field('ID','C'); w1.field('MEAS_TYPE', 'C'); w1.field('STN1','C'); w1.field('STN2','C'); w1.field('STN3','C'); w1.field('C','C');
    w1.field('VALUE', 'N', decimal=10); w1.field('ADJUSTED_VAL', 'N', decimal = 10); w1.field('CORRECTION', 'N', decimal = 10);
    w1.field('MEAS_SDEV', 'N', decimal=10); w1.field('ADJ_SDEV', 'N', decimal = 10); w1.field('RESIDUAL', 'N', decimal = 10);
    w1.field('N_STAT', 'N', decimal=10); w1.field('PELZER_REL', 'N', decimal = 10); w1.field('PRE_ADJ_CORR', 'N', decimal = 10); w1.field('OUTLIER', 'C');
    qry=cursor.execute(shp_sql[i]).fetchall()
    row_cnt=0
    for row in qry:
        if row[20] is None and row[22] is None:
            w1.point(row[18], row[17])
            w1.record(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13], row[14], row[15])
            row_cnt=row_cnt+1
    w1.close()
    if row_cnt==0:os.remove(f1+'.dbf');os.remove(f1+'.shp');os.remove(f1+'.shx')
    else:
        with open(f1 +'.prj', 'w') as prj:
            if shp_sql[i].find('ADJ1_POINTS')!=-1:rf=rfs[0]
            else:rf=rfs[1]
            if rf=='GDA2020\n':prj.write("GEOGCS['GDA2020',DATUM['D_GDA_2020',SPHEROID['GRS_1980',6378137,298.257222101]],PRIMEM['Greenwich',0],UNIT['Degree',0.017453292519943295]]")
            else: prj.write("GEOGCS['GCS_GDA_1994',DATUM['D_GDA_1994',SPHEROID['GRS_1980',6378137,298.257222101]],PRIMEM['Greenwich',0],UNIT['Degree',0.017453292519943295]]")
    i=i+1
    
f1 = r'./shp Layers/Stations'; w1 = shp.Writer(f1)
w1.field('STATION','C'); w1.field('CONST', 'C');w1.field('DESC', 'C');
w1.field('LAT1', 'N', decimal=10); w1.field('LONG1', 'N', decimal = 10); w1.field('H_ORTHO1', 'N', decimal = 4); w1.field('h_ELLIPSE1', 'N', decimal = 4);
w1.field('LAT2', 'N', decimal=10); w1.field('LONG2', 'N', decimal = 10); w1.field('H_ORTHO2', 'N', decimal = 4); w1.field('h_ELLIPSE2', 'N', decimal = 4);
w1.field('1_Hz_PosU', 'N', decimal=4); w1.field('1_Vt_PosU', 'N', decimal=4); w1.field('1_S_maj', 'N', decimal=4); w1.field('1_S_min', decimal = 4); w1.field('1_Orien', 'N', decimal = 4);
w1.field('2_Hz_PosU', 'N', decimal=4); w1.field('2_Vt_PosU', 'N', decimal=4); w1.field('2_S_maj', 'N', decimal=4); w1.field('2_S_min', decimal = 4); w1.field('2_Orien', 'N', decimal = 4);
w1.field('Mag_Coord_Chg', 'N', decimal=4);w1.field('Az_Coord_Chg', 'N', decimal=4); w1.field('Mag_PU_Chg', 'N', decimal = 4); w1.field('%_Coord_v_PU', decimal = 4); w1.field('%_PU_v_PU', 'N', decimal = 4);
w1.field('Mag_Vt_Chg', 'N', decimal=4);w1.field('Mag_vPU_Chg', 'N', decimal = 4); w1.field('%_Vt_v_PU', decimal = 4)
qry=cursor.execute("SELECT ADJ2_POINTS.STATION, ADJ2_POINTS.CONST, ADJ2_POINTS.DESC, \
                   ADJ1_POINTS.LATITUDE AS LATITUDE1, ADJ1_POINTS.LONGITUDE AS LONGITUDE2, ADJ1_POINTS.H_ORTHO AS H_ORTHO1, ADJ1_POINTS.h_ELLIPSE AS h_ELLIPSE1, \
                   ADJ2_POINTS.LATITUDE AS LATITUDE2, ADJ2_POINTS.LONGITUDE AS LONGITUDE2, ADJ2_POINTS.H_ORTHO AS H_ORTHO2, ADJ2_POINTS.h_ELLIPSE AS h_ELLIPSE2, \
                   APU1_POINTS.Hz_PosU AS APU1_Hz_PosU, APU1_POINTS.Vt_PosU AS APU1_Vt_PosU, APU1_POINTS.Semi_major AS APU1_Semi_major, APU1_POINTS.Semi_minor AS APU1_Semi_minor, APU1_POINTS.Orientation AS APU1_Orientation, \
                   APU2_POINTS.Hz_PosU AS APU2_Hz_PosU, APU2_POINTS.Vt_PosU AS APU2_Vt_PosU, APU2_POINTS.Semi_major AS APU2_Semi_major, APU2_POINTS.Semi_minor AS APU2_Semi_minor, APU2_POINTS.Orientation AS APU2_Orientation \
                   FROM ((ADJ2_POINTS LEFT JOIN APU1_POINTS ON ADJ2_POINTS.STATION = APU1_POINTS.STATION) LEFT JOIN APU2_POINTS ON ADJ2_POINTS.STATION = APU2_POINTS.STATION) LEFT JOIN ADJ1_POINTS ON ADJ2_POINTS.STATION = ADJ1_POINTS.STATION;").fetchall()
for row in qry:
    w1.point(row[8], row[7])
    Mag_Coord_Chg=[None,None]
    Mag_PU_Chg=None
    PER_Coord_v_PU=None
    PER_PU_v_PU=None
    Mag_Vt_Chg=None
    Mag_vPU_Chg=None
    PER_Vt_v_vPU=None
    if row[7]is not None and row[3]is not None and row[10]is not None and row[6]is not None:
        Mag_Coord_Chg=vincenty_inverse((row[7], row[8]),(row[3], row[4]))
        Mag_Vt_Chg=row[10]-row[6]
    if row[11]is not None and row[16]is not None and row[11]!=0:
        Mag_PU_Chg=row[16]-row[11]
        PER_Coord_v_PU=round((Mag_Coord_Chg[0]/row[11])*100,2)
        PER_PU_v_PU=round((Mag_PU_Chg/row[11])*100,2)
        Mag_Vt_Chg=row[10]-row[6]
        Mag_vPU_Chg=row[17]-row[12]
        PER_Vt_v_vPU=(Mag_Vt_Chg/row[12])*100
    w1.record(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], \
              row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20], \
              Mag_Coord_Chg[0], Mag_Coord_Chg[1], Mag_PU_Chg, PER_Coord_v_PU, PER_PU_v_PU,Mag_Vt_Chg,Mag_vPU_Chg,PER_Vt_v_vPU)
w1.close()
conn.close()
with open(f1 +'.prj', 'w') as prj:
    rf=rfs[1]
    if rf=='GDA2020\n':prj.write("GEOGCS['GDA2020',DATUM['D_GDA_2020',SPHEROID['GRS_1980',6378137,298.257222101]],PRIMEM['Greenwich',0],UNIT['Degree',0.017453292519943295]]")
    else: prj.write("GEOGCS['GCS_GDA_1994',DATUM['D_GDA_1994',SPHEROID['GRS_1980',6378137,298.257222101]],PRIMEM['Greenwich',0],UNIT['Degree',0.017453292519943295]]")
os.remove('Compare_Adj.db')
print ('...............................')
print( 'DONE')
print ('...............................')
