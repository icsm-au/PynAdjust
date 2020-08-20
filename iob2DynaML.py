# ----------------------------------------------------------------------
#                          iob2DynaML.py
# ----------------------------------------------------------------------
#  Author:  Kent Wheeler
#    Date:  17 March 2020
# Purpose:  Script to create DynaML stn and msr files from Geolab .iob file
# ----------------------------------------------------------------------
#   Usage:  cmd:\> python iob2DynaML.py -i <*.iob>
# ----------------------------------------------------------------------
#   Notes:  - If the input file contains the word final and the input contains
#             error ellipse information for fixed marks an extra stn and msr
#             file will be produced for doing a weighted adjustments
#             and 2-D propagation of uncertainty.
#           - If the iob contains start/finish time information of the
#             baselines, trivial baselines will be identified and exported
#             to kml.
# ----------------------------------------------------------------------
#  Update:  03 June 2020
#           - Refactor for Python 3 compliance, PEP8 styling
#           - Incorporate GeodePy functions as dependencies
# ----------------------------------------------------------------------

import os
import shutil, subprocess
import datetime as dt
import numpy as np
from math import sqrt, sin, cos, tan, radians, degrees
from geodepy.convert import hp2dec, llh2xyz, xyz2llh, dec2hp
from geodepy.geodesy import vincinv
from geodepy.constants import grs80
from datetime import datetime

class DnaStation:
    def __init__(self):
        self.name = ''
        self.Const = ''
        self.Type = ''
        self.XAxis = ''
        self.YAxis = ''
        self.Height = ''
        self.Desc = ''
        self.aAxis = 0
        self.bAxis = 0
        self.ErrAz = 0


class AdditionalInfoStn:
    def __init__(self):
        self.Ges_Name = ''
        self.Ges_Lat = ''
        self.Ges_Lng = ''
        self.Ges_Ht = ''
        self.Ges_Ht_Acc = ''
        self.Ges_Ht_Mth = ''
        self.Ges_Hz_Ord = ''
        self.Hz_Coord_Mth = ''
        self.rel_hz_acc = ''
        self.CE = ''
        self.NonGSNumber = ''
        self.SelectPoint = 'true'
        self.SelectRL = 'true'


class AdditionalInfoMsr:
    def __init__(self):
        # Additional information is included as a comment in the DynaML file.
        # This can be used for database import
        self.StartDateTime = dt.datetime(1994, 1, 1, 00, 00, 00)
        self.Duration = dt.datetime(1966, 1, 1, 00, 00, 00)
        self.TimeStatus = ''
        self.EphemerisType = ''
        self.AtReceiver = ''
        self.ToReceiver = ''
        self.FrequencyMode = ''
        self.SurveyTechnique = 'SLEV'
        self.Solution = ''
        self.EpochInterval = ''
        self.Class = 'LC'
        self.LevelDistance = '0.01'
        self.InstrumentModel = ''
        self.Derivation = 'MEAS'
        self.NonGSNumber = ''


class DnaMeasurement:
    def __init__(self):
        self.type = ''
        self.vscale = '1'
        self.pscale = '1'
        self.lscale = '1'
        self.hscale = '1'
        self.epoch = ''
        self.first = ''
        self.second = ''
        self.third = ''
        self.stddev = ''
        self.total = ''
        self.instheight = 0
        self.targheight = 0
        self.targets = []
        self.values = []
        self.targetstddevs = []
        self.dx = ''
        self.dy = ''
        self.dz = ''
        self.id = 0
        self.Vs = np.zeros([3, 3])

    def add_targets(self, target):
        self.targets.append(target)

    def add_values(self, value):
        self.values.append(value)

    def add_targetstddevs(self, targetstddev):
        self.targetstddevs.append(targetstddev)


class DeviceHeight:
    def __init__(self):
        # Device Height might be the height of instrument at a point
        # or Height of target
        self.StnName = []
        self.RefHeight = []


def add_device_height(self, stn, hgt):
    self.StnName.append(stn)
    self.RefHeight.append(hgt)


def hms2hp(hms_ang):
    # Input: HH MM SS.ssss used by Geolab
    # Output: HH.MMSSSsssss used by DynAdjust
    sign = 1
    if hms_ang.find('S') != -1 or hms_ang.find('-') != -1:
        sign = -1
    while hms_ang.find('  ') != -1:
        hms_ang = hms_ang.replace('  ', ' ')
    hms_ang = hms_ang.replace('S', '')
    hms_ang = hms_ang.replace('E', '')
    hms_ang = hms_ang.replace('.', ' ')
    ang_part = hms_ang.split()
    ang_part[0] = str(sign * abs(int(ang_part[0])))
    ang_part[1] = "%02d" % float(ang_part[1])
    ang_part[2] = "%02d" % float(ang_part[2])
    return ang_part[0] + '.' + ang_part[1] + ang_part[2] + ang_part[3]


def dd2hms(dd_ang,ds):
    #Input: DD.ddddd used by Geolab
    #Output: HH MM SS.sssss used by DynAdjust
    hp_ang = dec2hp(dd_ang)
    h=int(hp_ang)
    dec_m = (hp_ang-h)*100
    m=int(dec_m)
    s=(dec_m-m)*100
    return '{:} {:>2} {:>3}'.format(h, m, round(s,ds))


def find_job_num(strg):
    # search a string for 8 consecutive numbers, this is probably the Job Number
    job_num = ''
    i = 0
    while i+7 != len(strg):
        if strg[i:i + 8].isnumeric():
            job_num = strg[i:i + 8]
        i = i+1
    return job_num


def stn_xml_str(stn, stn_rec):
    # Output: String for printing to xml that is one complete station
    xml_str = '''<DnaStation>\n
    <Name>' + stn.Name + '</Name>\n
    <Constraints>' + stn.Const + '</Constraints>\n
    <Type>' + stn.Type + '</Type>\n
    <StationCoord>\n
    <Name>' + stn.Name + '</Name>\n
    <XAxis>' + str(stn.XAxis) + '</XAxis>\n
    <YAxis>' + str(stn.YAxis) + '</YAxis>\n
    <Height>' + stn.Height + '</Height>\n
    </StationCoord>\n
    <Description>' + stn.Desc + '</Description>\n
    <!--AdditionalInfoStn>\n
    <HorizCoordMethod>' + stn_rec.Hz_Coord_Mth + '</HorizCoordMethod>\n
    <RelativeHorizAccuracy>' + stn_rec.rel_hz_acc + '</RelativeHorizAccuracy>\n
    <NonGSNumber>' + stn_rec.NonGSNumber + '</NonGSNumber>\n
    <SelectPoint>' + stn_rec.SelectPoint + '</SelectPoint>\n
    <SelectRL>' + stn_rec.SelectRL + '</SelectRL>\n
    </AdditionalInfoStn-->\n
    </DnaStation>\n'''
    return xml_str


def msr_xml_str(msr, cntrl):
    # Output: xml string for printing to file.
    # Caters for type G, D, S, B, D, L, H
    xml_str = '<DnaMeasurement>\n'+ \
    '<Type>' + msr.type + '</Type>\n'+ \
    '<Ignore/>\n'
    if msr.type == 'G':
        xml_str = (xml_str + '<ReferenceFrame>' + 
                   gnss_date_2_ref(cntrl.StartDateTime) +
                   '</ReferenceFrame>\n'+ 
                   '<Epoch>' + 
                   cntrl.StartDateTime.strftime('%d.%m.%Y') + 
                   '</Epoch>\n' + 
                   '<Vscale>{0:.1f}</Vscale>\n'.format(float(msr.vscale)) 
                   + '<Pscale>{0:.1f}</Pscale>\n'.format(float(msr.pscale)) 
                   + '<Lscale>{0:.1f}</Lscale>\n'.format(float(msr.lscale)) 
                   + '<Hscale>{0:.1f}</Hscale>\n'.format(float(msr.hscale)))
    xml_str = xml_str + '<First>' + msr.first + '</First>\n'
    if msr.second != '':
        xml_str = xml_str + '<Second>' + msr.second + '</Second>\n'
    if msr.type != 'G' and msr.type != 'D':
        xml_str = (xml_str + '<Value>' + msr.values[0] + '</Value>\n' + 
             '<StdDev>' + msr.stddev + '</StdDev>\n')
    if msr.type == 'G':
        xml_str = (xml_str + '<GPSBaseline>\n'+
        '<X>' + str(msr.dx) + '</X>\n'+
        '<Y>' + str(msr.dy) + '</Y>\n'+
        '<Z>' + str(msr.dz) + '</Z>\n'+
        '<MeasurementID>' + str(msr.id) + '</MeasurementID>\n'+
        '<SigmaXX>' + str(msr.Vs[0, 0]) + '</SigmaXX>\n'+
        '<SigmaXY>' + str(msr.Vs[0, 1]) + '</SigmaXY>\n'+
        '<SigmaXZ>' + str(msr.Vs[0, 2]) + '</SigmaXZ>\n'+
        '<SigmaYY>' + str(msr.Vs[1, 1]) + '</SigmaYY>\n'+
        '<SigmaYZ>' + str(msr.Vs[1, 2]) + '</SigmaYZ>\n'+
        '<SigmaZZ>' + str(msr.Vs[2,2]) + '</SigmaZZ>\n'+
        '</GPSBaseline>\n '+
        '<!--AdditionalInfoMsrG>\n'+ 
        '<StartDateTime>' + 
        cntrl.StartDateTime.strftime('%Y-%m-%dT%H:%M:%S') + 
        '</StartDateTime>\n'+ 
        '<Duration>' + 
        'P' + str(cntrl.Duration.year - 1900) + 
        'Y' + str(cntrl.Duration.month - 1) + 
        'M' + str(cntrl.Duration.day - 1) + 
        'DT' + cntrl.Duration.strftime('%HH%MM%SS') + 
        '</Duration>\n'+ 
        '<TimeStatus>' + cntrl.TimeStatus + '</TimeStatus>\n'+ 
        '<EphemerisType>' + cntrl.EphemerisType + '</EphemerisType>\n'+ 
        '<AtReceiver>' + cntrl.AtReceiver + '</AtReceiver>\n'+ 
        '<ToReceiver>' + cntrl.ToReceiver + '</ToReceiver>\n'+ 
        '<FrequencyMode>' + cntrl.FrequencyMode + '</FrequencyMode>\n'+ 
        '<SurveyTechnique>' + cntrl.SurveyTechnique + '</SurveyTechnique>\n'+ 
        '<Solution>' + cntrl.Solution + '</Solution>\n'+ 
        '<EpochInterval>' + str(cntrl.EpochInterval) + '</EpochInterval>\n'+ 
        '<Class>' + cntrl.Class + '</Class>\n'+ 
        '<NonGSNumber>' + str(cntrl.NonGSNumber) + '</NonGSNumber>\n'+ 
        '</AdditionalInfoMsrG-->\n')
    
    if msr.type == 'L':
        xml_str = (xml_str + '<!--AdditionalInfoMsrL>\n'+ 
        '<SurveyTechnique>' + cntrl.SurveyTechnique + '</SurveyTechnique>\n'+ 
        '<LevelDistance>' + cntrl.LevelDistance + '</LevelDistance>\n'+ 
        '<ObsDate>' + cntrl.StartDateTime.strftime('%Y-%m-%d') + '</ObsDate>\n'+ 
        '<Derivation>' + cntrl.Derivation + '</Derivation>\n'+ 
        '<Class>' + cntrl.Class + '</Class>\n'+ 
        '<NonGSNumber>' + cntrl.NonGSNumber + '</NonGSNumber>\n'+ 
        '</AdditionalInfoMsrL-->\n')
    
    if msr.type == 'S':
        xml_str = (xml_str + '<InstHeight>' + str(msr.instheight) + '</InstHeight>\n'+ 
        '<TargHeight>' + str(msr.targheight) + '</TargHeight>\n'+ 
        '<!--AdditionalInfoMsrS>\n'+ 
        '<InstrumentModel>' + cntrl.InstrumentModel + '</InstrumentModel>\n'+ 
        '<ObsDate>' + cntrl.StartDateTime.strftime('%Y-%m-%d') + '</ObsDate>\n'+ 
        '<Derivation>' + cntrl.Derivation + '</Derivation>\n'+ 
        '<Class>' + cntrl.Class + '</Class>\n'+ 
        '<NonGSNumber>' + cntrl.NonGSNumber + '</NonGSNumber>\n'+ 
        '</AdditionalInfoMsrS-->\n')

    if msr.type == 'D':
        xml_str = (xml_str + '<Value>' + msr.values[0] + '</Value>\n'+ 
        '<StdDev>' + msr.targetstddevs[0] + '</StdDev>\n'+ 
        '<Total>' + str(msr.total - 1) + '</Total>\n')
        obs_num = 1
        while obs_num < msr.total:
            xml_str = (xml_str + '<Directions>\n'+ 
            '<Ignore/>\n'+ 
            '<Target>' + msr.targets[obs_num] + '</Target>\n'+ 
            '<Value>' + msr.values[obs_num] + '</Value>\n'+ 
            '<StdDev>' + msr.targetstddevs[obs_num] + '</StdDev>\n'+ 
            '</Directions>\n')
            obs_num = obs_num + 1
        xml_str = (xml_str + '<!--AdditionalInfoMsrD>\n'+ 
        '<InstrumentModel>' + cntrl.InstrumentModel + '</InstrumentModel>\n'+ 
        '<ObsDate>' + cntrl.StartDateTime.strftime('%Y-%m-%d') + '</ObsDate>\n'+ 
        '<Derivation>' + cntrl.Derivation + '</Derivation>\n'+ 
        '<Class>' + cntrl.Class + '</Class>\n'+ 
        '<NonGSNumber>' + cntrl.NonGSNumber + '</NonGSNumber>\n'+ 
        '</AdditionalInfoMsrD-->\n')
    xml_str = xml_str+'<Source></Source>\n'
    if msr.type != 'G':
        xml_str = xml_str+'<MeasurementID>' + str(msr.id) + '</MeasurementID>\n'
    xml_str = xml_str+'</DnaMeasurement>\n'
    
    return xml_str


def get_rl(p,l):
    rlat = radians(p)
    rlng = radians(l)
    rl = np.zeros([3, 3])
    rl[0, 0] = -sin(rlng)
    rl[0, 1] = -sin(rlat)*cos(rlng)
    rl[0, 2] = cos(rlat)*cos(rlng)
    rl[1, 0] = cos(rlng)
    rl[1, 1] = -sin(rlat)*sin(rlng)
    rl[1, 2] = cos(rlat)*sin(rlng)
    rl[2, 1] = cos(rlat)
    rl[2, 1] = sin(rlat)
    return rl


def err_ellip_2_ycluster(stn):
    # Input: Supply a station with coordinates
    #        and error ellipse for coordinate uncertainty
    # Output: xml string for  point cluster (Y-type observation)
    x, y, z = llh2xyz(float(stn.XAxis), float(stn.YAxis), float(stn.Height))
    
    a = stn.aAxis / 2.44774683068
    b = stn.bAxis / 2.44774683068
    az = 90 - stn.ErrAz
    
    r_az = radians(az)
    rl = get_rl(float(stn.XAxis),float(stn.YAxis))
    
    i_a = np.zeros([3, 3])
    i_a[0, 0] = (cos(r_az)*cos(r_az)*a*a)+(b*b*sin(r_az)*sin(r_az))
    i_a[0, 1] = (a*a-b*b)*cos(r_az)*sin(r_az)
    i_a[1, 0] = i_a[0, 1]
    i_a[1, 1] = (a*a*sin(r_az)*sin(r_az))+(b*b*cos(r_az)*cos(r_az))
    i_a[2, 2] = 0.000001
    
    wt = np.matmul(np.matmul(rl, i_a), rl.transpose())
    
    xml_str = ('<DnaMeasurement>\n'+ 
        '<Type>Y</Type>\n'+ 
        '<Ignore/>\n'+ 
        '<ReferenceFrame>GDA2020</ReferenceFrame>\n'+ 
        '<Epoch>01.01.2020</Epoch>\n'+ 
        '<Vscale>1.000</Vscale>\n'+ 
        '<Pscale>1.000</Pscale>\n'+ 
        '<Lscale>1.000</Lscale>\n'+ 
        '<Hscale>1.000</Hscale>\n'+ 
        '<Coords>XYZ</Coords>\n'+ 
        '<Total>1</Total>\n'+ 
        '<First>' + stn.Name + '</First>\n'+ 
        '<Clusterpoint>\n'+ 
        '<X>'+str(x)+'</X>\n'+ 
        '<Y>'+str(y)+'</Y>\n'+ 
        '<Z>'+str(z)+'</Z>\n'+ 
        '<SigmaXX>'+str(wt[0, 0])+'</SigmaXX>\n'+ 
        '<SigmaXY>'+str(wt[0, 1])+'</SigmaXY>\n'+ 
        '<SigmaXZ>'+str(wt[0, 2])+'</SigmaXZ>\n'+ 
        '<SigmaYY>'+str(wt[1, 1])+'</SigmaYY>\n'+ 
        '<SigmaYZ>'+str(wt[1, 2])+'</SigmaYZ>\n'+ 
        '<SigmaZZ>'+str(wt[2, 2])+'</SigmaZZ>\n'+ 
        '</Clusterpoint>\n'+ 
        '</DnaMeasurement>\n')
    
    return xml_str


def stn_header(datum):
    xml_str='<?xml version="1.0" encoding="utf-8"?>\n'
    epoch='01.01.2020'
    if datum.find('09')!=-1: epoch='01.01.1994'
    xml_str=(xml_str+'<DnaXmlFormat type="Station File" referenceframe="'+datum+
             '" epoch="'+epoch+
             '" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'+
             'xsi:noNamespaceSchemaLocation="DynaML.xsd">\n')
    return xml_str


def msr_header():
    xml_str = '<?xml version="1.0" encoding="utf-8"?>\n'+\
               '<DnaXmlFormat type="Measurement File" '+\
               'referenceframe="GDA94" epoch="01.01.1994" '+\
               'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '+\
               'xsi:noNamespaceSchemaLocation="DynaML.xsd">\n'
    return xml_str


def dml_footer():
    xml_str = '</DnaXmlFormat>\n'
    return xml_str


def gnss_date_2_ref(obs_date):
    # Use the date of GNSS baseline observation
    # to determine the reference frame used by broadcast ephemeris
    if dt.datetime(1900, 1, 1) <= obs_date < dt.datetime(1994, 1, 2):
        ref = 'ITRF1991'
    if dt.datetime(1994, 1, 2) <= obs_date < dt.datetime(1995, 1, 1):
        ref = 'ITRF1992'
    if dt.datetime(1995, 1, 1) <= obs_date < dt.datetime(1996, 6, 30):
        ref = 'ITRF1993'
    if dt.datetime(1996, 6, 30) <= obs_date < dt.datetime(1998, 3, 1):
        ref = 'ITRF1994'
    if dt.datetime(1998, 3, 1) <= obs_date < dt.datetime(1999, 8, 1):
        ref = 'ITRF1996'
    if dt.datetime(1999, 8, 1) <= obs_date < dt.datetime(2001, 12, 2):
        ref = 'ITRF1997'
    if dt.datetime(2001, 12, 2) <= obs_date < dt.datetime(2006, 11, 5):
        ref = 'ITRF2000'
    if dt.datetime(2006, 11, 5) <= obs_date < dt.datetime(2011, 4, 17):
        ref = 'ITRF2005'
    if dt.datetime(2011, 4, 17) <= obs_date < dt.datetime(2017, 1, 29):
        ref = 'ITRF2008'
    if obs_date >= dt.datetime(2017, 1, 29):
        ref = 'ITRF2014'
    return ref


def DynaGeol(adj_file):
    with open(adj_file) as f:
        first_line = f.readline()
    if first_line=='This DynAdjust file has been re-formatted by program DynaGEOL2020.\n':
        print(adj_file +' has already been re-formatted by program DynaGEOL2020.\n')
    else:
        try: os.remove('$'+adj_file)
        except: ''
        os.rename( adj_file, '$' + adj_file)
        f_in = open('$'+adj_file, 'r')
        f_out = open(adj_file, 'w')
        part_G_cnt=1; base_cnt=0; coords = {}; xyz_bases = {}; dah_bases = {}
        read_msr=False; read_coord=-1; angular_msr_type=''
        f_out.write('This GEOLAB file has been re-formatted by program DynaGEOL2020.\n')
        f_out.write('DYNADJUST ADJUSTMENT OUTPUT FILE\n')
        for linestr in f_in.readlines():
            f_out.write(linestr)
            if (linestr.find('Command line arguments:')!=-1 and 
                linestr.find('--angular-msr-type 1')!=-1):angular_msr_type='dd'
            if linestr.find('Station coordinate types:')!=-1: coord_types = linestr[35:].strip()
            if linestr.find('Adjusted Measurements')!=-1: read_msr=True
            if read_msr and linestr[0:2]=='G ':
                a_obs=linestr[67:].split()
                if part_G_cnt==1: dX=float(a_obs[0]); sdev_x=float(a_obs[3]); nstat_x=float(a_obs[6])
                if part_G_cnt==2: dY=float(a_obs[0]); sdev_y=float(a_obs[3]); nstat_y=float(a_obs[6])
                if part_G_cnt==3: 
                    dZ=float(a_obs[0]); sdev_z=float(a_obs[3]); nstat_z=float(a_obs[6])
                    xyz_bases[base_cnt] = {'first':linestr[2:22].strip(),'second':linestr[22:42].strip(),
                             'dX':dX, 'dY':dY, 'dZ':dZ,
                             'sdev_x':sdev_x, 'sdev_y':sdev_y, 'sdev_z':sdev_z,
                             'nstat_x':nstat_x, 'nstat_y':nstat_y, 'nstat_z':nstat_z}
                    f_out.write('\n')
                    base_cnt+=1
                    part_G_cnt=0
                part_G_cnt+=1
                
            if linestr.find('Adjusted Coordinates')!=-1:
                read_msr=False; read_coord=0
            if read_coord>4 and linestr!='\n':
                stn = linestr[0:20].strip()
                results = linestr[25:].split()
                r_count = 0
                for ct in coord_types:
                    if ct == 'E': E = float(results[r_count])
                    if ct == 'N': N = float(results[r_count])
                    if ct == 'z': z = int(results[r_count])
                    if ct == 'P':
                        if angular_msr_type!='dd':
                            P = hp2dec(float(results[r_count]))
                        else:
                            P = float(results[r_count])
                    if ct == 'L':
                        if angular_msr_type!='dd':
                            L = hp2dec(float(results[r_count]))
                        else:
                            L = float(results[r_count])
                    if ct == 'H': H = float(results[r_count])
                    if ct == 'h': h = float(results[r_count])
                    r_count += 1
                        
                coords[str(stn)] = {'E': E,'N': N,'Z': z,'LAT': P,'LON': L,'OHGT': H,'EHGT': h }
            if read_coord>=0:read_coord+=1
        f_in.close()
        f_out.close()
        os.remove('$'+adj_file)
    print('   Calculating baseline transformations ')
    for b in xyz_bases:
        f_x, f_y, f_z = llh2xyz(coords[xyz_bases[b]['first']]['LAT'],coords[xyz_bases[b]['first']]['LON'],coords[xyz_bases[b]['first']]['EHGT'])
        Calc_point_x = f_x+xyz_bases[b]['dX']
        Calc_point_y = f_y+xyz_bases[b]['dY']
        Calc_point_z = f_z+xyz_bases[b]['dZ']
        p, l, h = xyz2llh(Calc_point_x, Calc_point_y, Calc_point_z)
        e_dist, b_az, b_rev_az = vincinv(coords[xyz_bases[b]['first']]['LAT'],coords[xyz_bases[b]['first']]['LON'],p, l)
        adj_dist, adj_az, b_rev_az = vincinv(coords[xyz_bases[b]['first']]['LAT'],coords[xyz_bases[b]['first']]['LON'],coords[xyz_bases[b]['second']]['LAT'],coords[xyz_bases[b]['second']]['LON'])
        e_dist_res = adj_dist - e_dist
        b_az_res = adj_az - b_az
        e_ht_diff = h - coords[xyz_bases[b]['first']]['EHGT']
        n1=coords[xyz_bases[b]['first']]['OHGT']-coords[xyz_bases[b]['first']]['EHGT']
        n2=coords[xyz_bases[b]['second']]['OHGT']-coords[xyz_bases[b]['second']]['EHGT']
        o_ht_diff = -1*n1 + e_ht_diff + n2
        o_ht_diff_res =(coords[xyz_bases[b]['second']]['OHGT']-coords[xyz_bases[b]['first']]['OHGT']) - o_ht_diff
        
        #transfor the baseline to distance, azimuth, up
        #translated code from http://peas13.dli.wa.gov.au/svn-intsvc/Delphi/GesmarComponents/trunk/Geodetic.pas
        m_lat = (coords[xyz_bases[b]['first']]['LAT']+p)/2
        m_x =(f_x+Calc_point_x)/2
        m_y =(f_y+Calc_point_y)/2
        
        dR       = sqrt(m_x*m_x + m_y*m_y)
        dE2      = grs80.ecc1sq
        dQ       = radians(m_lat)
        sinQ     = sin(dQ)
        dTemp    = sqrt(1.0 - dE2 * sinQ * sinQ)
        dN       = grs80.semimaj / dTemp
        dM       = dN * ((1.0 - dE2) / (dTemp * dTemp))
        sinQ     = sin(dQ)
        cosQ     = cos(dQ)
        secQ     = 1.0/cosQ
        dF       = grs80.f
        dE2      = 2 * dF -(dF * dF)
        dA       = m_x * tan(dQ)/(dR*dR*(dE2 - secQ*secQ))
        dB       = 1.0/(dR*secQ*secQ - dE2*dN*cosQ)
        dC       = dR*sinQ/(cosQ*cosQ)-dN*dE2*sinQ*cosQ/(1.0-dE2*sinQ*sinQ)
        
        JMatrix      = np.zeros([3,3])
        JMatrix[0,0] = dM*dA
        JMatrix[0,1] = dM*(m_y*dA/m_x)
        JMatrix[0,2] = dM*dB
        JMatrix[1,0] = -(dN*cosQ)*m_y/(dR*dR)
        JMatrix[1,1] = (dN*cosQ)*m_x/(dR*dR)
        JMatrix[1,2] = 0.0
        JMatrix[2,0] = m_x/(dR*cosQ) + dA*dC
        JMatrix[2,1] = m_y/(dR*cosQ) + m_y*dA*dC/m_x
        JMatrix[2,2] = dB*dC
        
        b_var_matrix = np.zeros([3,3])
        b_var_matrix[0,0] = xyz_bases[b]['sdev_x']**2
        b_var_matrix[1,1] = xyz_bases[b]['sdev_y']**2
        b_var_matrix[2,2] = xyz_bases[b]['sdev_z']**2
                
        b_nst_matrix = np.zeros([3,3])
        b_nst_matrix[0,0] = xyz_bases[b]['nstat_x']**2
        b_nst_matrix[1,1] = xyz_bases[b]['nstat_y']**2
        b_nst_matrix[2,2] = xyz_bases[b]['nstat_z']**2

        cosAz        = cos(radians(b_az))
        sinAz        = sin(radians(b_az))
        AMatrix      = np.zeros([3,3])
        AMatrix[0,0] = cosAz
        AMatrix[0,1] = sinAz
        AMatrix[1,0] = degrees(sinAz/e_dist)
        AMatrix[1,1] = -degrees(cosAz/e_dist)

        GMatrix = np.matmul(np.matmul(JMatrix,b_var_matrix),JMatrix.transpose())
        dah_Matrix = np.matmul(np.matmul(AMatrix,GMatrix),AMatrix.transpose())

        n_GMatrix = np.matmul(np.matmul(JMatrix,b_nst_matrix),JMatrix.transpose())
        nstat_Matrix = np.matmul(np.matmul(AMatrix,n_GMatrix),AMatrix.transpose())

        dah_bases[b] = {'first':xyz_bases[b]['first'],'second':xyz_bases[b]['second'],
                        'dX':xyz_bases[b]['dX'], 'dY':xyz_bases[b]['dY'], 'dZ':xyz_bases[b]['dZ'],
                        'e_dist':e_dist, 'e_dist_res':e_dist_res,'e_dist_sdev':sqrt(abs(dah_Matrix[0,0])),'e_dist_nstat':sqrt(abs(nstat_Matrix[0,0])),
                        'b_az':b_az, 'b_az_res':b_az_res, 'b_az_sdev':sqrt(abs(dah_Matrix[1,1])),'b_az_nstat':sqrt(sin(radians(nstat_Matrix[1,1]))*e_dist),
                        'o_ht_diff':o_ht_diff, 'e_ht_diff':e_ht_diff, 'o_ht_diff_res':o_ht_diff_res,'e_ht_diff_sdev':sqrt(abs(GMatrix[2,2])),'e_ht_diff_nstat':sqrt(abs(n_GMatrix[2,2]))}
    return coords, dah_bases


def Create_DynaBAS(adj_file,bases):
    print('   writing: '+adj_file[:-3]+'BAS')
    f_bas = open(adj_file[:-3]+'BAS', 'w')
    f_bas.write( DynaBAS_header(adj_file))
    f_bas.write('                                       Ellipsoidal\n')
    f_bas.write('BASE                                   Distance                Standard\n')
    f_bas.write('LINE        From         To            Observation   Residual  Deviation     N-Stat\n')
    f_bas.write('----------- ------------ ------------  -----------  ---------  ---------   --------\n')
    for b in bases:
        f_bas.write('{:11s} {:12s} {:12s} {:>12} {:>10} {:>10} {:>10}'.format(
                str(b), bases[b]['first'], bases[b]['second'],
                '{:.3f}'.format(bases[b]['e_dist']), '{:.3f}'.format(bases[b]['e_dist_res']),
                '{:.3f}'.format(bases[b]['e_dist_sdev']), '{:.2f}'.format(bases[b]['e_dist_nstat'])) + '\n')
    
    f_bas.write('\n')
    f_bas.write( DynaBAS_header(adj_file))
    f_bas.write('\n')
    f_bas.write('BASE                                   Azimuth                  Standard\n')
    f_bas.write('LINE        From         To            Observation    Residual  Deviation     N-Stat\n')
    f_bas.write('----------- ------------ ------------  ------------  ---------  ---------   --------\n')
    for b in bases:
        f_bas.write('{:11s} {:12s} {:12s} {:>13s} {:>10} {:>10} {:>10}'.format(
                str(b), bases[b]['first'], bases[b]['second'], 
                dd2hms(bases[b]['b_az'],2), '{:.2f}'.format(bases[b]['b_az_res']*3600),
                '{:.2f}'.format(bases[b]['b_az_sdev']*3600),'{:.2f}'.format(bases[b]['b_az_nstat'])) + '\n')
    
    f_bas.write('\n')
    f_bas.write( DynaBAS_header(adj_file))
    f_bas.write('                                       Height Difference\n')
    f_bas.write('                                       Observation\n')
    f_bas.write('BASE                                   -----------------------             Standard\n')
    f_bas.write('LINE        From         To            Orthometric  Spheroidal   Residual  Deviation     N-Stat\n')
    f_bas.write('----------- ------------ ------------  ----------- -----------   --------  ---------   --------\n')
    for b in bases:
        f_bas.write('{:11s} {:12s} {:12s} {:>12} {:>11} {:>10} {:>10} {:>10}'.format(
                str(b), bases[b]['first'], bases[b]['second'], '{:.3f}'.format(
                        bases[b]['o_ht_diff']), '{:.3f}'.format(bases[b]['e_ht_diff']), '{:.3f}'.format(bases[b]['o_ht_diff_res']),
                        '{:.3f}'.format(bases[b]['e_ht_diff_sdev']),'{:.2f}'.format(bases[b]['e_ht_diff_nstat'])) + '\n')
    
    f_bas.close()


def DynaBAS_header(adj_file):
    header_str='================================================================================\n'
    header_str=header_str + adj_file.center(80)+'\n'
    header_str=header_str + 'GRS80\n'
    header_str=header_str + '================================================================================\n'
    return header_str


if __name__ == "__main__":
    rel_uncert_sum = 0
    avg_rel_uncert = 0
    stn_cnt = 0
    xsd_dir='Z:\AssetsPCsAndSoftware\Software\DynAdjust'
    xsd_file = os.path.join(xsd_dir, 'DynaML.xsd')
    script_path = os.path.abspath(os.path.realpath(__file__))
    script_dir, script_name = os.path.split(script_path)
    os.chdir(script_dir)
    shutil.copy(xsd_file, os.getcwd())
    
    files = os.listdir(os.getcwd())
    for df in files:
        if df.endswith(".iob"):
            # Open, run through and close the file
            # for initial information on the adjustment
            print('Reading the Geolab (.iob) File...')
            adjustment_name=df.replace('.iob','')
            with open(df, 'r') as f:
                GNSSmarksStr = ';'
                FloatmarksStr = ';'
                w_MksWithObs = ';'
                for ln in f.readlines():
                    if ln[0:4] == ' PLO' or ln[0:4] == ' PLH':
                        if (ln[72:len(ln)].strip() != ''
                                and ln.find('ppm') != -1):
                            stnRec = ln[72:len(ln)].strip().replace(';', '|').split('|')
                            rel_uncert_sum = (rel_uncert_sum +
                                              float(stnRec[7].replace('ppm', '')))
                            stn_cnt = stn_cnt + 1
                    if ln[0:6]==' GFIL ': 
                        geoid=ln[6:].strip()
                        geoid_path='\\'.join(geoid.split('\\')[0:-1])
                        gfiles = os.listdir(geoid_path)
                        datum='GDA2020'
                        if geoid_path.find('09')!=-1:datum='GDA94'
                        for g in gfiles:
                            if g.endswith('.gsb'):geoid=geoid_path +'\\' + g
                    if (ln[0:5] == ' OHDF' or
                        ln[0:5] == ' GAZI' or
                        ln[0:5] == ' AZIM' or
                        ln[0:5] == ' DIST' or
                        ln[0:4] == ' DIR' or
                            ln[0:5] == ' DXYZ'):
                        CurrentMsr = DnaMeasurement()
                        CurrentMsr.first = ln[10:23].strip()
                        CurrentMsr.second = ln[23:35].strip()
                        if ln[0:5] == ' DXYZ':
                            GNSSmarksStr = (GNSSmarksStr + ';' +
                                            CurrentMsr.first + ';' + CurrentMsr.second)
                        if (FloatmarksStr.find(';' + CurrentMsr.first+';') != -1 or
                                FloatmarksStr.find(';'+CurrentMsr.second + ';') != -1):
                            w_MksWithObs = (w_MksWithObs + CurrentMsr.first + ';' +
                                            CurrentMsr.second + ';')
        
            # bin the guessing of GDA94 relative uncertainty for the float marks
            if stn_cnt != 0:
                avg_rel_uncert = rel_uncert_sum / stn_cnt
                if avg_rel_uncert < 3:
                    avg_rel_uncert = 3
                if 3 < avg_rel_uncert <= 7:
                    avg_rel_uncert = 7.5
                if 7.5 < avg_rel_uncert <= 10:
                    avg_rel_uncert = 10
                if 10 < avg_rel_uncert <= 20:
                    avg_rel_uncert = 20
                if 20 < avg_rel_uncert <= 30:
                    avg_rel_uncert = 30
                if 30 < avg_rel_uncert <= 50:
                    avg_rel_uncert = 50
                if avg_rel_uncert > 50:
                    avg_rel_uncert = int(avg_rel_uncert / 10) * 10
            with open(df, 'r') as f:
                stnout = open(df.replace('.iob', '.stn.xml'), 'w')
                msrout = open(df.replace('.iob', '.msr.xml'), 'w')
                stnout.write(stn_header(datum))
                msrout.write(msr_header())
        
                # Run through each line of the input file and extract the relevant lines
                lineCount = 0; msr_cnt=0
                InstHts = DeviceHeight()
                TgtHts = DeviceHeight()
                ControlRec = AdditionalInfoMsr()
                for ln in f.readlines():
                    if ln[0:5] == ' TITL':
                        jobNumber = find_job_num(ln)
                        if jobNumber == '':
                            jobNumber = find_job_num(os.getcwd())
                    if ln[0:4] == ' PLO' or ln[0:4] == ' PLH':
                        CurrentStn = DnaStation()
                        CurrentStn.Name = ln[10:23].strip()
                        CurrentStn.Const = ln[6:9].replace('1', 'C')
                        CurrentStn.Const = CurrentStn.Const.replace('0', 'F')
                        CurrentStn.Const = CurrentStn.Const.strip()
                        if ln[0:4] == ' PLO':
                            CurrentStn.Type = 'LLH'
                        if ln[0:4] == ' PLH':
                            CurrentStn.Type = 'LLh'
                        CurrentStn.XAxis = hms2hp(ln[23:41].strip())
                        CurrentStn.YAxis = hms2hp(ln[41:59].strip())
                        CurrentStn.Height = ln[59:72].strip()
                        CurrentStn.Desc = ln[72:len(ln)].strip()
                        stn_addn_rec = AdditionalInfoStn()
                        stn_addn_rec.NonGSNumber = 'E' + jobNumber
                        if CurrentStn.Desc != '':
                            stnRec = CurrentStn.Desc.replace(';', '|').split('|')
                            stn_addn_rec.Ges_Name = stnRec[0]
                            stn_addn_rec.Ges_Lat = stnRec[1]
                            stn_addn_rec.Ges_Lng = stnRec[2]
                            stn_addn_rec.Ges_Ht = stnRec[3]
                            stn_addn_rec.Ges_Ht_Acc = stnRec[4]
                            stn_addn_rec.Ges_Ht_Mth = stnRec[5]
                            stn_addn_rec.Ges_Hz_Ord = stnRec[6]
                            stn_addn_rec.rel_hz_acc = stnRec[7]
                            stn_addn_rec.Hz_Coord_Mth = stnRec[8]
                            if stnRec[10] != '':
                                CurrentStn.CE = stnRec[9]
                                CurrentStn.aAxis = float(stnRec[10])
                                CurrentStn.bAxis = float(stnRec[11])
                                CurrentStn.ErrAz = float(stnRec[12])
                            stn_addn_rec.SelectPoint = 'true'
                            stn_addn_rec.SelectRL = 'true'
        
                        if CurrentStn.Const[0:2] == 'FF' and avg_rel_uncert != 0:
                            stn_addn_rec.rel_hz_acc = str(avg_rel_uncert) + 'ppm'
                            if GNSSmarksStr.find(';'+CurrentStn.Name+';'):
                                stn_addn_rec.Hz_Coord_Mth = 'GNSS'
        
                        stnout.write(stn_xml_str(CurrentStn, stn_addn_rec))
        
                    if ln[0:5] == ' HI  ':
                        add_device_height(InstHts, ln[10:23].strip(), ln[23:33].strip())
        
                    if ln[0:5] == ' HT  ':
                        add_device_height(TgtHts, ln[10:23].strip(), ln[23:33].strip())
        
                    if (ln[0:5] == ' OHDF' or
                        ln[0:5] == ' OHGT' or
                        ln[0:5] == ' GAZI' or
                        ln[0:5] == ' AZIM' or
                        ln[0:5] == ' DIST' or
                        ln[0:5] == ' DSET' or
                        ln[0:5] == ' GRP '):
                            msr_cnt+=1
                            CurrentMsr = DnaMeasurement()
                            CurrentMsr.id = msr_cnt
                            
                    if ln[0:5] == ' OHGT':
                        CurrentMsr.type = 'H'
                        CurrentMsr.first = ln[10:23].strip()
                        CurrentMsr.add_values(ln[36:65].strip())
                        CurrentMsr.stddev = ln[65:76].strip()
                        msrout.write(msr_xml_str(CurrentMsr, ControlRec))
        
                    if ln[0:5] == ' OHDF':
                        CurrentMsr.type = 'L'
                        CurrentMsr.first = ln[10:23].strip()
                        CurrentMsr.second = ln[23:35].strip()
                        CurrentMsr.add_values(ln[50:65].strip())
                        CurrentMsr.stddev = ln[65:76].strip()
                        ControlRec.LevelDistance = ln[36:50].strip()
                        msrout.write(msr_xml_str(CurrentMsr, ControlRec))
        
                    if ln[0:5] == ' GAZI':
                        CurrentMsr.type = 'B'
                        CurrentMsr.first = ln[10:23].strip()
                        CurrentMsr.second = ln[23:35].strip()
                        CurrentMsr.add_values(hms2hp(ln[36:65].strip()))
                        CurrentMsr.stddev = ln[65:76].strip()
                        msrout.write(msr_xml_str(CurrentMsr, ControlRec))
        
                    if ln[0:5] == ' AZIM':
                        CurrentMsr.type = 'K'
                        CurrentMsr.first = ln[10:23].strip()
                        CurrentMsr.second = ln[23:35].strip()
                        CurrentMsr.add_values(hms2hp(ln[36:65].strip()))
                        CurrentMsr.stddev = ln[65:76].strip()
                        msrout.write(msr_xml_str(CurrentMsr, ControlRec))
        
                    if ln[0:5] == ' DIST':
                        CurrentMsr.type = 'S'
                        CurrentMsr.first = ln[10:23].strip()
                        CurrentMsr.second = ln[23:35].strip()
                        CurrentMsr.add_values(ln[36:65].strip())
                        rw = 0
                        for Stn in InstHts.StnName:
                            if Stn == CurrentMsr.first:
                                CurrentMsr.instheight = InstHts.RefHeight[rw]
                            rw = rw+1
                        rw = 0
                        for Stn in TgtHts.StnName:
                            if Stn == CurrentMsr.first:
                                CurrentMsr.targheight = TgtHts.RefHeight[rw]
                            rw = rw+1
                        CurrentMsr.stddev = ln[65:76].strip()
                        msrout.write(msr_xml_str(CurrentMsr, ControlRec))
                        
                    if ln[0:5] == ' DSET':
                        CurrentMsr.type = 'D'
                        lineCount = 0
                    if ln[0:4] == ' DIR' and lineCount == 1 and CurrentMsr.type == 'D':
                        CurrentMsr.first = ln[10:23].strip()
                        CurrentMsr.second = ln[23:35].strip()
                        CurrentMsr.add_targets(ln[23:35].strip())
                        CurrentMsr.add_values(hms2hp(ln[36:65].strip()))
                        CurrentMsr.add_targetstddevs(ln[65:76].strip())
                        CurrentMsr.total = lineCount
                    if ln[0:4] == ' DIR' and lineCount > 1 and CurrentMsr.type == 'D':
                        CurrentMsr.add_targets(ln[23:35].strip())
                        CurrentMsr.add_values(hms2hp(ln[36:65].strip()))
                        CurrentMsr.add_targetstddevs(ln[65:76].strip())
                        CurrentMsr.total = lineCount
                    if CurrentMsr.type == 'D' and ln[0:4] != ' DIR' and lineCount > 1:
                        msrout.write(msr_xml_str(CurrentMsr, ControlRec))
                        CurrentMsr = DnaMeasurement()
        
                # Scrape information from Landgate GESMAR control Records
                # eg.*CONTROL;GPS;201810010007;012057;E;B;TRIM;TRIM;D;ST  ;FX;015;N
                # eg.*CONTROL;OHDF;20181001;SLEV;LC;2.71;MEAS
                # eg.*CONTROL;DIS;20181213;TS 16;C;MEAS
                # eg.*CONTROL;ANG;20181213;TS16;C;MEAS
                    if ln[0:9] == '*CONTROL;':
                        ControlRec = AdditionalInfoMsr()
                        alinestr = ln.split(';')
                        stdatetimestr = alinestr[2] + '0000'
                        yr = int(stdatetimestr[0:4])
                        mth = int(stdatetimestr[4:6])
                        ddy = int(stdatetimestr[6:8])
                        hr = int(stdatetimestr[8:10])
                        mn = int(stdatetimestr[10:12])
                        ControlRec.StartDateTime = dt.datetime(yr, mth, ddy, hr, mn, 00)
                        ControlRec.NonGSNumber = 'E' + jobNumber
                        if ln[0:13] == '*CONTROL;DIS;':
                            ControlRec.InstrumentModel = alinestr[2].strip()
                            ControlRec.Class = alinestr[3].strip()
                            ControlRec.Derivation = alinestr[4].strip()
                        if ln[0:13] == '*CONTROL;ANG;':
                            ControlRec.InstrumentModel = alinestr[2].strip()
                            ControlRec.Class = alinestr[3].strip()
                            ControlRec.Derivation = alinestr[4].strip()
                        if ln[0:14] == '*CONTROL;OHDF;':
                            ControlRec.SurveyTechnique = alinestr[3].strip()
                            ControlRec.Class = alinestr[4].strip()
                            ControlRec.LevelDistance = alinestr[5].strip()
                            ControlRec.Derivation = alinestr[6].strip()
                        if ln[0:13] == '*CONTROL;GPS;':
                            durationstr = alinestr[3]
                            hr = int(durationstr[0:2])
                            mn = int(durationstr[2:4])
                            sec = int(durationstr[4:6])
                            ControlRec.Duration = (dt.datetime(1900, 1, 1, 0, 0, 0) +
                                                   dt.timedelta(hours=hr,
                                                                minutes=mn,
                                                                seconds=sec))
                            ControlRec.TimeStatus = alinestr[4].strip()
                            ControlRec.EphemerisType = alinestr[5].strip()
                            ControlRec.AtReceiver = alinestr[6].strip()
                            ControlRec.ToReceiver = alinestr[7].strip()
                            ControlRec.FrequencyMode = alinestr[8].strip()
                            ControlRec.SurveyTechnique = alinestr[9].strip()
                            ControlRec.Solution = alinestr[10].strip()
                            ControlRec.EpochInterval = int(alinestr[11])
                            ControlRec.Class = alinestr[12].strip()
        
                    if ln[0:4] == ' GRP':
                        CurrentMsr.type = 'G'
                    if ln[0:5] == ' DXYZ':
                        CurrentMsr.first = ln[10:23].strip()
                        CurrentMsr.second = ln[23:35].strip()
                        CurrentMsr.dx = ln[36:50].strip()
                        CurrentMsr.dy = ln[50:64].strip()
                        CurrentMsr.dz = ln[64:78].strip()
                    if ln[0:6] == ' COV  ' or ln[0:6] == ' CORR ':
                        m_type = ln[0:13]
                        lineCount = 0
                        if ln[26:36].strip()!='': CurrentMsr.vscale = ln[26:36].strip()
                        diagscale=1
                        if ln[48:58].strip()!='': diagscale = float(ln[48:58].strip())
                        if diagscale == 0: diagscale=1
                    if ln[0:5] == ' ELEM' and CurrentMsr.type == 'G' and lineCount == 1:
                        CurrentMsr.Vs[0, 0] = float(ln[7:30].strip())*diagscale
                        if m_type== ' COV  LG DIAG':
                            CurrentMsr.Vs[1, 1] = ln[30:54].strip()
                            CurrentMsr.Vs[2, 2] = ln[54:78].strip()
                            msrout.write(msr_xml_str(CurrentMsr, ControlRec))
                        else:    
                            CurrentMsr.Vs[0, 1] = ln[30:54].strip()
                            CurrentMsr.Vs[0, 2] = ln[54:78].strip()
                    if ln[0:5] == ' ELEM' and CurrentMsr.type == 'G' and lineCount == 2:
                        CurrentMsr.Vs[1, 1] = float(ln[7:30].strip())*diagscale
                        CurrentMsr.Vs[1, 2] = ln[30:54].strip()
                    if ln[0:5] == ' ELEM' and CurrentMsr.type == 'G' and lineCount == 3:
                        CurrentMsr.Vs[2, 2] = float(ln[7:30].strip())*diagscale
                        if m_type== ' COV  CT UPPR':
                            msrout.write(msr_xml_str(CurrentMsr, ControlRec))
                    if ln[0:5] == ' ELEM' and CurrentMsr.type == 'G' and lineCount == 4:
                        Ds=np.zeros([3,3])
                        Ds[0,0] = ln[7:30].strip()
                        Ds[1,1] = ln[30:54].strip()
                        Ds[2,2] = ln[54:78].strip()
                        CurrentMsr.Vs[1, 0] =CurrentMsr.Vs[0, 1]
                        CurrentMsr.Vs[2, 0] =CurrentMsr.Vs[0, 2]
                        CurrentMsr.Vs[2, 1] =CurrentMsr.Vs[1, 2]
                        CurrentMsr.Vs=np.matmul(np.matmul(Ds,CurrentMsr.Vs),Ds)
                        msrout.write(msr_xml_str(CurrentMsr, ControlRec))
        
                    lineCount = lineCount+1
        
                # Write footers, close
                stnout.write(dml_footer())
                msrout.write(dml_footer())
                stnout.close()
                msrout.close()

                ####################################
                ### Run the DynAdjust Adjustment ###
                ####################################
                print('  Adjusting: ' + adjustment_name)
                subprocess.run("import -n " + adjustment_name + " " + adjustment_name + ".msr.xml "+ adjustment_name + ".stn.xml --flag-unused-stations --remove-ignored-msr -r " + datum)
                subprocess.run("reftran " + adjustment_name + " -r " + datum)
                subprocess.run("geoid " + adjustment_name + " -g \""+geoid+"\" --convert-stn-hts")
                subprocess.run("segment " + adjustment_name + " --min-inner-stns 500 --max-block-stns 500")
                subprocess.run("adjust " + adjustment_name + " --staged --create-stage-files --output-adj-msr --output-pos-uncertainty --output-adj-gnss-units 0 --max-iterations 20 --free-stn-sd 5 --iteration-threshold 0.0005 --stn-coord-types ENzPLHhXYZ")

                coords, bases = DynaGeol(adjustment_name + '.phased-stage.adj')
                if len(bases) > 0: Create_DynaBAS(adjustment_name + '.phased-stage.adj',bases)
                
    files = os.listdir(os.getcwd())
    for f in files:
        if f.endswith('.aml'):os.remove(f)
        if f.endswith('.asl'):os.remove(f)
        if f.endswith('.bms'):os.remove(f)
        if f.endswith('.bst'):os.remove(f)
        if f.endswith('.dbid'):os.remove(f)
        if f.endswith('.dnaproj'):os.remove(f)
        if f.endswith('.imp'):os.remove(f)
        if f.endswith('.map'):os.remove(f)
        if f.endswith('.mtx'):os.remove(f)
        if f.endswith('.seg'):os.remove(f)
        if f.endswith('.xyz'):os.remove(f)
        if f.endswith('.db'):os.remove(f)
        if f.endswith('.adj') or f.endswith('.apu'):
            shutil .move(f, f.replace('.phased-stage',''))
    
    print('\n============== Complete ==============')
