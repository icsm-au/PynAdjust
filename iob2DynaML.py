from math import sqrt
import math
import numpy
from numpy import matmul
from math import sin, cos, radians
import sys,datetime, os

class DnaStation:
    def __init__(self):
        self.name = ''
        self.Constraint=''
        self.Type=''
        self.XAxis=''
        self.YAxis=''
        self.Height=''
        self.Description=''
        self.aAxis=0
        self.bAxis=0
        self.ErrAz=0
class AdditionalInfoStn:
    def __init__(self):
        self.HorizCoordMethod=''
        self.RelativeHorizAccuracy=''
        self.NonGSNumber=''
        self.SelectPoint='true'
        self.SelectRL='true'
class AdditionalInfoMsr:
    def __init__(self):
        #Additional information is included as a comment in the DynaML file. This can be used for database import
        self.StartDateTime=datetime.datetime(1994, 1, 1, 00, 00,00)
        self.Duration=datetime.datetime(1966, 1, 1, 00, 00,00)
        self.TimeStatus=''
        self.EphemerisType=''
        self.AtReceiver=''
        self.ToReceiver=''
        self.FrequencyMode=''
        self.SurveyTechnique='SLEV'
        self.Solution=''
        self.EpochInterval=''
        self.Class='LC'
        self.LevelDistance='0.01'
        self.InstrumentModel=''
        self.Derivation='MEAS'
        self.NonGSNumber=''
class DnaMeasurement:
    def __init__(self):
        self.type = ''
        self.vscale='1'
        self.pscale='1'
        self.lscale='1'
        self.hscale='1'
        self.first=''
        self.second=''
        self.third=''
        self.stddev=''
        self.total=''
        self.instheight=0
        self.targheight=0
        self.targets=[]
        self.values=[]
        self.targetstddevs=[]
        self.dx=''
        self.dy=''
        self.dz=''
        self.Vs=numpy.zeros([3,3])
    def add_targets(self,target):
        self.targets.append(target)
    def add_values(self,Value):
        self.values.append(Value)
    def add_targetstddevs(self,targetstddev):
        self.targetstddevs.append(targetstddev)
        
class DeviceHeight:
    def __init__(self):
        #Device Height might be the height of instrument at a point or Height of target
        self.StnName=[]
        self.RefHeight=[]
def add_DeviceHeight(self,Stn,Hgt):
    self.StnName.append(Stn)
    self.RefHeight.append(Hgt)
def hms2hp(HMS_Ang):
    #Input: HH MM SS.ssss used by Geolab
    #Output: HH.MMSSSsssss used by DynAdjust
    sign=1
    if HMS_Ang.find('S')<>-1 or HMS_Ang.find('-')<>-1:
        sign=-1
    while HMS_Ang.find('  ')<>-1:
        HMS_Ang=HMS_Ang.replace('  ',' ')
    HMS_Ang=HMS_Ang.replace('S','')
    HMS_Ang=HMS_Ang.replace('E','')
    HMS_Ang=HMS_Ang.replace('.',' ')
    aAng=HMS_Ang.split()
    aAng[0]=str(sign*abs(int(aAng[0])))
    aAng[1]="%02d" % float(aAng[1])
    aAng[2]="%02d" % float(aAng[2])
    return aAng[0] + '.' + aAng[1] + aAng[2]+ aAng[3]

def hp2dec(hp):
    #Input: HH.MMSSsss
    #Output: dd.dddddd
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    dec = degree + (minute / 60) + (second / 360)
    return dec if hp >= 0 else -dec

def FindJobNumber(strg):
    #search a string for 8 consecutive numbers, this is probably the Job Number
    JN=''
    i=0
    while i+7<>len(strg):
        if unicode(strg[i:i+8]).isnumeric()==True:
            JN=strg[i:i+8]
        i=i+1
    return JN
def Stn_xml_str(Stn,stnAdditionalRec):
    #Output: String for printing to xml that is one complete station
    xml_str='<DnaStation>\n'
    xml_str=xml_str+'<Name>' + Stn.Name + '</Name>\n'
    xml_str=xml_str+'<Constraints>' + Stn.Constraint + '</Constraints>\n'
    xml_str=xml_str+'<Type>' + Stn.Type + '</Type>\n'
    xml_str=xml_str+'<StationCoord>\n'
    xml_str=xml_str+'<Name>' + Stn.Name + '</Name>\n'
    xml_str=xml_str+'<XAxis>' + str(Stn.XAxis) + '</XAxis>\n'
    xml_str=xml_str+'<YAxis>' + str(Stn.YAxis) + '</YAxis>\n'
    xml_str=xml_str+'<Height>' + Stn.Height + '</Height>\n'
    xml_str=xml_str+'</StationCoord>\n'
    xml_str=xml_str+'<Description>'+ Stn.Description +'</Description>\n'
    xml_str=xml_str+'<!--AdditionalInfoStn>\n'
    xml_str=xml_str+'<HorizCoordMethod>' + stnAdditionalRec.HorizCoordMethod + '</HorizCoordMethod>\n'
    xml_str=xml_str+'<RelativeHorizAccuracy>' + stnAdditionalRec.RelativeHorizAccuracy + '</RelativeHorizAccuracy>\n'
    xml_str=xml_str+'<NonGSNumber>' + stnAdditionalRec.NonGSNumber + '</NonGSNumber>\n'
    xml_str=xml_str+'<SelectPoint>' + stnAdditionalRec.SelectPoint + '</SelectPoint>\n'
    xml_str=xml_str+'<SelectRL>' + stnAdditionalRec.SelectRL + '</SelectRL>\n'
    xml_str=xml_str+'</AdditionalInfoStn-->\n'    
    xml_str=xml_str+'</DnaStation>\n'
    return xml_str

def Msr_xml_str(Msr, ControlRec):
    #Output: xml string for printing to file. Caters for type G, D, S, B, D, L, H
    xml_str='<DnaMeasurement>\n'
    xml_str=xml_str+'<Type>' + Msr.type + '</Type>\n'
    xml_str=xml_str+'<Ignore/>\n'
    if Msr.type == 'G':
        xml_str=xml_str+'<ReferenceFrame>' + GNSSdate2Ref(ControlRec.StartDateTime) + '</ReferenceFrame>\n'
        xml_str=xml_str+'<Epoch>' + ControlRec.StartDateTime.strftime('%d.%m.%Y') + '</Epoch>\n'
        xml_str=xml_str+'<Vscale>' + Msr.vscale + '</Vscale>\n'
        xml_str=xml_str+'<Pscale>' + Msr.pscale + '</Pscale>\n'
        xml_str=xml_str+'<Lscale>' + Msr.lscale + '</Lscale>\n'
        xml_str=xml_str+'<Hscale>' + Msr.hscale + '</Hscale>\n'
    xml_str=xml_str+'<First>' + Msr.first + '</First>\n'
    if Msr.second <> '':
        xml_str=xml_str+'<Second>' + Msr.second + '</Second>\n'
    if Msr.type <> 'G' and Msr.type <> 'D':
        xml_str=xml_str+'<Value>' + Msr.values[0] + '</Value>\n'
        xml_str=xml_str+'<StdDev>' + Msr.stddev + '</StdDev>\n'
    if Msr.type == 'G':
        xml_str=xml_str+'<GPSBaseline>\n'
        xml_str=xml_str+'<X>' + Msr.dx + '</X>\n'
        xml_str=xml_str+'<Y>' + Msr.dy + '</Y>\n'
        xml_str=xml_str+'<Z>' + Msr.dz + '</Z>\n'
        xml_str=xml_str+'<SigmaXX>' + str(Msr.Vs[0,0]) + '</SigmaXX>\n'
        xml_str=xml_str+'<SigmaXY>' + str(Msr.Vs[0,1]) + '</SigmaXY>\n'
        xml_str=xml_str+'<SigmaXZ>' + str(Msr.Vs[0,2]) + '</SigmaXZ>\n'
        xml_str=xml_str+'<SigmaYY>' + str(Msr.Vs[1,0]) + '</SigmaYY>\n'
        xml_str=xml_str+'<SigmaYZ>' + str(Msr.Vs[1,1]) + '</SigmaYZ>\n'
        xml_str=xml_str+'<SigmaZZ>' + str(Msr.Vs[2,0]) + '</SigmaZZ>\n'
        xml_str=xml_str+'</GPSBaseline>\n'
        xml_str=xml_str+'<!--AdditionalInfoMsrG>\n'
        xml_str=xml_str+'<StartDateTime>'+ControlRec.StartDateTime.strftime('%Y-%m-%dT%H:%M:%S')+'</StartDateTime>\n'
        xml_str=xml_str+'<Duration>P'+str(ControlRec.Duration.year-1900)+'Y'+str(ControlRec.Duration.month-1)+'M'+str(ControlRec.Duration.day-1)+'DT'+ControlRec.Duration.strftime('%HH%MM%SS')+'</Duration>\n'
        xml_str=xml_str+'<TimeStatus>'+ControlRec.TimeStatus+'</TimeStatus>\n'
        xml_str=xml_str+'<EphemerisType>'+ControlRec.EphemerisType+'</EphemerisType>\n'
        xml_str=xml_str+'<AtReceiver>'+ControlRec.AtReceiver+'</AtReceiver>\n'
        xml_str=xml_str+'<ToReceiver>'+ControlRec.ToReceiver+'</ToReceiver>\n'
        xml_str=xml_str+'<FrequencyMode>'+ControlRec.FrequencyMode+'</FrequencyMode>\n'
        xml_str=xml_str+'<SurveyTechnique>'+ControlRec.SurveyTechnique+'</SurveyTechnique>\n'
        xml_str=xml_str+'<Solution>'+ControlRec.Solution+'</Solution>\n'
        xml_str=xml_str+'<EpochInterval>'+str(ControlRec.EpochInterval)+'</EpochInterval>\n'
        xml_str=xml_str+'<Class>'+ControlRec.Class+'</Class>\n'
        xml_str=xml_str+'<NonGSNumber>'+ControlRec.NonGSNumber+'</NonGSNumber>\n'
        xml_str=xml_str+'</AdditionalInfoMsrG-->\n'
    
    if Msr.type == 'L':
        xml_str=xml_str+'<!--AdditionalInfoMsrL>\n'
        xml_str=xml_str+'<SurveyTechnique>'+ControlRec.SurveyTechnique+'</SurveyTechnique>\n'
        xml_str=xml_str+'<LevelDistance>'+ ControlRec.LevelDistance +'</LevelDistance>\n'
        xml_str=xml_str+'<ObsDate>'+ControlRec.StartDateTime.strftime('%Y-%m-%d')+'</ObsDate>\n'
        xml_str=xml_str+'<Derivation>'+ControlRec.Derivation+'</Derivation>\n'
        xml_str=xml_str+'<Class>'+ControlRec.Class+'</Class>\n'
        xml_str=xml_str+'<NonGSNumber>'+ControlRec.NonGSNumber+'</NonGSNumber>\n'
        xml_str=xml_str+'</AdditionalInfoMsrL-->\n'
    
    if Msr.type == 'S':
        xml_str=xml_str+'<InstHeight>' + str(Msr.instheight) + '</InstHeight>\n'
        xml_str=xml_str+'<TargHeight>' + str(Msr.targheight) + '</TargHeight>\n'
        xml_str=xml_str+'<!--AdditionalInfoMsrS>\n'
        xml_str=xml_str+'<InstrumentModel>'+ControlRec.InstrumentModel+'</InstrumentModel>\n'
        xml_str=xml_str+'<ObsDate>'+ControlRec.StartDateTime.strftime('%Y-%m-%d')+'</ObsDate>\n'
        xml_str=xml_str+'<Derivation>'+ControlRec.Derivation+'</Derivation>\n'
        xml_str=xml_str+'<Class>'+ControlRec.Class+'</Class>\n'
        xml_str=xml_str+'<NonGSNumber>'+ControlRec.NonGSNumber+'</NonGSNumber>\n'
        xml_str=xml_str+'</AdditionalInfoMsrS-->\n'

    if Msr.type == 'D':
        xml_str=xml_str+'<Value>' + Msr.values[0] + '</Value>\n'
        xml_str=xml_str+'<StdDev>' + Msr.targetstddevs[0] + '</StdDev>\n'
        xml_str=xml_str+'<Total>' + str(Msr.total-1) + '</Total>\n'
        ObsNumber=1
        while ObsNumber<Msr.total:
            xml_str=xml_str+'<Directions>\n'
            xml_str=xml_str+'<Ignore/>\n'
            xml_str=xml_str+'<Target>' + Msr.targets[ObsNumber] + '</Target>\n'
            xml_str=xml_str+'<Value>' + Msr.values[ObsNumber] + '</Value>\n'
            xml_str=xml_str+'<StdDev>' + Msr.targetstddevs[ObsNumber] + '</StdDev>\n'
            xml_str=xml_str+'</Directions>\n'
            ObsNumber=ObsNumber+1
        xml_str=xml_str+'<!--AdditionalInfoMsrD>\n'
        xml_str=xml_str+'<InstrumentModel>'+ControlRec.InstrumentModel+'</InstrumentModel>\n'
        xml_str=xml_str+'<ObsDate>'+ControlRec.StartDateTime.strftime('%Y-%m-%d')+'</ObsDate>\n'
        xml_str=xml_str+'<Derivation>'+ControlRec.Derivation+'</Derivation>\n'
        xml_str=xml_str+'<Class>'+ControlRec.Class+'</Class>\n'
        xml_str=xml_str+'<NonGSNumber>'+ControlRec.NonGSNumber+'</NonGSNumber>\n'
        xml_str=xml_str+'</AdditionalInfoMsrD-->\n'
    xml_str=xml_str+'<Source></Source>\n'
    xml_str=xml_str+'</DnaMeasurement>\n'
    return xml_str

c_vac = 299792.458
k_0 = 0.9996

# Ellipsoid Constants
class Ellipsoid(object):
    def __init__(self, semimaj, inversef):
        self.semimaj = semimaj
        self.inversef = inversef
        self.f = 1 / self.inversef
        self.semimin = float(self.semimaj * (1 - self.f))
        self.ecc1sq = float(self.f * (2 - self.f))
        self.ecc2sq = float(self.ecc1sq / (1 - self.ecc1sq))
        self.ecc1 = sqrt(self.ecc1sq)
        self.n = float(self.f / (2 - self.f))
        self.n2 = self.n ** 2

# Geodetic Reference System 1980
grs80 = Ellipsoid(6378137, 298.25722210088)

def llh2xyz(lat, lng, ellht, ellipsoid=grs80):
    # Add input for ellipsoid (default: grs80)
    # Convert lat & long to radians
    lat = radians(hp2dec(float(lat)))
    lng = radians(hp2dec(float(lng)))
    ellht=float(ellht)
    # Calculate Ellipsoid Radius of Curvature in the Prime Vertical - nu
    if lat == 0:
        nu = grs80.semimaj
    else:
        nu = ellipsoid.semimaj/(sqrt(1 - ellipsoid.ecc1sq * (sin(lat)**2)))
    # Calculate x, y, z
    x = (nu + ellht) * cos(lat) * cos(lng)
    y = (nu + ellht) * cos(lat) * sin(lng)
    z = ((ellipsoid.semimin**2 / ellipsoid.semimaj**2) * nu + ellht) * sin(lat)
    return x, y, z

def ErrEllip2Ycluster(Stn):
    #Input: Supply a station with coordinates and error ellipse for coordinate uncertainty
    #Output: xml string for  point cluster (Y-type observation)
    x, y, z = llh2xyz(Stn.XAxis, Stn.YAxis, Stn.Height)
    
    a=Stn.aAxis/2.44774683068
    b=Stn.bAxis/2.44774683068
    Az=90-Stn.ErrAz
    
    rAz=math.radians(Az)
    rlat=math.radians(float(Stn.XAxis))
    rlng=math.radians(float(Stn.YAxis))
    
    rl=numpy.zeros([3,3])
    rl[0,0]=-sin(rlng)
    rl[0,1]=-sin(rlat)*cos(rlng)
    rl[0,2]=cos(rlat)*cos(rlng)
    rl[1,0]=cos(rlng)
    rl[1,1]=-sin(rlat)*sin(rlng)
    rl[1,2]=cos(rlat)*sin(rlng)
    rl[2,1]=cos(rlat)
    rl[2,2]=sin(rlat)

    iA=numpy.zeros([3,3])
    iA[0,0]=(cos(rAz)*cos(rAz)*a*a)+(b*b*sin(rAz)*sin(rAz))
    iA[0,1]=(a*a-b*b)*cos(rAz)*sin(rAz)
    iA[1,0]=iA[0,1]
    iA[1,1]=(a*a*sin(rAz)*sin(rAz))+(b*b*cos(rAz)*cos(rAz))
    iA[2,2]=0.000001
    
    Wt=matmul(matmul(rl,iA),rl.transpose())
    
    xml_str='<DnaMeasurement>\n'
    xml_str=xml_str+'<Type>Y</Type>\n'
    xml_str=xml_str+'<Ignore/>\n'
    xml_str=xml_str+'<ReferenceFrame>GDA2020</ReferenceFrame>\n'
    xml_str=xml_str+'<Epoch>01.01.2020</Epoch>\n'
    xml_str=xml_str+'<Vscale>1.000</Vscale>\n'
    xml_str=xml_str+'<Pscale>1.000</Pscale>\n'
    xml_str=xml_str+'<Lscale>1.000</Lscale>\n'
    xml_str=xml_str+'<Hscale>1.000</Hscale>\n'
    xml_str=xml_str+'<Coords>XYZ</Coords>\n'
    xml_str=xml_str+'<Total>1</Total>\n'
    xml_str=xml_str+'<First>' + Stn.Name + '</First>\n'
    xml_str=xml_str+'<Clusterpoint>\n'
    xml_str=xml_str+'<X>'+str(x)+'</X>\n'
    xml_str=xml_str+'<Y>'+str(y)+'</Y>\n'
    xml_str=xml_str+'<Z>'+str(z)+'</Z>\n'
    xml_str=xml_str+'<SigmaXX>'+str(Wt[0,0])+'</SigmaXX>\n'
    xml_str=xml_str+'<SigmaXY>'+str(Wt[0,1])+'</SigmaXY>\n'
    xml_str=xml_str+'<SigmaXZ>'+str(Wt[0,2])+'</SigmaXZ>\n'
    xml_str=xml_str+'<SigmaYY>'+str(Wt[1,1])+'</SigmaYY>\n'
    xml_str=xml_str+'<SigmaYZ>'+str(Wt[1,2])+'</SigmaYZ>\n'
    xml_str=xml_str+'<SigmaZZ>'+str(Wt[2,2])+'</SigmaZZ>\n'
    xml_str=xml_str+'</Clusterpoint>\n'
    xml_str=xml_str+'</DnaMeasurement>\n'
    
    return xml_str

def stn_header():
    xml_str='<?xml version="1.0" encoding="utf-8"?>\n'
    xml_str=xml_str+'<DnaXmlFormat type="Station File" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="DynaML.xsd">\n'
    return xml_str
def msr_header():
    xml_str='<?xml version="1.0" encoding="utf-8"?>\n'
    xml_str=xml_str+'<DnaXmlFormat type="Measurement File" referenceframe="GDA94" epoch="01.01.1994" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="DynaML.xsd">\n'
    return xml_str
def dML_footer():
    xml_str='</DnaXmlFormat>\n'
    return xml_str
def GNSSdate2Ref(obsDate):
    #Use the date of GNSS baseline obsedrvation to determine the reference frame used by broadcast ephemeris
    if obsDate >= datetime.datetime(1900, 1, 1) and obsDate < datetime.datetime(1994, 1, 2):    
        Ref = 'ITRF1991'
    if obsDate >= datetime.datetime(1994, 1, 2) and obsDate < datetime.datetime(1995, 1, 1):    
        Ref = 'ITRF1992'
    if obsDate >= datetime.datetime(1995, 1, 1) and obsDate < datetime.datetime(1996, 6, 30):    
        Ref = 'ITRF1993'
    if obsDate >= datetime.datetime(1996, 6, 30) and obsDate < datetime.datetime(1998, 3, 1):    
        Ref = 'ITRF1994'
    if obsDate >= datetime.datetime(1998, 3, 1) and obsDate < datetime.datetime(1999, 8, 1):    
        Ref = 'ITRF1996'
    if obsDate >= datetime.datetime(1999, 8, 1) and obsDate < datetime.datetime(2001, 12, 2):    
        Ref = 'ITRF1997'
    if obsDate >= datetime.datetime(2001, 12, 2) and obsDate < datetime.datetime(2006, 11, 5):    
        Ref = 'ITRF2000'
    if obsDate >= datetime.datetime(2006, 11, 5) and obsDate < datetime.datetime(2011, 4, 17):    
        Ref = 'ITRF2005'
    if obsDate >= datetime.datetime(2011, 4, 17):    
        Ref = 'ITRF2008'
    return Ref
#####################################################################################
#### Input:Geolab *.iob file                                                    #####
#### Output: DynaML stn and msr file, If the input file contains the word final #####
####         and the inpuut contains error ellipse information for fixed marks  #####
####         an extra stn and msr file will be produced for doing a weighted    #####
####         adjustmanets and 2-D propogation of uncertainty.                   #####
#####################################################################################

filename = '20192277_3.Final_Adjustment.iob'
if len(sys.argv) >1: filename = sys.argv[1]
f = open(filename, 'r')
SumRelativeUncertainties=0
AvgRelativeUncertainties=0
stnCnt=0
#Open,run through and close the file for initial informatiokn on the adjustment
GNSSmarksStr=';'
FloatmarksStr=';'
w_MksWithObs=';'
for linestr in f.readlines():
    if linestr[0:4]==' PLO' or linestr[0:4] == ' PLH':
        if linestr[72:len(linestr)].strip()<>'' and linestr.find('ppm')<>-1:
            stnRec=linestr[72:len(linestr)].strip().replace(';','|').split('|')
            SumRelativeUncertainties=SumRelativeUncertainties+float(stnRec[7].replace('ppm',''))
            stnCnt=stnCnt+1
        if linestr[6:8]=='00' and filename.lower().find('final')<>-1:FloatmarksStr=FloatmarksStr + linestr[10:23].strip() +';'
    if linestr[0:5] == ' OHDF' or linestr[0:5] == ' GAZI' or linestr[0:5] == ' AZIM' or linestr[0:5] == ' DIST' or linestr[0:5] == ' DSET' or linestr[0:4] == ' DIR' or linestr[0:5] == ' DXYZ':
        CurrentMsr=DnaMeasurement()
        CurrentMsr.first=linestr[10:23].strip()
        CurrentMsr.second=linestr[23:35].strip()
        if linestr[0:5] == ' DXYZ':        
            GNSSmarksStr=GNSSmarksStr+';'+CurrentMsr.first+';'+CurrentMsr.second
        if FloatmarksStr.find(';' + CurrentMsr.first+';')<>-1 or FloatmarksStr.find(';' + CurrentMsr.second+';')<>-1:
             w_MksWithObs=w_MksWithObs+CurrentMsr.first+';'+CurrentMsr.second+';'
f.close

#bin the guessing of GDA94 relative uncertainty for the float marks
if stnCnt<>0:
    AvgRelativeUncertainties=SumRelativeUncertainties/stnCnt
    if AvgRelativeUncertainties<3:
        AvgRelativeUncertainties=3
    if AvgRelativeUncertainties>3 and AvgRelativeUncertainties<=7:
        AvgRelativeUncertainties=7.5
    if AvgRelativeUncertainties>7.5 and AvgRelativeUncertainties<=10:
        AvgRelativeUncertainties=10
    if AvgRelativeUncertainties>10 and AvgRelativeUncertainties<=20:
        AvgRelativeUncertainties=20
    if AvgRelativeUncertainties>20 and AvgRelativeUncertainties<=30:
        AvgRelativeUncertainties=30
    if AvgRelativeUncertainties>30 and AvgRelativeUncertainties<=50:
        AvgRelativeUncertainties=50
    if AvgRelativeUncertainties>50:
        AvgRelativeUncertainties=int(AvgRelativeUncertainties/10)*10
f = open(filename, 'r')
stnout = open(filename.replace('.iob', '.stn.xml'), 'w')
msrout = open(filename.replace('.iob', '.msr.xml'), 'w')
stnout.write(stn_header())
msrout.write(msr_header())

if filename.lower().find('final')<>-1:
    w_stnout = open(filename.replace('.iob', '_W.stn.xml'), 'w')
    w_msrout = open(filename.replace('.iob', '_W.msr.xml'), 'w')
    w_stnout.write(stn_header())
    w_msrout.write(msr_header())
# Run through each line of the input file and extract the relevant lines
lineCount=0
InstHts=DeviceHeight()
TgtHts=DeviceHeight()
CurrentMsr=DnaMeasurement()
ControlRec=AdditionalInfoMsr()
for linestr in f.readlines():
    print(linestr)
    if linestr[0:5]==' TITL':
        jobNumber = FindJobNumber(linestr)
        if jobNumber=='':
            jobNumber=FindJobNumber(os.getcwd())
    if linestr[0:4]==' PLO' or linestr[0:4] == ' PLH':
        CurrentStn=DnaStation()
        CurrentStn.Name=linestr[10:23].strip()
        CurrentStn.Constraint=linestr[6:9].replace('1','C')
        CurrentStn.Constraint=CurrentStn.Constraint.replace('0','F')
        CurrentStn.Constraint=CurrentStn.Constraint.strip()
        if linestr[0:4] == ' PLO':
            CurrentStn.Type='LLH'
        if linestr[0:4] == ' PLH':
            CurrentStn.Type='LLh'
        CurrentStn.XAxis=hms2hp(linestr[23:41].strip())
        CurrentStn.YAxis=hms2hp(linestr[41:59].strip())
        CurrentStn.Height=linestr[59:72].strip()
        CurrentStn.Description=linestr[72:len(linestr)].strip()
        stnAdditionalRec=AdditionalInfoStn()
        stnAdditionalRec.NonGSNumber='E'+jobNumber
        if CurrentStn.Description<>'':
            stnRec=CurrentStn.Description.replace(';','|').split('|')
            stnAdditionalRec.HorizCoordMethod=stnRec[8]
            stnAdditionalRec.RelativeHorizAccuracy=stnRec[7]
            if stnRec[10]!='':
                CurrentStn.aAxis=float(stnRec[10])
                CurrentStn.bAxis=float(stnRec[11])
                CurrentStn.ErrAz=float(stnRec[12])
            stnAdditionalRec.SelectPoint='true'
            stnAdditionalRec.SelectRL='false'

        if CurrentStn.Constraint[0:2]=='FF' and AvgRelativeUncertainties !=0:
            stnAdditionalRec.RelativeHorizAccuracy=str(AvgRelativeUncertainties)+'ppm'
            if GNSSmarksStr.find(';'+CurrentStn.Name+';'):
                stnAdditionalRec.HorizCoordMethod='GNSS'
                
        stnout.write(Stn_xml_str(CurrentStn,stnAdditionalRec))
        if filename.lower().find('final')<>-1:
            w_CurrentStn=CurrentStn
            if w_MksWithObs.find(';' + w_CurrentStn.Name+';')<>-1:
                if CurrentStn.Constraint[0:2]<>'FF'and  CurrentStn.aAxis<>0:
                    w_CurrentStn.Constraint='FF' + CurrentStn.Constraint[2:3]
                    w_msrout.write(ErrEllip2Ycluster(w_CurrentStn))
                w_stnout.write(Stn_xml_str(w_CurrentStn,stnAdditionalRec))
                
    if linestr[0:5] == ' HI  ':
        add_DeviceHeight(InstHts,linestr[10:23].strip(),linestr[23:33].strip())

    if linestr[0:5] == ' HT  ':
        add_DeviceHeight(TgtHts,linestr[10:23].strip(),linestr[23:33].strip())

    if linestr[0:5] == ' OHGT':
        CurrentMsr=DnaMeasurement()
        CurrentMsr.type='H'
        CurrentMsr.first=linestr[10:23].strip()
        CurrentMsr.add_values(linestr[36:65].strip())
        CurrentMsr.stddev=linestr[65:76].strip()
        msrout.write(Msr_xml_str(CurrentMsr,ControlRec))
        
    if linestr[0:5] == ' OHDF':
        CurrentMsr=DnaMeasurement()
        CurrentMsr.type='L'
        CurrentMsr.first=linestr[10:23].strip()
        CurrentMsr.second=linestr[23:35].strip()
        CurrentMsr.add_values(linestr[50:65].strip())
        CurrentMsr.stddev=linestr[65:76].strip()
        ControlRec.LevelDistance=linestr[36:50].strip()
        msrout.write(Msr_xml_str(CurrentMsr,ControlRec))
        if FloatmarksStr.find(';' + CurrentMsr.first+';')<>-1 or FloatmarksStr.find(';' + CurrentMsr.second+';')<>-1:
            w_msrout.write(Msr_xml_str(CurrentMsr,ControlRec))

    if linestr[0:5] == ' GAZI':
        CurrentMsr=DnaMeasurement()
        CurrentMsr.type='B'
        CurrentMsr.first=linestr[10:23].strip()
        CurrentMsr.second=linestr[23:35].strip()
        CurrentMsr.add_values(hms2hp(linestr[36:65].strip()))
        CurrentMsr.stddev=linestr[65:76].strip()
        msrout.write(Msr_xml_str(CurrentMsr,ControlRec))
        if FloatmarksStr.find(';' + CurrentMsr.first+';')<>-1 or FloatmarksStr.find(';' + CurrentMsr.second+';')<>-1:
            w_msrout.write(Msr_xml_str(CurrentMsr,ControlRec))

    if linestr[0:5] == ' AZIM':
        CurrentMsr=DnaMeasurement()
        CurrentMsr.type='K'
        CurrentMsr.first=linestr[10:23].strip()
        CurrentMsr.second=linestr[23:35].strip()
        CurrentMsr.add_values(hms2hp(linestr[36:65].strip()))
        CurrentMsr.stddev=linestr[65:76].strip()
        msrout.write(Msr_xml_str(CurrentMsr,ControlRec))
        if FloatmarksStr.find(';' + CurrentMsr.first+';')<>-1 or FloatmarksStr.find(';' + CurrentMsr.second+';')<>-1:
            w_msrout.write(Msr_xml_str(CurrentMsr,ControlRec))

    if linestr[0:5] == ' DIST':
        CurrentMsr=DnaMeasurement()
        CurrentMsr.type='S'
        CurrentMsr.first=linestr[10:23].strip()
        CurrentMsr.second=linestr[23:35].strip()
        CurrentMsr.add_values(linestr[36:65].strip())
        rw=0
        for Stn in InstHts.StnName:
            if Stn==CurrentMsr.first:
                CurrentMsr.instheight=InstHts.RefHeight[rw]
            rw=rw+1
        rw=0
        for Stn in TgtHts.StnName:
            if Stn==CurrentMsr.first:
                CurrentMsr.targheight=TgtHts.RefHeight[rw]
            rw=rw+1
        CurrentMsr.stddev=linestr[65:76].strip()
        msrout.write(Msr_xml_str(CurrentMsr,ControlRec))
        if FloatmarksStr.find(';' + CurrentMsr.first+';')<>-1 or FloatmarksStr.find(';' + CurrentMsr.second+';')<>-1:
            w_msrout.write(Msr_xml_str(CurrentMsr,ControlRec))
            print FloatmarksStr, CurrentMsr.first, CurrentMsr.second
    if linestr[0:5] == ' DSET':
        CurrentMsr=DnaMeasurement()
        CurrentMsr.type='D'
        lineCount=0
    if linestr[0:4] == ' DIR' and lineCount==1 and CurrentMsr.type=='D':
        CurrentMsr.first=linestr[10:23].strip()
        CurrentMsr.second=linestr[23:35].strip()
        CurrentMsr.add_targets(linestr[23:35].strip())
        CurrentMsr.add_values(hms2hp(linestr[36:65].strip()))
        CurrentMsr.add_targetstddevs(linestr[65:76].strip())
        CurrentMsr.total=lineCount
    if linestr[0:4] == ' DIR' and lineCount>1 and CurrentMsr.type=='D':        
        CurrentMsr.add_targets(linestr[23:35].strip())
        CurrentMsr.add_values(hms2hp(linestr[36:65].strip()))
        CurrentMsr.add_targetstddevs(linestr[65:76].strip())
        CurrentMsr.total=lineCount
    if CurrentMsr.type=='D' and linestr[0:4] <> ' DIR' and lineCount>1:
        msrout.write(Msr_xml_str(CurrentMsr,ControlRec))
        if FloatmarksStr.find(';' + CurrentMsr.first+';')<>-1 or FloatmarksStr.find(';' + CurrentMsr.second+';')<>-1:
            w_msrout.write(Msr_xml_str(CurrentMsr,ControlRec))
        CurrentMsr=DnaMeasurement()

# Scrap information from Landgate GESMAR control Records
# eg.*CONTROL;GPS;201810010007;012057;E;B;TRIM;TRIM;D;ST  ;FX;015;N
# eg.*CONTROL;OHDF;20181001;SLEV;LC;2.71;MEAS
# eg.*CONTROL;DIS;20181213;TS 16;C;MEAS
# eg.*CONTROL;ANG;20181213;TS16;C;MEAS
    if linestr[0:9]=='*CONTROL;':        
        ControlRec=AdditionalInfoMsr()
        alinestr=linestr.split(';')
        stdatetimestr=alinestr[2] + '0000'
        yr=int(stdatetimestr[0:4])
        mth=int(stdatetimestr[4:6])
        ddy=int(stdatetimestr[6:8])
        hr=int(stdatetimestr[8:10])
        mn=int(stdatetimestr[10:12])
        ControlRec.StartDateTime=datetime.datetime(yr, mth, ddy, hr, mn,00)
        ControlRec.NonGSNumber='E' + jobNumber
        if linestr[0:13]=='*CONTROL;DIS;':
            ControlRec.InstrumentModel=alinestr[2].strip()
            ControlRec.Class=alinestr[3].strip()
            ControlRec.Derivation=alinestr[4].strip()
        if linestr[0:13]=='*CONTROL;ANG;':
            ControlRec.InstrumentModel=alinestr[2].strip()
            ControlRec.Class=alinestr[3].strip()
            ControlRec.Derivation=alinestr[4].strip()
        if linestr[0:14]=='*CONTROL;OHDF;':
            ControlRec.SurveyTechnique=alinestr[3].strip()
            ControlRec.Class=alinestr[4].strip()
            ControlRec.LevelDistance=alinestr[5].strip()
            ControlRec.Derivation=alinestr[6].strip()
        if linestr[0:13]=='*CONTROL;GPS;':
            durationstr=alinestr[3]
            hr=int(durationstr[0:2])
            mn=int(durationstr[2:4])
            sec=int(durationstr[4:6])
            ControlRec.Duration=datetime.datetime(1900,1,1,0,0,0) + datetime.timedelta(hours=hr, minutes= mn, seconds=sec)
            ControlRec.TimeStatus=alinestr[4].strip()
            ControlRec.EphemerisType=alinestr[5].strip()
            ControlRec.AtReceiver=alinestr[6].strip()
            ControlRec.ToReceiver=alinestr[7].strip()
            ControlRec.FrequencyMode=alinestr[8].strip()
            ControlRec.SurveyTechnique=alinestr[9].strip()
            ControlRec.Solution=alinestr[10].strip()
            ControlRec.EpochInterval=int(alinestr[11])
            ControlRec.Class=alinestr[12].strip()
        
    if linestr[0:4] == ' GRP':
        CurrentMsr=DnaMeasurement()
        CurrentMsr.type='G'
    if linestr[0:5] == ' DXYZ':        
        CurrentMsr.first=linestr[10:23].strip()
        CurrentMsr.second=linestr[23:35].strip()
        CurrentMsr.dx=linestr[36:50].strip()
        CurrentMsr.dy=linestr[50:64].strip()
        CurrentMsr.dz=linestr[64:78].strip()
    if linestr[0:13] == ' COV  CT UPPR':        
        lineCount=0
        CurrentMsr.vscale=linestr[26:36].strip()
    if linestr[0:5] == ' ELEM' and CurrentMsr.type=='G' and lineCount==1:        
        CurrentMsr.Vs[0,0]=linestr[7:30].strip()
        CurrentMsr.Vs[0,1]=linestr[30:54].strip()
        CurrentMsr.Vs[0,2]=linestr[54:78].strip()
    if linestr[0:5] == ' ELEM' and CurrentMsr.type=='G' and lineCount==2:        
        CurrentMsr.Vs[1,0]=linestr[7:30].strip()
        CurrentMsr.Vs[1,1]=linestr[30:54].strip()
    if linestr[0:5] == ' ELEM' and CurrentMsr.type=='G' and lineCount==3:        
        CurrentMsr.Vs[2,0]=linestr[7:30].strip()
        msrout.write(Msr_xml_str(CurrentMsr,ControlRec))
        if FloatmarksStr.find(';' + CurrentMsr.first+';')<>-1 or FloatmarksStr.find(';' + CurrentMsr.second+';')<>-1:
            w_msrout.write(Msr_xml_str(CurrentMsr,ControlRec))

    lineCount=lineCount+1
# Close the files
f.close
stnout.write(dML_footer())
msrout.write(dML_footer())
stnout.close()
msrout.close()
if filename.lower().find('final')<>-1:
    w_stnout.write(dML_footer())
    w_msrout.write(dML_footer())
    w_stnout.close()
    w_msrout.close()
print('Done :)')
