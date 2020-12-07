###############################################################################
#                           DynaStdDevStyles.py
#
# ----------------------------------------------------------------------
#   Usage:  cmd:\> python ApplyStdDevStyles.py <*msr.xml>
# ----------------------------------------------------------------------
#
# This code allows the applying of centring errors for instrument setups  
# and allows recalculation of the observation Std Dev based on ppm + Const
# also allows to test the original Std Dev against the newly calculated
# standard deviation to either:
#     - only Increase Original Std Devs that are less than calculated
#     - only Reduce Original Std Devs that are greater than calculated
#     - modify the Std Dev to calculated regardless of original
#
# Code applies only to G Type observations
#
# Style Guides have been added to the msr.xml file as comments.  
# The 'StdDevObsStyle' will overwrite the original vcv matrix
# The 'StdDevSetupStyle' will be added to the vcv matrix.
# The 'IncreaseOrReduce' field is optional.
#
#	<!--StdDevSetupStyle>
#		<StyleName>StdTribrach</StyleName>
#		<CentringStdDev>0.002</CentringStdDev>
#		<VtStdDev>0.003</VtStdDev>
#	</StdDevSetupStyle-->
#	<!--StdDevSetupStyle>
#		<StyleName>Pillar</StyleName>
#		<CentringStdDev>0.0005</CentringStdDev>
#		<VtStdDev>0.0002</VtStdDev>
#	</StdDevSetupStyle-->
#	<!--StdDevObsStyle>
#		<StyleName>NSW</StyleName>
#		<HzPPM>0.7</HzPPM>
#		<HzConstant>0.005</HzConstant>
#		<VtPPM>2.0</VtPPM>
#		<VtConstant>0.015</VtConstant>
#		<IncreaseOrReduce>Increase</IncreaseOrReduce>
#	</StdDevObsStyle-->
#
# These Std Dev Styles are applied to observations by including reference to 
# them in the 'DnaMeasurement'
#
#		<First>BYUP</First>
#		<!--FirstStdDevSetupStyle>StdTribrach</FirstStdDevSetupStyle-->
#		<Second>BYUP_T</Second>
#		<!--SecondStdDevSetupStyle>Pillar</SecondStdDevSetupStyle-->
#		<!--StdDevObsStyle>NSW</StdDevObsStyle-->
#		<GPSBaseline>
#
from geodepy.convert import hp2dec
from geodepy.statistics import vcv_local2cart, vcv_cart2local
from geodepy.geodesy import vincinv
import xml.dom.minidom
import xmltodict
import sys, os
import numpy as np

#'''''''''''''''''''''''''''''''''''''''''''''''
#'''' Assign file Names and Check for errors '''
#'''''''''''''''''''''''''''''''''''''''''''''''
msr_xml='BYUP1.msr.xml'
if len(sys.argv)>1:msr_xml=sys.argv[1]
sd_xml=msr_xml.replace('msr.','sd_msr.')
stn_xml=msr_xml.replace('msr.','stn.')

if msr_xml.endswith('msr.xml')==False: 
    print('Error - Input file must be DynaML msr.xml file')
if os.path.exists(stn_xml)==False:
    print('Error - Input file must have a corresponding DynaML stn.xml file')

# Put the stations and coordintates into a dictionary
with open(stn_xml, 'r') as x:
    xml_s = xmltodict.parse(x.read())
xml_s = xml_s['DnaXmlFormat']['DnaStation']
if type(xml_s)!=list: s=[xml_s]
stn={}
for s in xml_s:
    stn[s['Name']]={
            'P':float(s['StationCoord']['XAxis']),
            'L':float(s['StationCoord']['YAxis'])}


# Put the StdDev Styles into into dictionaries
with open(msr_xml, 'r') as x:
    xml_m = (x.read()
        .replace('<!--StdDev','<StdDev')
        .replace('Style-->','Style>')
        .replace('<!--FirstStdDevSetupStyle','<FirstStdDevSetupStyle')
        .replace('<!--SecondStdDevSetupStyle','<SecondStdDevSetupStyle'))
xml_m = xmltodict.parse(xml_m)

if 'StdDevSetupStyle' in xml_m['DnaXmlFormat']:
    x = xml_m['DnaXmlFormat']['StdDevSetupStyle']
    if type(x)!=list: x=[x]
    sd_styl={}
    for s in x:
        sd_styl[s['StyleName']]={
                'CentringStdDev':float(s['CentringStdDev']),
                'VtStdDev':float(s['VtStdDev'])}

if 'StdDevObsStyle' in xml_m['DnaXmlFormat']:    
    x = xml_m['DnaXmlFormat']['StdDevObsStyle']
    if type(x)!=list: x=[x]
    obs_styl={}
    for s in x:
        obs_styl[s['StyleName']]={
                'HzPPM':float(s['HzPPM']),
                'HzConstant':float(s['HzConstant']),
                'VtPPM':float(s['VtPPM']),
                'VtConstant':float(s['VtConstant'])}

# Put the msr.xml file into a list ... Including comments
xml_m = xml.dom.minidom.parse(msr_xml)
xml_m = xml_m.toprettyxml()
xml_m = [s for s in xml_m.split('\n') if s.strip()!='']


i=0
# Run through the msr file and reprint Std Dev Styled Baselines
with open(sd_xml, 'w') as f_out:
    while i != len(xml_m):
        msr_out = xml_m[i]
        # Stripout the observations or print the line.
        if xml_m[i].find('<DnaMeasurement>')!=-1:
            msr_out=''
            while xml_m[i].find('</DnaMeasurement>')==-1:
                msr_out = msr_out + xml_m[i] +'\n'
                i+=1
            msr_out = msr_out + xml_m[i]
            
            # Test if the measurement needs changing
            m = xmltodict.parse(msr_out
                    .replace('<!--','<')
                    .replace('-->','>'))['DnaMeasurement']
            if (m['Type']=='G' 
              and ('FirstStdDevSetupStyle' in m 
                   or 'SecondStdDevSetupStyle' in m
                   or 'StdDevObsStyle' in m)):
                #find the coordinates of the at and to station
                at_lat   = hp2dec(stn[m['First']]['P'])
                at_lng   = hp2dec(stn[m['First']]['L'])
                to_lat   = hp2dec(stn[m['Second']]['P'])
                to_lng   = hp2dec(stn[m['Second']]['L'])
                
                #Initialise the required matricies
                at_sd    = np.zeros([3, 3])
                to_sd    = np.zeros([3, 3])
                orig_vcv = np.zeros([3, 3])
                cal_vcv  = np.zeros([3, 3])
                
                orig_vcv[0,0] = float(m['GPSBaseline']['SigmaXX'])
                orig_vcv[0,1] = float(m['GPSBaseline']['SigmaXY'])
                orig_vcv[0,2] = float(m['GPSBaseline']['SigmaXZ'])
                orig_vcv[1,1] = float(m['GPSBaseline']['SigmaYY'])
                orig_vcv[1,2] = float(m['GPSBaseline']['SigmaYZ'])
                orig_vcv[2,2] = float(m['GPSBaseline']['SigmaZZ'])
                orig_enu  = vcv_cart2local(orig_vcv, at_lat, at_lng)
                
                # find and apply the centring error styling of the at and to station
                # transform these from enu to xyz
                if 'FirstStdDevSetupStyle' in m:
                    c = sd_styl[m['FirstStdDevSetupStyle']]['CentringStdDev']
                    v = sd_styl[m['FirstStdDevSetupStyle']]['VtStdDev']
                    at_sd[0,0] = (c)**2 / float(m['Lscale'])
                    at_sd[1,1] = (c)**2 / float(m['Pscale'])
                    at_sd[2,2] = (v)**2 / float(m['Hscale'])
                    at_sd = at_sd / float(m['Vscale'])

                if 'SecondStdDevSetupStyle' in m:
                    c = sd_styl[m['SecondStdDevSetupStyle']]['CentringStdDev']
                    v = sd_styl[m['SecondStdDevSetupStyle']]['VtStdDev']
                    to_sd[0,0] = (c)**2 / float(m['Lscale'])
                    to_sd[1,1] = (c)**2 / float(m['Pscale'])
                    to_sd[2,2] = (v)**2 / float(m['Hscale'])
                    to_sd = to_sd / float(m['Vscale'])

                # Sum the 3 matricies and test against orig               
                cal_enu = orig_enu + at_sd + to_sd
                
                # find and apply any StdDevObsStyle styling templates
                s=''
                if 'StdDevObsStyle' in m:
                    dis=vincinv(at_lat,at_lng,to_lat,to_lng)
                    enu=np.zeros([3, 3])
                    s = obs_styl[m['StdDevObsStyle']]
                    enu[0,0]=((s['HzConstant'])+dis[0]*(s['HzPPM'])*10E-6)**2
                    enu[1,1]=((s['HzConstant'])+dis[0]*(s['HzPPM'])*10E-6)**2
                    enu[2,2]=((s['VtConstant'])+dis[0]*(s['VtPPM'])*10E-6)**2
                    enu  = enu / float(m['Vscale'])
                    enu[0,0] = enu[0,0] / float(m['Lscale'])
                    enu[1,1] = enu[1,1] / float(m['Pscale'])
                    enu[2,2] = enu[2,2] / float(m['Hscale'])

                    if 'IncreaseReduce' in s:
                        if (s['IncreaseOrReduce'] == 'Increase'
                            and enu[0,0] > orig_enu[0,0] 
                            and enu[1,1] > orig_enu[1,1] 
                            and enu[2,2] > orig_enu[2,2]):
                                cal_enu = enu + at_sd + to_sd
                        elif (s['IncreaseOrReduce'] == 'Reduce'
                            and enu[0,0] < orig_enu[0,0] 
                            and enu[1,1] < orig_enu[1,1] 
                            and enu[2,2] < orig_enu[2,2]):
                                cal_enu = enu + at_sd + to_sd
                    
                if cal_enu.any != orig_enu.any:
	                #transform and apply to output file
	                vcv = vcv_local2cart(cal_enu, at_lat, at_lng)
	                new_vcv = (
	                    '\t\t\t<SigmaXX>'+ str(vcv[0,0]) + '</SigmaXX>\n' +
	                    '\t\t\t<SigmaXY>'+ str(vcv[0,1]) + '</SigmaXY>\n' +
	                    '\t\t\t<SigmaXZ>'+ str(vcv[0,2]) + '</SigmaXZ>\n' +
	                    '\t\t\t<SigmaYY>'+ str(vcv[1,1]) + '</SigmaYY>\n' +
	                    '\t\t\t<SigmaYZ>'+ str(vcv[1,2]) + '</SigmaYZ>\n' +
	                    '\t\t\t<SigmaZZ>'+ str(vcv[2,2]) + '</SigmaZZ>\n'
	                        )
	
	                msr_out = (msr_out
                            .replace('<Sigma','<!--Sigma')
                            .replace('</SigmaXX>','</SigmaXX-->')
                            .replace('</SigmaXY>','</SigmaXY-->')
                            .replace('</SigmaXZ>','</SigmaXZ-->')
                            .replace('</SigmaYY>','</SigmaYY-->')
                            .replace('</SigmaYZ>','</SigmaYZ-->')
                            .replace('</SigmaZZ>','</SigmaZZ-->')
                            .replace('</GPSBaseline>',
                                 '<!-- VCV Rescaled with Std Dev Styles-->\n'
                                 + new_vcv + '\t\t</GPSBaseline>' ))

        f_out.write(msr_out+'\n')
        i+=1  
