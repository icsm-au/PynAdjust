# -*- coding: utf-8 -*-

import os
from geodepy.geodesy import vincinv
from geodepy.convert import xyz2llh, hp2dec
from math import tan, radians, sqrt

adj = [ f for f in os.listdir(os.path.join(os.getcwd(),'adjustments')) if f.endswith('simult.adj')]

geojson = open('HzAnalysis2.geojson', 'w')
geojson.write('{\n')
geojson.write('"type": "FeatureCollection",\n')
geojson.write('"name": "NGCA_HzAnalysis",\n')
geojson.write('"crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::4019" } },\n')
geojson.write('"features": [\n')
geojson_str =''  

  
for f in adj:
    ###########################################################################
    #        Open the adj file and find the best apref station                #
    ###########################################################################
    best_nstat = 100000000000000
    part_G_cnt=1; o = 1; c = 1
    angular_msr_type=''; coords = {}; xyz_bases = {}; apref={}; ngca = {}
    with open(os.path.join('adjustments',f), 'r') as f_in:
        for l in f_in.readlines():
            #Read the baselines
            if (l.find('Command line arguments:')!=-1 and 
                l.find('--angular-msr-type 1')!=-1):angular_msr_type='dd'
            if l.find('Station coordinate types:')!=-1: coord_types = l[35:].strip()
            if l.find('Adjusted Measurements')!=-1: o = -4
            if o == 0 and l.startswith('X '):
                a_obs=l[67:].split()
                if part_G_cnt==1: dX=float(a_obs[0]); nstat_x=float(a_obs[6])
                if part_G_cnt==2: dY=float(a_obs[0]); nstat_y=float(a_obs[6])
                if part_G_cnt==3: 
                    dZ=float(a_obs[0]); nstat_z=float(a_obs[6])
                    first_stn = l[2:22].strip()
                    xyz_bases[l[22:42].strip()] = {
                             'dX':dX, 'dY':dY, 'dZ':dZ,
                             'nstat_x':nstat_x, 'nstat_y':nstat_y, 'nstat_z':nstat_z}
                    part_G_cnt=0
                part_G_cnt+=1
            if o <0: o+=1
            #Read the Coordinates
            if l.find('Adjusted Coordinates')!=-1:
                c=-5; o=1
            if c==0 and l!='\n':
                stn = l[0:20].strip()
                results = l[25:].split()
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
                    if ct == 'X': X = float(results[r_count])
                    if ct == 'Y': Y = float(results[r_count])
                    if ct == 'Z': Z = float(results[r_count])

                    r_count += 1
                
                if l[20:25].strip() == 'CCC':
                    apref[str(stn)] = {'P': P,'L': L,'OHGT': H,'EHGT': h,
                                      'X': X,'Y': Y,'Z': Z}
                if str(stn)==first_stn:
                    for b in xyz_bases:
                        x = X + float(xyz_bases[b]['dX'])
                        y = Y + float(xyz_bases[b]['dY'])
                        z = Z + float(xyz_bases[b]['dZ'])
                        P,L,h = xyz2llh(x,y,z)
                        ngca[str(b)] = {'P': P,'L': L,'EHGT': h,
                                         'X': x, 'Y': y, 'Z': z}
            if c<0:c+=1

    ###########################################################################
    #               Analyse vector agreement of NGCA Vs APREF                 #
    ###########################################################################
    best_sum = 0
    prev_best = 1000000000000000000
    all_joins = []
    for s1 in apref:
        for s2 in apref:
            if s1!=s2:
                NGCA_dis,NGCA_az,RvAz = vincinv(ngca[s1]['P'],ngca[s1]['L'],
                                                ngca[s2]['P'],ngca[s2]['L'])
                APREF_dis,APREF_az,RvAz = vincinv(apref[s1]['P'],apref[s1]['L'],
                                                  apref[s2]['P'],apref[s2]['L'])
                dif_dis = NGCA_dis-APREF_dis
                dif_Az  = tan(radians(NGCA_az-APREF_az))*APREF_dis
                dif_ht = ((apref[s2]['EHGT'] - apref[s1]['EHGT']) 
                        - (ngca[s2]['EHGT'] - ngca[s1]['EHGT']))
                best_sum = best_sum + dif_dis
                all_joins = all_joins + [[
                        s1,s2,apref[s1]['P'],apref[s1]['L'],
                        apref[s2]['P'],apref[s2]['L']
                        ,round(dif_dis,4),round(dif_Az,4),round(dif_Az,4)]]
        if best_sum < prev_best:
            prev_best = best_sum
            best_mark = s1
    print('Best Mark: ' + best_mark)
    best_joins = [ i for i in all_joins if i[0]==best_mark]
    
    for j in best_joins:
        geojson_str = geojson_str + ('{ "type": "Feature", "properties": { ' +
        '"CLUSTER": "' + f.replace('.SNX.AUS','') + '", ' +
        '"FROM": "' + str(j[0]) + '", ' +
        '"TO": "' + str(j[1]) + '", ' +
        '"DIST(m)": ' + str(j[6]) + ', ' +
        '"ROTATION(m)": ' + str(j[7]) + ', ' +
        '"Height(m)": ' + str(j[8]) + ', ' + 
        '"Combined(m)": ' + str(round(sqrt(j[6]**2+j[7]**2+j[8]**2),4)) + ' ' + '}, '+
        '"geometry": { "type": "LineString", "coordinates": [ ' +
        '[' + str(j[3]) + ', ' + str(j[2]) + '],' +
        '[' + str(j[5]) + ', ' + str(j[4]) + '] ] } }' +
        ',\n')
geojson.write(geojson_str[:-2] + '\n')
geojson.write(']\n')
geojson.write('}\n')
geojson.close()
