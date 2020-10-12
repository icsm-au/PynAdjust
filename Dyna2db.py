import os, sqlite3
import xmltodict

def delete_DynaML_db(db):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    cursor = conn.cursor()
    cursor.execute('DROP TABLE IF EXISTS FILES')
    cursor.execute('DROP TABLE IF EXISTS DNA_STATION')
    cursor.execute('DROP TABLE IF EXISTS DNA_MEASUREMENT')
    cursor.execute('DROP TABLE IF EXISTS DNA_BASELINES')
    cursor.execute('DROP TABLE IF EXISTS DNA_DIRECTIONS')
    cursor.execute('DROP TABLE IF EXISTS ADJ_MEASUREMENT')
    cursor.execute('DROP TABLE IF EXISTS X_G_Y_COMPONENTS')
    cursor.execute('DROP TABLE IF EXISTS APU_STATION')
    conn.commit()
    conn.close()

def create_DynaML_db(db):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    cursor.execute('''CREATE TABLE IF NOT EXISTS FILES (
                  ID integer PRIMARY KEY);''')
    cursor.execute('''CREATE TABLE IF NOT EXISTS DNA_STATION (
                  ID integer PRIMARY KEY);''')
    cursor.execute('''CREATE TABLE IF NOT EXISTS DNA_MEASUREMENT (
                  ID integer PRIMARY KEY);''')
    cursor.execute('''CREATE TABLE IF NOT EXISTS DNA_BASELINES (
                  ID integer PRIMARY KEY);''')
    cursor.execute('''CREATE TABLE IF NOT EXISTS DNA_DIRECTIONS (
                  ID integer PRIMARY KEY);''')
    cursor.execute('''CREATE TABLE IF NOT EXISTS ADJ_MEASUREMENT (
                  ID integer PRIMARY KEY);''')
    cursor.execute('''CREATE TABLE IF NOT EXISTS X_G_Y_COMPONENTS (
                  ID integer PRIMARY KEY);''')
    cursor.execute('''CREATE TABLE IF NOT EXISTS APU_STATION (
                  ID integer PRIMARY KEY);''')

    conn.commit()
    conn.close()


def add_tbl_clm(tbl,clms,conn):
    cursor = conn.cursor()
    e_clms=cursor.execute('PRAGMA table_info('+tbl+')').fetchall()
    for c in clms:
        add_c=False
        for e in e_clms:
            if c in e: add_c=True
        if add_c==False:
            sql='ALTER TABLE ' + tbl + ' ADD '
            sql= sql + c
            if (c =='M' 
                or c.startswith('Station') 
                or c =='C' or c =='Outlier?' 
                or c =='Const' 
                or c == 'Description'):
                sql = sql + ' short text;'
            else:
                sql = sql + ' double;'
            cursor.execute(sql)
            conn.commit()
        
        
def fmt4sql(stg):
    stg=(stg.replace('*',' ')
            .replace('?',' ')
            .replace('.','')
            .replace('-','_')
            .replace('(','_')
            .replace(')',' ')
            .replace('M Station 1','M  Station 1')
            .upper())
    clms=[s.strip().replace(' ','_') for s in stg.split('  ') if s]
    return clms


def list2sql (clms):
    sql = ' ('
    v=''
    for c in clms:
        sql=sql + c + ', '
        v = v + '?,'
    sql = sql[:-2] + ') VALUES ('+v[:-1] +')'     
    return sql



def add_file_ref(f,tbl,conn):
    cursor = conn.cursor()
    add_tbl_clm(tbl,['FILES_ID'],conn)
    add_tbl_clm('FILES',['FILE_NAME'],conn)
    fle_cnt= cursor.execute("SELECT MAX(ID) FROM FILES").fetchall()[0][0]
    if fle_cnt==None: fle_cnt=0
    cursor.execute('UPDATE '+tbl+ \
                   ' set FILES_ID='+ str(fle_cnt+1) + \
                   ' WHERE FILES_ID is null')
    cursor.execute('''INSERT OR IGNORE INTO FILES (FILE_NAME)
                    VALUES (?)''', [f])
    conn.commit()    
    
def stn_xml2db (f,db):  
    ######## Open the station file 
    ##         create table columns
    ##         add all the station fields to one table
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    with open(f, 'r') as x:
        xml = xmltodict.parse(x.read())

    clms=[]
    for stn in xml['DnaXmlFormat']['DnaStation']:
        dta_stn=[]
        stn.update(stn['StationCoord'])
        for k in stn:
            if k!='StationCoord':
                dta_stn=dta_stn+[stn[k]]
                if k not in clms:
                    add_tbl_clm('DNA_STATION',[k],conn)
                    clms=clms+[k]
        sql =list2sql(clms)
        cursor.execute('INSERT INTO DNA_STATION '+sql, dta_stn)
    conn.commit()
    # Add the file name as a reference
    add_file_ref(f,'DNA_MEASUREMENT',conn)
    conn.close()
    

def msr_xml2db(f,db):
    ######## Open the msr file 
    ##         create table columns
    ##         every measurement is given a measurement number
    ##         every measurement is added to the DNA_MEASUREMENT table
    ##         the components of Directions and GNSS baselines are added to the DNA_DIRECTIONS and DNA_BASELINES tables
    with open(f, 'r') as x:
        xml = xmltodict.parse(x.read())

    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    msr_cnt= cursor.execute("SELECT MAX(ID) FROM DNA_MEASUREMENT").fetchall()[0][0]
    if msr_cnt==None: msr_cnt=0
    add_tbl_clm('DNA_BASELINES',['MEASUREMENT_ID'],conn)
    add_tbl_clm('DNA_DIRECTIONS',['MEASUREMENT_ID'],conn)
    for msr in xml['DnaXmlFormat']['DnaMeasurement']:
        clms=[]
        dta=[]
        msr_cnt+=1
        for k1 in msr:
            if k1=='Directions':
                add_tbl_clm('DNA_DIRECTIONS',list(msr[k1].keys()),conn)
                sql =list2sql(['MEASUREMENT_ID']+list(msr[k1].keys()))
                cursor.execute('INSERT INTO DNA_DIRECTIONS '+
                               sql, [msr_cnt]+list(msr[k1].values()))
            
            elif k1=='GPSBaseline':
                add_tbl_clm('DNA_BASELINES',list(msr[k1].keys()),conn)
                sql =list2sql(['MEASUREMENT_ID']+list(msr[k1].keys()))
                cursor.execute('INSERT INTO DNA_BASELINES '+
                               sql, [msr_cnt]+list(msr[k1].values()))
                
            else:
                if k1 not in clms: 
                    add_tbl_clm('DNA_MEASUREMENT',[k1],conn)
                    clms=clms+[k1]
                dta=dta + [msr[k1]]
        sql =list2sql(clms)
        cursor.execute('INSERT INTO DNA_MEASUREMENT '+
                       sql, dta)
        
    conn.commit()
    # Add the file name as a reference
    add_file_ref(f,'DNA_MEASUREMENT',conn)
    conn.close()

def import_adj(f,db):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    obs_id = cursor.execute("SELECT MAX(ID) FROM ADJ_MEASUREMENT").fetchall()[0][0]
    if obs_id==None: obs_id=0
    obs_id=obs_id+1
    gnss_p=[]
    c_ln=-1
    with open(f, 'r') as adj_f:    
        for ln in adj_f.readlines():
            if ln.strip()=='/n' or ln.strip()=='': continue
            if ln.startswith('Adjusted Coordinates'): c_ln=-1
            if ln.startswith('Station             Const      Latitude      Longitude   H(Ortho) h(Ellipse)'):
                clms=fmt4sql(ln)
                add_tbl_clm('ADJ' ,clms,conn)
                c_ln=2
            if c_ln == 0 and clms[1]=='CONST':
                dta=[(ln[:20].strip())]
                dta = dta + ln[20:].split()
                sql =list2sql(['FILE_ID'] + clms[:len(dta)-1])
                cursor.execute('INSERT INTO ADJ '+sql, dta)                    
            
            if ln.startswith('M Station 1           Station 2           Station 3'):
                clms=fmt4sql(ln)
                add_tbl_clm('ADJ_MEASUREMENT' ,clms,conn)
                add_tbl_clm('X_G_Y_COMPONENTS' ,fmt4sql(ln[65:]),conn)
                c_ln=2
            if c_ln == 0 and clms[0]=='M':
                dta=[ln[:2].strip()]
                dta.append(ln[2:20].strip())    #First
                dta.append(ln[22:40].strip())   #second
                dta.append(ln[42:60].strip())   #third
                dta.append(ln[62:67].strip())   #C
                for e in [s.strip() for s in ln[68:].split('  ') if s]:
                    try: dta.append(float(e))
                    except ValueError: dta.append(e)
                if dta[0] == 'X' or dta[0] == 'Y' or dta[0] == 'G':
                    xgy=[obs_id]+ln[61:].split()                   
                    add_tbl_clm('X_G_Y_COMPONENTS',['OBSERVATION_ID'],conn)
                    sql =list2sql(['OBSERVATION_ID']+clms[4:len(xgy)+3])
                    cursor.execute('INSERT INTO X_G_Y_COMPONENTS '+sql, xgy)
                    gnss_p.append(dta)
                    if len(gnss_p)==3:
                        dta[4]='Avg'
                        for c in range(5,7):
                            dta[c]=''
                        for c in range(7,len(gnss_p[0])-2):
                            avg_xgy=sum(row[c] for row in gnss_p)/3
                            dta[c]=round(avg_xgy,4)
                        gnss_p=[]

                if len(gnss_p)==0:
                    sql = list2sql(clms[:len(dta)])
                    cursor.execute('INSERT INTO ADJ_MEASUREMENT '+sql, dta)                    
                
            if c_ln > 0: c_ln-=1
    conn.commit()
    # Add the file name as a reference
    add_file_ref(f,'ADJ_MEASUREMENT',conn)
    conn.close()

def import_apu(f,db):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    c_ln=-1
    with open(f, 'r') as apu_f:    
        for ln in apu_f.readlines():
            if ln=='/n': c_ln=-1
            if ln.startswith('Station                     Latitude      Longitude    Hz PosU    Vt PosU   Semi-major   Semi-minor  Orientation        Variance(X)        Variance(Y)        Variance(Z)'):
                clms=fmt4sql(ln)
                add_tbl_clm('APU_STATION' ,clms,conn)
                c_ln=2
            if c_ln == 0:
                stn = [ln[:20].strip()]         #Station
                stn = stn + ln[20:].split()     #Latitude...  etc
                sql =list2sql(clms)
                cursor.execute('INSERT INTO APU_STATION '+sql, stn)
                c_ln=3
            if c_ln > 0: c_ln-=1
    conn.commit()

    # Add the file name as a reference
    add_file_ref(f,'APU_STATION',conn)           
    conn.close()

if __name__ == "__main__":
    networks =[s.replace('.stn.xml','') for s in os.listdir('.') if s.endswith('.stn.xml')]
    ######### Create a database and empty tables ############
    for network in networks:
        create_DynaML_db(network+'.db')
        
        if os.path.exists(network+'.stn.xml'):
            stn_xml2db (network+'.stn.xml',network+'.db')
      
        if os.path.exists(network+'.msr.xml'):
            msr_xml2db (network+'.msr.xml',network+'.db')
    
        if os.path.exists(network+'.phased-stage.adj'):
            import_adj (network+'.phased-stage.adj',network+'.db')

        if os.path.exists(network+'.phased-stage.apu'):
            import_apu (network+'.phased-stage.apu',network+'.db')    
