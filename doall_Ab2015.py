from astropy.table import Table
import numpy as np
from astroquery.utils.tap.core import TapPlus

def doquery_from_table(table):
    #Execute a TAP+ query using information from a single row of an astropy table
    server = TapPlus(url = table['queryurl'], \
                     default_protocol_is_https = True, \
                     verbose = True)
    if table['upload_table_name'] == '':
        result = server.launch_job_async(query = table['query'], verbose = True, \
                                         dump_to_file = True, output_file = table['outfile'])
    else:
        result = server.launch_job_async(query = table['query'], verbose = True, \
                                         dump_to_file = True, output_file = table['outfile'], \
                                         upload_table_name = table['upload_table_name'], \
                                         upload_resource = table['upload_resource'])
    return result

def store_queries(querytable = 'queries.xml'):
    #for coloured messages
    class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        def disable(self):
            self.HEADER = ''
            self.OKBLUE = ''
            self.OKGREEN = ''
            self.WARNING = ''
            self.FAIL = ''
            self.ENDC = ''
    #The default TAP server url is the one for Vizier.
    queryurl = "http://tapvizier.u-strasbg.fr/TAPVizieR/tap"

    #query records
    url1 = queryurl
    ur1 = ''
    un1 = ''
    out1 = 'Ab2015.xml'
    q1 = """
    --Abrahamyan et al., 2015
    select a."IRAS-PSC" as IRAS_PSC_ID, "IRAS-FSC" as IRAS_FSC_ID, "Akari-IRC" as AKARI_IRC_ID, "Akari-FIS" as AKARI_FIS_ID,
    WISE as WISE_AllSky_ID, RApsc as IRAS_PSC_RA, DEpsc as IRAS_PSC_DEC, pscMaj as IRAS_PSC_e_Maj, pscMin as IRAS_PSC_e_Min, pscPA as IRAS_PSC_e_PA,
    RAfsc as IRAS_FSC_RA, DEfsc as IRAS_FSC_DEC, fscMaj as IRAS_FSC_e_Maj, fscMin as IRAS_FSC_e_Min, fscPA as IRAS_FSC_e_PA,
    RAir as IRAS_best_pos_RA, e_RAir as IRAS_best_pos_e_RA, DEir as IRAS_best_pos_DEC, e_DEir as IRAS_best_pos_e_DEC,
    irMaj as IRAS_best_pos_e_Maj, irMin as IRAS_best_pos_e_Min, irPA as IRAS_best_pos_e_PA, "P/F" as IRAS_comb_coord_flag,
    "d(ir)" as IRAS_PSC_FSC_dist, "k(ir)" as IRAS_PSC_FSC_rel_dist, "d(irc)" as IRAS_AKARI_IRC_dist, "k(irc)" as IRAS_AKARI_IRC_rel_dist,
    "d(fis)" as IRAS_AKARI_FIS_dist, "k(fis)" as IRAS_AKARI_FIS_rel_dist, "d(wis)" as IRAS_AKARI_WISE_AllSky_dist,
    flag as IRAS_AKARI_WISE_best_coord_flag
    from "II/338/catalog" a
    """
    
    url2 = queryurl
    ur2 = ''
    un2 = ''
    out2 = 'Ab2015_x_IRASPSC.xml'
    q2 = """
    --Abrahamyan et al., 2015 x IRAS PSC
    --NOTE THAT THIS IS AN INNER JOIN
    select a."IRAS-PSC" as IRAS_PSC_ID, a."IRAS-FSC" as IRAS_FSC_ID,
    b.Fnu_12 as IRAS_PSC_Fnu_12, b.e_Fnu_12 as IRAS_PSC_e_Fnu_12, b.q_Fnu_12 as IRAS_PSC_q_Fnu_12, b.Fnu_25 as IRAS_PSC_Fnu_25,
    b.e_Fnu_25 as IRAS_PSC_e_Fnu_25, b.q_Fnu_25 as IRAS_PSC_q_Fnu_25, b.LRSChar as IRAS_PSC_LRS_Char
    from "II/338/catalog" a join "II/125/main" b
    on a."IRAS-PSC" = b.IRAS
    """
    
    url3 = queryurl
    ur3 = ''
    un3 = ''
    out3 = 'Ab2015_x_IRAS.xml'
    q3 = """
    --Abrahamyan et al., 2015 x IRAS FSC
    --NOTE THAT THIS IS AN INNER JOIN
    select a."IRAS-PSC" as IRAS_PSC_ID, a."IRAS-FSC" as IRAS_FSC_ID,
    b.Fnu12 as IRAS_FSC_Fnu12, b.e_Fnu12 as IRAS_FSC_e_Fnu12, b.q_Fnu12 as IRAS_FSC_q_Fnu12, b.Fnu25 as IRAS_FSC_Fnu25,
    b.e_Fnu25 as IRAS_FSC_e_Fnu25, b.q_Fnu25 as IRAS_FSC_q_Fnu25
    from "II/338/catalog" a join "II/156A/main" b
    on a."IRAS-FSC" = b.IRAS
    """
    
    url4 = queryurl
    ur4 = ''
    un4 = ''
    out4 = 'Ab2015_IRAS_x_AKARIIRC.xml'
    q4 = """
    --Abrahamyan et al., 2015 x AKARI IRC
    --NOTE THAT THIS IS AN INNER JOIN
    select a."IRAS-PSC" as IRAS_PSC_ID, a."IRAS-FSC" as IRAS_FSC_ID,
    b.RAJ2000 as AKARI_IRC_RA, b.DEJ2000 as AKARI_IRC_DEC, b.S09 as AKARI_IRC_S09, b.e_S09 as AKARI_IRC_e_S09, b.q_S09 as AKARI_IRC_q_S09,
    b.S18 as AKARI_IRC_S18, b.e_S18 as AKARI_IRC_e_S18, b.q_S18 as AKARI_IRC_q_S18
    from "II/338/catalog" a join "II/297/irc" b
    on a."Akari-IRC" = b.objname
    """
    
    url5 = queryurl
    ur5 = ''
    un5 = ''
    out5 = 'Ab2015_IRAS_AKARI.xml'
    q5 = """
    --Abrahamyan et al., 2015 x AKARI IRC
    --NOTE THAT THIS IS AN INNER JOIN
    select a."IRAS-PSC" as IRAS_PSC_ID, a."IRAS-FSC" as IRAS_FSC_ID,
    b.RAJ2000 as AKARI_FIS_RA, b.DEJ2000 as AKARI_FIS_DEC, b.S65 as AKARI_FIS_S65, b.e_S65 as AKARI_FIS_e_S65, b.q_S65 as AKARI_FIS_q_S65,
    b.S90 as AKARI_FIS_S90, b.e_S90 as AKARI_FIS_e_S90, b.q_S90 as AKARI_FIS_q_S90,
    b.S140 as AKARI_FIS_S140, b.e_S140 as AKARI_FIS_e_S140, b.q_S140 as AKARI_FIS_q_S140,
    b.S160 as AKARI_FIS_S160, b.e_S160 as AKARI_FIS_e_S160, b.q_S160 as AKARI_FIS_q_S160
    from "II/338/catalog" a join "II/298/fis" b
    on a."Akari-FIS" = b.objNAME
    """
    
    url6 = "http://dc.zah.uni-heidelberg.de/__system__/tap/run/tap"
    ur6 = 'Ab2015_PSCFSCWISEID.xml'
    un6 = ur6[:-4]
    out6 = 'Ab2015_IRAS_AKARI_x_WISEAllSky.xml'
    q6 = """
    --This query requires the generation of Ab2015_PSCFSCWISEID.xml, which is done in runAb2015queries below.
    --Abrahamyan et al. 2015 x WISE All-Sky
    --NOTE THAT THIS IS AN INNER JOIN
    select a.IRAS_PSC_ID, a.IRAS_FSC_ID,
    b.RAJ2000 as WISE_AllSky_RA, b.DEJ2000 as WISE_AllSky_DEC,
    b.sigra as WISE_AllSky_E_RA, b.sigdec as WISE_AllSky_E_DEC,
    b.w1mpro as WISE_AllSky_W1_MPRO, b.w1sigmrpo as WISE_AllSky_W1_SIGMPRO, b.w1mag as WISE_AllSky_W1_MAG, b.w1sigm as WISE_AllSky_W1_SIGM,
    b.w2mpro as WISE_AllSky_W2_MPRO, b.w2sigmrpo as WISE_AllSky_W2_SIGMPRO, b.w2mag as WISE_AllSky_W2_MAG, b.w2sigm as WISE_AllSky_W2_SIGM,
    b.w3mpro as WISE_AllSky_W3_MPRO, b.w3sigmrpo as WISE_AllSky_W3_SIGMPRO, b.w3mag as WISE_AllSky_W3_MAG, b.w3sigm as WISE_AllSky_W3_SIGM,
    b.w4mpro as WISE_AllSky_W4_MPRO, b.w4sigmrpo as WISE_AllSky_W4_SIGMPRO, b.w4mag as WISE_AllSky_W4_MAG, b.w4sigm as WISE_AllSky_W4_SIGM,
    b.tmass_key as WISE_AllSky_2MASS_KEY
    from tap_upload.""" + un6.strip() + """ a join wise.main b
    on a.WISE_AllSky_ID = b.designation
    """
    
    url7 = "http://dc.zah.uni-heidelberg.de/__system__/tap/run/tap"
    ur7 = 'Ab2015_PSCFSCWISE2MASSID.xml'
    un7 = ur7[:-4]
    out7 = 'Ab2015_IRAS_AKARI_WISEAllSky_x_2MASS.xml'
    q7 = """
    --Abrahamyan et al. 2015 x WISE All-Sky x 2MASS (2MASS counterparts to WISE All-Sky counterparts to Ab2015)
    --NOTE THAT THIS IS AN INNER JOIN
    select a.IRAS_PSC_ID, a.IRAS_FSC_ID,
    b.raj2000 as WISE_AllSky_2MASS_RA, b.dej2000 as WISE_AllSky_2MASS_DEC, b.errmaj as WISE_AllSky_2MASS_ERR_MAJ, 
    b.errmin as WISE_AllSky_2MASS_ERR_MIN, b.errpa as WISE_AllSky_2MASS_ERR_PA, b.mainid as WISE_AllSky_2MASS_DESIGNATION, 
    b.jmag as WISE_AllSky_2MASS_JMAG, b.e_jmag as WISE_AllSky_2MASS_E_JMAG, b.hmag as WISE_AllSky_2MASS_HMAG, 
    b.e_hmag as WISE_AllSky_2MASS_E_HMAG, b.kmag as WISE_AllSky_2MASS_KMAG, b.e_kmag as WISE_AllSky_2MASS_E_KMAG
    from tap_upload. """ + un7.strip()+ """ a, b
    where a.WISE_AllSky_2MASS_KEY = b.PTS_KEY
    """
    
    #Combine everything to generate data columns for the astropy table
    q = [q1, q2, q3, q4, q5, q6, q7]
    ur = [ur1, ur2, ur3, ur4, ur5, ur6, ur7]
    un = [un1, un2, un3, un4, un5, un6, un7]
    url = [url1, url2, url3, url4, url5, url6, url7]
    out = [out1, out2, out3, out4, out5, out6, out7]

    t = Table([url, q, ur, un, out], names = ('queryurl', 'query', \
                                              'upload_resource', 'upload_table_name', 'outfile'), \
              meta = {'ID': 'Queries using Abrahamyan et al. 2015'})
    t.write('Queries_Ab2015.xml', overwrite = True, format = 'votable')

def runAb2015queries():
    t = Table.read('Queries_Ab2015.xml', format = 'votable')
    r = []
    for i in range(len(t)):
        print(bcolors.HEADER + "-------------------------------------------------" + bcolors.ENDC)
        print(bcolors.HEADER + "------Running query #{}--------------------------".format(i+1) + bcolors.ENDC)
        try:
            #required for query 6
            if i == 5:
                t = Table.read('Ab2015.xml', format = 'votable')
                s = Table([t['IRAS_PSC_ID'], t['IRAS_FSC_ID'], t['WISE_AllSky_ID']])
                s.write('Ab2015_PSCFSCWISEID.xml', format = 'votable', overwrite = True)
                print(bcolors.WARNING + "This may take a while. Go do something else." + bcolors.ENDC)
            #required for query 7
            if i == 6:
                t = Table.read('Ab2015_IRAS_AKARI_x_WISEAllSky.xml', format = 'votable')
                s = Table([t['IRAS_PSC_ID'], t['IRAS_FSC_ID'], t['WISE_AllSky_2MASS_KEY']])
                s.write('Ab2015_PSCFSCWISE2MASSID.xml', overwrite = True, format = 'votable')
                print(bcolors.WARNING + "This may take a while. Go do something else." + bcolors.ENDC)
            #rr = doquery2(query = t['query'][i], outfile = t['outfile'][i])
            rr = doquery_from_table(t[i])
            r.append(rr)
            print(bcolors.OKGREEN + "------End query #{}--------------------------".format(i+1) + bcolors.ENDC)
            print(bcolors.OKGREEN + "1111111111111111111111111111111111111111111111111" + bcolors.ENDC)
        except:
            print(bcolors.FAIL + "Uh oh, something went wrong. Continuing to next query" + bcolors.ENDC)
            print(bcolors.FAIL + "0000000000000000000000000000000000000000000000000" + bcolors.ENDC)

def combinetables():
    #combine the results of Ab2015queries into one giant table.
    t1 = Table.read('Ab2015.xml', format = 'votable')
    t2 = Table.read('Ab2015_x_IRASPSC.xml', format = 'votable')
    tt = join()

