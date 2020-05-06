from astropy.table import Table, join, vstack
import numpy as np
import pandas as pd
import os, subprocess, time
from astroquery.utils.tap.core import TapPlus
import requests as r
from os.path import isfile

def doquery_from_table(t):
    table = t[0]
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

def combineVOTables(infiles = [], outfile = ''):
    #Combine a list of VOTable files into a single VOTable file using vstack
    nfiles = len(infiles)
    if nfiles == 1:
        print("Only one file in list. Nothing to be done. Exiting.")
        return 0
    else:
        tlist = []
        for i in range(nfiles):
            tlist.append(Table.read(infiles[i], format = 'votable'))
        vstack(tlist).write(outfile, format = 'votable', overwrite = True)

def splitVOTable(intable = None, infile = '', n_outfiles = None):
    #Split an input VOTable into n_outfiles output files
    #   If intable is set, then astropy table from that variable is
    #   split-written over n_outfiles files.
    #   infile is a necessary input, as it is used to determine
    #   output file names.
    if (n_outfiles == 1) or (n_outfiles is None):
        print("splitVOTable: Number of output files either not set or set to 1. Exiting.")
        return 0
    else:
        if intable is None:
            t = Table.read(infile, format = 'votable')
        else:
            t = intable
        if infile == '':
            print("splitVOTable: Name of input file must be specified. Exiting.")
            return 0
        else:
            name, ext = tuple(infile.split('.'))
            outfiles = [name + '_part' + str(j + 1) + '.' + ext for j in range(n_outfiles)]
            chunk = int(np.round(len(t)/n_outfiles))
            for i in range(n_outfiles - 1):
                t[i*chunk: (i+1)*chunk].write(outfiles[i], format = 'votable', overwrite = True)
            t[(i+1)*chunk:].write(outfiles[i+1], format = 'votable', overwrite = True)

def getAb2015():
    #Download the Abrahamyan+ 2015 catalog
    query = """select * from "II/338/catalog" """
    queryurl = 'http://tapvizier.u-strasbg.fr/TAPVizieR/tap'
    upload_resource = ''
    upload_table_name = ''
    outfile = 'Ab2015.xml'
    q = Table([[queryurl], [query], [upload_resource], [upload_table_name], [outfile]], \
              names = ('queryurl', 'query', 'upload_resource', 'upload_table_name', 'outfile'))
    r = doquery_from_table(q)
    #Change column name 'recno' to 'IRAS_CNTR' and sort by 'IRAS_CNTR'
    t = Table.read('Ab2015.xml', format = 'votable')
    t.rename_column('recno', 'IRAS_CNTR')
    k = np.argsort(t['IRAS_CNTR'])
    #Now, set the matching radius according to Abrahamyan+ 2015:
    #3" if AKARI-IRC, 15" if AKARI-FIS, 20" if IRAS PSC/FSC combined, 30" if either FSC or PSC.
    WISE_MATCH_RADIUS = np.repeat(30., len(t))
    WISE_MATCH_RADIUS[t['flag'] >= 5] = 3.0
    WISE_MATCH_RADIUS[t['flag'] == 4] = 15.0
    WISE_MATCH_RADIUS[t['flag'] == 3] = 20.0
    WISE_MATCH_RADIUS[t['flag'] <= 2] = 30.0
    t['WISE_MATCH_RADIUS'] = WISE_MATCH_RADIUS
    #Overwrite
    t[k].write('Ab2015.xml', format = 'votable', overwrite = True)

def bestcoord_IRAS_AKARI():
    """Starting with the Ab2015.xml table, store the best coordinate pair
    determined hierarchically ONLY FROM IRAS and AKARI data thusly:
    AKARI IRC > AKARI FIS > IRAS combined > IRAS PSC > IRAS FSC.
    These coordinates are used to match to the AllWISE database.
    """
    t = Table.read('Ab2015.xml', format = 'votable')

    ##First, generate a table with only IRAS PSC/FSC best coordinates
    s = Table([t['IRAS_CNTR'], t['RAir'], t['DEir']])
    #Split into pieces to allow resulting match tables to be manageable in size
    splitVOTable(intable = s, infile = 'Ab2015_IRAS_coord.xml', n_outfiles = 3)
    #s[:int(np.round(len(s)/3))].write('Ab2015_IRAS_coord_part1.xml', \
    #                                  format = 'votable', overwrite = True)
    #s[int(np.round(len(s)/3)):2*int(np.round(len(s)/3))].write('Ab2015_IRAS_coord_part2.xml', \
    #                                  format = 'votable', overwrite = True)
    #s[2*int(np.round(len(s)/3)):].write('Ab2015_IRAS_coord_part3.xml', \
    #                                  format = 'votable', overwrite = True)

    #If best coordinate is WISE, find the next best thing.
    #   If the best coordinate isn't WISE, accept it as is.
    BEST_RA = np.where(t['flag'] != 6, t['RAJ2000'], \
                       np.where(t['Akari-IRC'] != b'', t['RAirc'], \
                                np.where(t['Akari-FIS'] != b'', t['RAfis'], t['RAir'])
                                )
                       )
    BEST_DEC = np.where(t['flag'] != 6, t['DEJ2000'], \
                       np.where(t['Akari-IRC'] != b'', t['DEirc_deg'], \
                                np.where(t['Akari-FIS'] != b'', t['DEfis'], t['DEir'])
                                )
                       )

    #Add BEST_RA, BEST_DEC, and WISE_MATCH_RADIUS to the table
    t['BEST_RA'] = BEST_RA
    t['BEST_DEC'] = BEST_DEC
    #Define new sparse table with only the useful columns
    s = Table([t['IRAS_CNTR'], t['flag'], t['BEST_RA'], t['BEST_DEC']])
    #Split into two pieces to allow resulting match tables to be manageable in size
    splitVOTable(intable = s, infile = 'Ab2015_IRAS_AKARI_bestcoord.xml', n_outfiles = 3)
    #s[:int(np.round(len(s)/3))].write('Ab2015_IRAS_AKARI_bestcoord_part1.xml', \
    #                                  format = 'votable', overwrite = True)
    #s[int(np.round(len(s)/3)):2*int(np.round(len(s)/3))].write('Ab2015_IRAS_AKARI_bestcoord_part2.xml', \
    #                                  format = 'votable', overwrite = True)
    #s[2*int(np.round(len(s)/3)):].write('Ab2015_IRAS_AKARI_bestcoord_part3.xml', \
    #                                  format = 'votable', overwrite = True)

def run_IPAC_query(infiles, outfiles):
    """Given an input table containing the IRAS_PSC_ID and IRAS_FSC_ID and the corresponding
    AllWISE neighbour's unique long-integer ID, find the corresponding TMASS_KEY by querying
    the IPAC TAP server.
    infiles is a list of input file names for each of which the query is executed, and the output
    is stored in the corresponding element of outfiles.
    Upon submission, the program determines the url location ('Location') to check query status.
    When the status is 'COMPLETED', the output download is requested."""
    """NOTE: currently implemented using a call to curl via subprocess. Would be nice to use
    python's requests package, but there seems to be no direct analogue for the --dump-header option,
    which is needed to read in the url location ('Location') to query for output."""

    #First, submit both queries.
    url = [] #will store the location to query for status and output.
    for i in range(len(infiles)):
        filestring = '"table=@' + infiles[i] + '"'
        curl_cmd = ['curl', '--dump-header', 'header', '-F', '"UPLOAD=my_table,param:table"', \
                    '-F', filestring, \
                    #'-F', '"QUERY=SELECT TAP_UPLOAD.my_table.IRAS_PSC_ID as IRAS_PSC_ID, ' + \
                    #'TAP_UPLOAD.my_table.IRAS_FSC_ID as IRAS_FSC_ID, ' + \
                    '-F', '"QUERY=SELECT TAP_UPLOAD.my_table.TEMP_CNTR as TEMP_CNTR, ' + \
                    'allwise_p3as_psd.designation as AllWISE_ID, allwise_p3as_psd.cntr as AllWISE_CNTR, ' + \
                    'allwise_p3as_psd.tmass_key as TMASS_KEY ' + \
                    'FROM allwise_p3as_psd, TAP_UPLOAD.my_table ' + \
                    #'WHERE TAP_UPLOAD.my_table.ID = allwise_p3as_psd.cntr"', \
                    'WHERE TAP_UPLOAD.my_table.ID = allwise_p3as_psd.cntr and allwise_p3as_psd.cntr IS NOT NULL"', \
                    'https://irsa.ipac.caltech.edu/TAP/async']
        #The following doesn't work. Until it's fixed, have to force user to paste command
        #   into the shell.
        #print("Running command: " + ' '.join(curl_cmd))
        #subprocess.run(curl_cmd)
        #replace following three lines with the two preceding lines when fixed.
        print("Run the following command in the shell:")
        print(' '.join(curl_cmd))
        j = input("Press ENTER when done: ")

        with open('header') as f:
            lines = f.readlines()
            print(lines)
        subprocess.run(['rm', 'header'])
        url.append([(line.strip('\n').split(' '))[1] for line in lines if 'Location:' in line][0])
        print("The request is being processed at {}".format(url[i]))

    #Next, wait until queries are COMPLETED and download data.
    #The code checks every ten minutes since these files are huge.
    #In a separate loop because submission takes much less time than execution, so
    #   queries after the first one keep running in the background while we check
    #   for success of the first one.
    for i in range(len(infiles)):
        while True:
            op = subprocess.check_output(['curl', url[i] + '/phase'])
            if op == b'COMPLETED':
                print("Query #{} successful! Downloading output to {}".format(i+1, outfiles[i]))
                subprocess.run(['curl', '-o' + outfiles[i], url[i] + '/results/result'])
                break
            elif b'ERROR' in op or b'ABORTED' in op:
                print("Uh oh, query status: " + op)
                break
            else:
                #check again in 10 min.
                time.sleep(600)

    #Combine the two output tables from the above into an input table to obtain 2MASS information
    #   by retaining only the rows that have non-trivial TMASS_KEY values
    t1 = Table.read(outfiles[0], format = 'votable')
    t2 = Table.read(outfiles[1], format = 'votable')
    t = vstack([t1[~t1['tmass_key'].mask], t2[~t2['tmass_key'].mask]])

    return 0

#def CDSxmatch(table1, table2, radius, outfile):
#    """
#    Format: table names are provided as either 'filename.ext (upload)' or 'CDS_table_reference (CDS/Vizier)'.
#    For example: 'Ab2015_IRAS_AKARI_combined_part1.xml (upload)' or 'II/338/catalog (CDS/Vizier)'
#    """
#    class bcolors:
#        HEADER = '\033[95m'
#        OKBLUE = '\033[94m'
#        OKGREEN = '\033[92m'
#        WARNING = '\033[93m'
#        FAIL = '\033[91m'
#        ENDC = '\033[0m'
#        def disable(self):
#            self.HEADER = ''
#            self.OKBLUE = ''
#            self.OKGREEN = ''
#            self.WARNING = ''
#            self.FAIL = ''
#            self.ENDC = ''
#
#    string = """
#    Go to http://cdsxmatch.u-strasbg.fr/ and perform a cross-match between the following two tables:
#    """ + table1 + """ and """ + table2 + """ with a radius of """ + str(radius) + """ arcsec.
#    Download the output in CSV format into """ + outfile + """."""
#
#    print(bcolors.WARNING + string + bcolors.ENDC)
#
#    input("Press ENTER when done: ")

def cdsxmatch(cat1 = None, cat2 = None, match_radius = 5,
              outformat = 'votable', outfile = 'out', overwrite = True,
              colRA1 = None, colDEC1 = None, colRA2 = None, colDEC2 = None,
              cols1 = None, cols2 = None,
              selection = 'all'):
    """
    DESCRIPTION: Find matches between two catalogues with a specified radius.
    Command-line access to the CDS xmatch utility. This script is based on
    http://cdsxmatch.u-strasbg.fr/xmatch/doc/API-calls.html and Peter Scicluna's
    original interpretation of the description on this webpage.
    INPUTS:
    - cat1 and cat2 can either be strings containing vizier catalog names (e.g., 'vizier:II/328/allwise')
    or file objects for uploads (e.g., open('upload_filename', 'r')).
    - match_radius is the desired match radius in arcsec.
    - outformat and outfile, respectively, specify the format and name of the output file. 
    outformat can be 'votable' or 'csv' or 'ascii'
    The extension for outfile is determined automatically from outformat:
    outformat = 'votable' --> extension = '.vot'
    outformat = 'csv' --> extension = '.csv'
    outformat = 'ascii' --> extension = '.txt'
    - colRA* and colDEC* are strings identifying the columns in the catalogues that contain the RA/DEC.
    This is recommended when uploading tables.
    - cols1 and cols2 are strings consisting of comma-separated names of columns in the two catalogues. If
    specified, only these columns are extracted.
    - selection can either be 'all' (in which case all neighbours within match_radius are extracted) or
    'best' (in which case only the nearest neighbour for each row is extracted).
    """
    #Initialise
    data = {'request': 'xmatch', 'distMaxArcsec': match_radius, 'RESPONSEFORMAT': outformat,
            selection: selection}
    files = {}
    #Depending on whether cat1 and cat2 are strings containing Vizier catalogue names or file objects
    #   containing upload file information, append to one of the two above dictionaries.
    if isinstance(cat1, str):
        data['cat1'] = cat1
    else:
        files['cat1'] = cat1
    if isinstance(cat2, str):
        data['cat2'] = cat2
    else:
        files['cat2'] = cat2
    #Coordinate identification
    if colRA1 is not None:
        data['colRA1'] = colRA1
        data['colDec1'] = colDEC1
    if colRA2 is not None:
        data['colRA2'] = colRA2
        data['colDec2'] = colDEC2
    if cols1 is not None:
        data['cols1'] = cols1
    if cols2 is not None:
        data['cols2'] = cols2

    #
    if files != {}:
        res = r.post('http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync',
                     data = data, files = files)
    else:
        res = r.post('http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync',
                     data = data)

    #Output file and format resolution
    fmt = ['votable', 'csv', 'ascii']
    ext = ['.vot', '.csv', '.txt']
    indx = 0
    try:
        indx = fmt.index(outformat)
    except:
        print("Warning: invalid output file format specified! Saving results to VOTable {}.vot.".\
              format(outfile))
    #Output
    outf = outfile + ext[indx]
    if (isfile(outf)) & ~(overwrite):
        choice = input("Output file exists! Overwrite [Y]/N: ")
        if (choice != '') & (choice.upper() != 'Y'):
            print("Writing to file output" + ext[indx])
            outf = 'output' + ext[indx]
    f = open(outfile + ext[indx], 'w')
    f.write(res.text, overwrite = True)
    f.close()

def getWISEAllSkyAllWISEnnb(radius = 60):
    t = Table.read('Ab2015.xml', format = 'votable')
    k = np.nonzero(t['WISE'] != '')[0]
    s = Table([t[k]['IRAS_CNTR'], t[k]['RAwis'], t[k]['DEwis']])
    splitVOTable(intable = s, infile = 'Ab2015_WISEAllSky.xml', n_outfiles = 3)
    for i in range(3):
        cdsxmatch(cat1 = open('Ab2015_WISEAllSky_part' + str(i+1) + '.xml', 'r'), cat2 = 'vizier:II/328/allwise',
                  match_radius = radius, colRA1 = 'RAwis', colDEC1 = 'DEwis',
                  outfile = 'Ab2015_WISEAllSky_part' + str(i+1) + 'AllWISEnnb_60arcsec', outformat = 'votable', \
                  overwrite = True)
    combineVOTables(infiles = ['Ab2015_WISEAllSky_part' + str(i+1) + 'AllWISEnnb_60arcsec.vot' for i in range(3)],
                    outfile = 'Ab2015_WISEAllSky_AllWISEnnb_60arcsec.vot')
    t = Table.read('Ab2015_WISEAllSky_AllWISEnnb_60arcsec.vot', format = 'votable')
    lk = np.lexsort([t['IRAS_CNTR'], t['angDist']])
    _, lu = np.unique(t[lk]['IRAS_CNTR'], return_index = True)
    t[lk[lu]].write('Ab2015_WISEAllSky_AllWISEnnb_60arcsec.vot', format = 'votable', overwrite = True)

def getWISEAllSkyAllWISEMSXnnb(radius = 60):
    t = Table.read('MSXonly.vot', format = 'votable')
    k = np.nonzero(t['WISE'] != '')[0]
    s = Table([t[k]['IRAS_CNTR'], t[k]['RAwis'], t[k]['DEwis']])
    splitVOTable(intable = s, infile = 'Ab2015_WISEAllSky.xml', n_outfiles = 3)
    for i in range(3):
        cdsxmatch(cat1 = open('Ab2015_WISEAllSky_part' + str(i+1) + '.xml', 'r'), cat2 = 'vizier:II/328/allwise',
                  match_radius = radius, colRA1 = 'RAwis', colDEC1 = 'DEwis',
                  outfile = 'Ab2015_WISEAllSky_part' + str(i+1) + 'AllWISEnnb_60arcsec', outformat = 'votable', \
                  overwrite = True)
    combineVOTables(infiles = ['Ab2015_WISEAllSky_part' + str(i+1) + 'AllWISEnnb_60arcsec.vot' for i in range(3)],
                    outfile = 'Ab2015_WISEAllSky_AllWISEnnb_60arcsec.vot')
    t = Table.read('Ab2015_WISEAllSky_AllWISEnnb_60arcsec.vot', format = 'votable')
    lk = np.lexsort([t['IRAS_CNTR'], t['angDist']])
    _, lu = np.unique(t[lk]['IRAS_CNTR'], return_index = True)
    t[lk[lu]].write('Ab2015_WISEAllSky_AllWISEnnb_60arcsec.vot', format = 'votable', overwrite = True)

def getAllWISE_nnb():
    #Combine AllWISE matches from the three methods and grab only the nearest neighbour from each method
    #   for a given IRAS source.
    f = ['WISEAllSky_AllWISEnnb', 'IRAS_coord_x_AllWISE', 'IRAS_AKARI_bestcoord_x_AllWISE']
    files = ['Ab2015_' + x.strip() + '_60arcsec.vot' for x in f]
    l = []
    for i in range(3):
        t = Table.read(files[i], format = 'votable')
        lk = np.lexsort([t['angDist'], t['IRAS_CNTR']])
        _, uk = np.unique(t[lk]['IRAS_CNTR'], return_index = True)
        l.append(t[lk[uk]])
    #aa = Table.read('Ab2015_WISEAllSky_AllWISEnnb_60arcsec.vot', format = 'votable')
    ## 
    #b = Table.read('Ab2015_IRAS_coord_x_AllWISE_60arcsec.vot', format = 'votable')
    ##Get nearest neighbour only
    #lk = np.lexsort([b['angDist'], b['IRAS_CNTR']]); _, uk = np.unique(b[lk]['IRAS_CNTR'], return_index = True)
    #bb = b[lk[uk]].copy(); b = 0
    ## 
    #c = Table.read('Ab2015_IRAS_AKARI_bestcoord_x_AllWISE_60arcsec.vot', format = 'votable')
    ##Get nearest neighbour only
    #lk = np.lexsort([c['angDist'], c['IRAS_CNTR']]); _, uk = np.unique(c[lk]['IRAS_CNTR'], return_index = True)
    #cc = c[lk[uk]].copy(); c = 0
    ## 
    i = Table.read('Ab2015.xml', format = 'votable')
    j = Table([i['IRAS_CNTR']]); i = 0
    # 
    p = join(join(join(j, a[0], keys = 'IRAS_CNTR', join_type = 'left'), a[1], keys = 'IRAS_CNTR', join_type = 'left'), \
             a[2], keys = 'IRAS_CNTR', join_type = 'left') 
    p.write('Ab2015_AllWISE_combined.xml', format = 'votable', overwrite = True)
    #len(np.nonzero((p['ID_1'] == p['ID_2']) & (p['ID_1'] == p['ID']))[0])
    #
    s1 = Table([p['IRAS_CNTR'], p['ID_1'], p['AllWISE_1']], names = ('IRAS_CNTR', 'AllWISE_ID', 'AllWISE_Designation'))
    k1 = np.nonzero(~s1['AllWISE_ID'].mask)[0]
    s1[k1].write('Ab2015_WISEAllSky_AllWISE_nnb.xml', format = 'votable', overwrite = True)
    s2 = Table([p['IRAS_CNTR'], p['ID_2'], p['AllWISE_2']], names = ('IRAS_CNTR', 'AllWISE_ID', 'AllWISE_Designation'))
    k2 = np.nonzero(~s2['AllWISE_ID'].mask)[0]
    s2[k2].write('Ab2015_IRAS_AllWISE_nnb.xml', format = 'votable', overwrite = True)
    s3 = Table([p['IRAS_CNTR'], p['ID'], p['AllWISE']], names = ('IRAS_CNTR', 'AllWISE_ID', 'AllWISE_Designation'))
    k3 = np.nonzero(~s3['AllWISE_ID'].mask)[0]
    s3[k3].write('Ab2015_IRAS_AKARI_AllWISE_nnb.xml', format = 'votable', overwrite = True)

def getAllWISE_unique():
    """Have to make a decision. Do we retain all the WISE matches and proceed to find matches to all of them in various
    catalogs, or do we just save the nearest neighbours?"""
    #t = Table.read('Ab2015_AllWISE_GaiaDR2_combined.xml', format = 'votable')
    #a = np.concatenate([t[~t['GaiaDR2_ID_1'].mask]['GaiaDR2_ID_1'].data.data, \
    #                    t[~t['GaiaDR2_ID_2'].mask]['GaiaDR2_ID_2'].data.data, \
    #                    t[~t['GaiaDR2_ID'].mask]['GaiaDR2_ID'].data.data])
    #_, u = np.unique(a, return_index = True)
    #q = Table([Column(a[u], name = 'GaiaDR2_ID')])
    #q.write('Ab2015_unique_GaiaDR2_IDs.xml', format = 'votable')


def getAllWISE():
    """Three methods to find AllWISE matches within 60":
    (1) Use the WISE AllSky coordinate.
    (2) Use the IRAS best coordinate.
    (3) Use the IRAS/AKARI best coordinate.
    """
    #Use the WISE AllSky coordinate
    getWISEAllSkyAllWISEnnb(radius = 60)
    #Use the IRAS best coordinate
    for i in range(3):
        cdsxmatch(cat1 = open('Ab2015_IRAS_coord_part' + str(i+1) + '.xml', 'r'), cat2 = 'vizier:II/328/allwise',
                  match_radius = 60, colRA1 = 'RAir', colDEC1 = 'DEir',
                  outfile = 'Ab2015_IRAS_coord_part' + str(i+1) + '_x_AllWISE_60arcsec', outformat = 'votable', \
                  overwrite = True)
    combineVOTables(infiles = ['Ab2015_IRAS_coord_part' + str(i+1) + '_x_AllWISE_60arcsec.vot' for i in range(3)],
                    outfile = 'Ab2015_IRAS_coord_x_AllWISE_60arcsec.vot')
    #Use the IRAS/AKARI best coordinate
    for i in range(3):
        cdsxmatch(cat1 = open('Ab2015_IRAS_AKARI_bestcoord_part' + str(i+1) + '.xml', 'r'), cat2 = 'vizier:II/328/allwise',
                  match_radius = 60, colRA1 = 'BEST_RA', colDEC1 = 'BEST_DEC', 
                  outfile = 'Ab2015_IRAS_AKARI_bestcoord_part' + str(i+1) + '_x_AllWISE_60arcsec', outformat = 'votable', \
                  overwrite = True)
    combineVOTables(infiles = ['Ab2015_IRAS_AKARI_bestcoord_part' + str(i+1) + '_x_AllWISE_60arcsec.vot' for i in range(3)],
                    outfile = 'Ab2015_IRAS_AKARI_bestcoord_x_AllWISE_60arcsec.vot')
    #
    #Combine AllWISE IDs from the three tables, save unique IDs.
    getAllWISE_unique()
    #Combine results of the three tables and only pick the nearest neighbour from each method.
    getAllWISE_nnb()

def getAllWISEMSX():
    """ Simpler alternative to the above function based on improved positional accuracy of MSX compared to IRAS/AKARI

Three methods to find AllWISE matches within 60":
    (1) Use the WISE AllSky coordinate.
    (2) Use the IRAS best coordinate.
    (3) Use the IRAS/AKARI best coordinate.
    """
    #Use the WISE AllSky coordinate
    getWISEAllSkyAllWISEMSXnnb(radius = 3)
    

def getGaiaDR2():
    queryurl = 'http://gea.esac.esa.int/tap-server/tap'
    un = ['WISEAllSky_AllWISE', 'IRAS_AllWISE', 'IRAS_AKARI_AllWISE']
    ur = ['Ab2015_' + item.strip() + '_nnb.xml' for item in un]
    out = ['Ab2015_' + item.strip() + '_nnb_GaiaDR2.xml' for item in un]
    for i in range(3):
        query = """select IRAS_CNTR, AllWISE_ID, AllWISE_Designation, source_id as GaiaDR2_ID, 
        best_neighbour_multiplicity as GaiaDR2_mult, angular_distance as GaiaDR2_AllWISE_dist
        from tap_upload.""" + un[i].strip()+ """ a left join gaiadr2.allwise_best_neighbour b
        on a.AllWISE_Designation = b.original_ext_source_id
        """
        q = Table([[queryurl], [query], [ur[i]], [un[i]], [out[i]]], \
                  names = ('queryurl', 'query', 'upload_resource', 'upload_table_name', 'outfile'))
        r = doquery_from_table(q)

    #Rename columns, combine into one entry per IRAS source.
    a = []
    for i in range(3):
        t = Table.read(out[i], format = 'votable')
        try:
            t.rename_columns(['iras_cntr', 'allwise_id', 'allwise_designation', 
                              'gaiadr2_id', 'gaiadr2_mult', 'gaiadr2_allwise_dist'], 
                             ['IRAS_CNTR', 'AllWISE_ID', 'AllWISE_Designation', 
                              'GaiaDR2_ID', 'GaiaDR2_mult', 'GaiaDR2_AllWISE_dist'])
            t.write(out[i], format = 'votable', overwrite = True)
        except:
            print("getGaiaDR2: Columns have already been renamed, doing nothing.")
        lk = np.lexsort([t['GaiaDR2_AllWISE_dist'], t['IRAS_CNTR']])
        _, uk = np.unique(t[lk]['IRAS_CNTR'], return_index = True)
        a.append(t[lk[uk]])
    i = Table.read('Ab2015.xml', format = 'votable'); j = Table([i['IRAS_CNTR']]); i = 0
    p = join(join(join(j, a[0], keys = 'IRAS_CNTR', join_type = 'left'),
                  a[1], keys = 'IRAS_CNTR', join_type = 'left'), 
             a[2], keys = 'IRAS_CNTR', join_type = 'left')
    p.write('Ab2015_AllWISE_GaiaDR2_combined.xml', format = 'votable', overwrite = True)
    #extract array containing unique instances of GaiaDR2_ID after combining the three columns from three
    #   matching methods
    t = Table.read('Ab2015_AllWISE_GaiaDR2_combined.xml', format = 'votable')
    a = np.concatenate([t[~t['GaiaDR2_ID_1'].mask]['GaiaDR2_ID_1'].data.data, \
                        t[~t['GaiaDR2_ID_2'].mask]['GaiaDR2_ID_2'].data.data, \
                        t[~t['GaiaDR2_ID'].mask]['GaiaDR2_ID'].data.data])
    _, u = np.unique(a, return_index = True)
    q = Table([Column(a[u], name = 'GaiaDR2_ID')])
    q.write('Ab2015_unique_GaiaDR2_IDs.xml', format = 'votable')

def getGaiaDR2_othercats():
    #Use the GaiaDR2 IDs to get neighbours in other catalogues.
    queryurl = 'http://gea.esac.esa.int/tap-server/tap'
    upload_table_name = 'unique_GaiaDR2_IDs'
    upload_resource = 'Ab2015_unique_GaiaDR2_IDs.xml'
    #Matches in the PANSTARRS catalogue
    outfile = 'Ab2015_unique_GaiaDR2_IDs_x_PANSTARRS1.xml'
    query = """select a.source_id as GaiaDR2_ID, a.original_ext_source_id as GaiaDR2_PANSTARRS1_ID, 
    a.angular_distance as GaiaDR2_PANSTARRS1_dist
    from tap_upload.""" + upload_table_name.strip() + """ c left join gaiadr2.panstarrs1_best_neighbour a
    on c.GaiaDR2_ID = a.source_id"""
    q = Table([[queryurl], [query], [upload_resource], [upload_table_name], [outfile]], \
              names = ('queryurl', 'query', 'upload_resource', 'upload_table_name', 'outfile'))
    r = doquery_from_table(q)

    #Matches in the TMASS catalogue
    outfile = 'Ab2015_unique_GaiaDR2_IDs_x_TMASS.xml'
    query = """select a.source_id as GaiaDR2_ID, a.original_ext_source_id as GaiaDR2_TMASS_ID, 
    a.angular_distance as GaiaDR2_TMASS_dist
    from tap_upload.""" + upload_table_name.strip() + """ c left join gaiadr2.tmass_best_neighbour a
    on c.GaiaDR2_ID = a.source_id"""
    q = Table([[queryurl], [query], [upload_resource], [upload_table_name], [outfile]], \
              names = ('queryurl', 'query', 'upload_resource', 'upload_table_name', 'outfile'))
    r = doquery_from_table(q)


def doMSX():
    """ 
    This function is intended to replicate the current functionality of 
    doall but modified to meet our needs to start from MSX instead of IRAS. 
    For MSX, we can trust the coordinates relatively well, unlike IRAS, thanks
    to the improved angular resolution. Hence, instead of trying to use
    AKARI to improve the coordinate precision, we just have to match MSX to 
    AKARI to get the fluxes, then propagate that through to (ALL)WISE/2MASS/Gaia
    just the same as with the IRAS Bestpos, but using the MSX coordinates as 
    the Bestpos.
    """

    #Open the table of MSX sources with no IRAS counterparts
    t = Table.read("MSXonly.vot", format="votable")

    #Now, we don't need to update the bestcoord, but we do need to get AKARI information


    #Goddamnit Sundar, I have to re-write your above functions for the MSX matches because the filenames are hardcoded.
    getAllWISEMSX()
    pass

def doall():
    #Run query to download the Abrahamyan+ 2015 catalogue.
    getAb2015()
    #Generate table of best coordinates from IRAS/AKARI to get AllWISE matches.
    bestcoord_IRAS_AKARI()
    #Get AllWISE matches with radius 60".
    getAllWISE()
    #Get Gaia DR2 matches to these three sets of AllWISE matches
    getGaiaDR2()
    #
    #
    #
    #extract matching info for the NESS subset.
    t = Table.read('Ab2015_AllWISE_GaiaDR2_combined.xml', format = 'votable')
    i = Table.read('Ab2015.xml', format = 'votable')
    n = Table.read('../NESScatalog.xml', format = 'votable')
    j = Table([i['IRAS-PSC'], i['IRAS_CNTR']], names = ('IRASPSC', 'IRAS_CNTR'))
    nn = join(join(n, j, keys = 'IRASPSC', join_type = 'left'), t, keys = 'IRAS_CNTR', join_type = 'left')
    #
    #Query to select the highest proper motion objects from the Gaia catalog that also have an AllWISE neighbour.
    #   the proper motion lower limit is set to 350 mas/yr.
    #select a.source_id as GaiaDR2_ID, ra as GaiaDR2_RA, dec as GaiaDR2_DEC, 
    #original_ext_source_id as AllWISE_Designation, angular_distance as GaiaDR2_AllWISE_dist
    #from gaiadr2.gaia_source a left join gaiadr2.allwise_best_neighbour b
    #on a.source_id = b.source_id
    #where abs(pmra) > 350 or abs(pmdec) > 350




    #For this table, use the WISE_MATCH_RADIUS to trim anything outside it. 

    #Get Gaia DR2 matches to those AllWISE matches.
    ##Split the output CSV fro above into five pieces, to limit the size of the resulting xmatched files.
    t1 = Table.from_pandas(pd.read_csv('/Users/sundar/Desktop/temp2/Ab2015_IRAS_AKARI_bestcoord_part1_x_AllWISE_60arcsec.csv', delimiter = ','))
    t2 = Table.from_pandas(pd.read_csv('/Users/sundar/Desktop/temp2/Ab2015_IRAS_AKARI_bestcoord_part2_x_AllWISE_60arcsec.csv', delimiter = ','))
    x = vstack([t1, t2])
    nt = len(x)
    ##Preserve only relevant columns
    s = Table([t['IRAS_CNTR'], t['ID'], t['RAJ2000'], t['DEJ2000']], names = ('IRAS_CNTR', 'ALLWISE_CNTR', 'ALLWISE_RA', 'ALLWISE_DEC'))
    s[:int(np.round(nt/5))].write('Ab2015_IRAS_AKARI_bestcoord_x_AllWISE_60arcsec_part1.xml', format = 'votable', overwrite = True)
    s[int(np.round(nt/5)):2*int(np.round(nt/5))].write('Ab2015_IRAS_AKARI_bestcoord_x_AllWISE_60arcsec_part2.xml', format = 'votable', overwrite = True)
    s[2*int(np.round(nt/5)):3*int(np.round(nt/5))].write('Ab2015_IRAS_AKARI_bestcoord_x_AllWISE_60arcsec_part3.xml', format = 'votable', overwrite = True)
    s[3*int(np.round(nt/5)):4*int(np.round(nt/5))].write('Ab2015_IRAS_AKARI_bestcoord_x_AllWISE_60arcsec_part4.xml', format = 'votable', overwrite = True)
    s[4*int(np.round(nt/5)):].write('Ab2015_IRAS_AKARI_bestcoord_x_AllWISE_60arcsec_part5.xml', format = 'votable', overwrite = True)

    #For subset of Abrahamyan+ 2015 table that have WISE All-Sky, get Gaia DR2 matches to those positions.
    #For this table, use the WISE_MATCH_RADIUS to trim anything outside it. 
    #I/345/gaia2

    #########
    #get nearest neighbour
    #lk = np.lexsort([t['angDist'], t['IRAS_CNTR']])
    #_, lu = np.unique(t[lk]['IRAS_CNTR'])
    #isolate rows where AllWISE IDs agree
    #k = np.nonzero(t['ID_1'] == t['ID_2'])

