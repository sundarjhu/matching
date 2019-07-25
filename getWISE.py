from astropy.table import Table
import pandas as pd
from bestcoord import * #this imports write_IPAC_VOTable
import os, subprocess, time

def getWISE():
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
                        '-F', '"QUERY=SELECT TAP_UPLOAD.my_table.IRAS_PSC_ID, TAP_UPLOAD.my_table.IRAS_FSC_ID, ' + \
                        'allwise_p3as_psd.tmass_key FROM allwise_p3as_psd, TAP_UPLOAD.my_table ' + \
                        'WHERE TAP_UPLOAD.my_table.ID = allwise_p3as_psd.cntr"', \
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
                    subprocess.run(['curl', '-o ' + outfiles[i], url[i] + '/results/result'])
                    break
                elif b'ERROR' in op or b'ABORTED' in op:
                    print("Uh oh, query status: " + op)
                    break
                else:
                    #check again in 10 min.
                    time.sleep(600)

        return 0

    #bestcoord_IRAS_AKARI()
    infiles = ['Ab2015_IRAS_AKARI_combined_part1_x_AllWISE_60arcsec.csv', \
               'Ab2015_IRAS_AKARI_combined_part2_x_AllWISE_60arcsec.csv']
    outfiles = ['Ab2015_IRAS_AKARI_combined_part1_x_AllWISE_60arcsec_TMASS_KEY.xml', \
                'Ab2015_IRAS_AKARI_combined_part2_x_AllWISE_60arcsec_TMASS_KEY.xml']
    print(bcolors.WARNING + """You will have to perform a cross-match between the Abrahamyan et al. 2015 table and
    the AllWISE catalog. Instructions follow:
    1) Go to http://cdsxmatch.u-strasbg.fr/
    2) Upload tables Ab2015_IRAS_AKARI_combined_bestcoord_part*.csv
    3) Type "AllWISE" in the right-hand tab.
    4) Expand "Show options" and set a radius of 60".
    5) Click "Begin the X-Match". This takes about 20 min and generates a >300 MB file.
    6) Download the output in CSV format (because we have issues with their VOT format) 
          into Ab2015_IRAS_AKARI_combined_part*_x_AllWISE_60arcsec.csv""" + bcolors.ENDC)
    input("Press ENTER when done: ")
    for i in range(len(infiles)):
        t = Table.from_pandas(pd.read_csv(infiles[i], delimiter = ','))
        s = Table([t['IRAS_PSC_ID'], t['IRAS_FSC_ID'], t['ID']])
        infiles[i] = infiles[i].replace('csv', 'vo')
        write_IPAC_VOTable(s, outfile = infiles[i])

    infiles = [infile.replace('csv', 'vo') for infile in infiles]
    #for i in range(len(infiles)):
    #    subprocess.run(['cp', infiles[i], 'in' + str(i+1) + '.vo'])
    #    infiles[i] = 'in' + str(i+1) + '.vo'

    r = run_IPAC_query(infiles, outfiles)
