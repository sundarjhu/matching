import numpy as np
import requests as r
from astropy.table import Table
from astropy.io.votable import parse, parse_single_table #, VOTableFile
#from astroquery.xmatch import XMatch
from astroquery.vizier import Vizier

if __name__ == "__main__":
    repeat_query=False
    if repeat_query:

        res = r.post('http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync/tables',
                     data={'action':'getColList',
                           'tabname':'vizier:II/338/catalog',
                           'RESPONSEFORMAT': 'votable'}
        )
        print(res.text)
        res = r.post('http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync/tables',
                     data={'action':'getColList',
                           'tabname':'vizier:II/328/allwise',
                           'RESPONSEFORMAT': 'votable'}
        )
        print(res.text)
        #exit()
    
        res = r.post(
            'http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync',
            data={'request': 'xmatch',
                  'distMaxArcsec': 30,
                  'RESPONSEFORMAT': 'votable',
                  'cat1': 'vizier:II/338/catalog',
                  'cat2': 'vizier:II/328/allwise'#,
                  #'colRA1': 'RAJ2000',
                  #'colDec1': 'DEJ2000',
                  #'cols2':'2Mkey,d2M'
            },
            #files={'cat1': open('posList.vot', 'r')}
        )

        #print(res)
        ''' write out the return so we can test reading it back in without waiting forever for the query every time '''
        f = open("output.vot",'w')
        f.write(res.text)
        f.close()
    
    votable = parse('output.vot') #res.text)
    print(votable)
    #help(votable)
    tablist = []
    for resource in votable.resources:
        for table in resource.tables:
            tablist.append(table.to_table())
            tablist[-1].pprint()
            tablist[-1].write(str(table.ID)+'.vot',format="votable",overwrite=True)
    #tab.pprint()



    ''' 
        Not all columns of vizier tables are available to xMatch, 
        so we have to now query vizier by AllWISE ID to get the 2Mkey, 
        which gives us the link to 2MASS. 
    '''

    
