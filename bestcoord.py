import numpy as np
from astropy.table import Table
from write_IPAC_VOTable import *

def bestcoord_IRAS_AKARI():
    """Given the IRAS/AKARI combined VOTable, store the best coordinate pair
    obtained hierarchically: AKARI IRC > AKARI FIS > IRAS PSC > IRAS FSC.
    The WISE and 2MASS data are currently ignored, as they are obtained in
    subsequent queries."""
    t = Table.read('Ab2015_IRAS_AKARI_combined.xml', format = 'votable')
    if 'RA_BEST' in list(t.columns):
        print("Best coordinate fields already exist in table Ab2015_IRAS_AKARI_combined.xml, modifying...")
    ndata = len(t)
    ra_best = t['AKARI_IRC_RA'].copy(); dec_best = t['AKARI_IRC_DEC'].copy(); best_coord = np.repeat('AKARI_IRC', ndata)
    ra_best[list(t['AKARI_IRC_RA'].mask)] = t['AKARI_FIS_RA'][list(t['AKARI_IRC_RA'].mask)].copy()
    dec_best[list(t['AKARI_IRC_DEC'].mask)] = t['AKARI_FIS_DEC'][list(t['AKARI_IRC_DEC'].mask)].copy()
    best_coord[list(t['AKARI_IRC_RA'].mask)] = 'AKARI_FIS'
    ra_best[list(t['AKARI_FIS_RA'].mask)] = t['IRAS_PSC_RA'][list(t['AKARI_FIS_RA'].mask)].copy()
    dec_best[list(t['AKARI_FIS_DEC'].mask)] = t['IRAS_PSC_DEC'][list(t['AKARI_FIS_DEC'].mask)].copy()
    best_coord[list(t['AKARI_IRC_RA'].mask)] = 'IRAS_PSC'
    ra_best[list(t['IRAS_PSC_RA'].mask)] = t['IRAS_FSC_RA'][list(t['IRAS_PSC_RA'].mask)].copy()
    dec_best[list(t['IRAS_PSC_DEC'].mask)] = t['IRAS_FSC_DEC'][list(t['IRAS_PSC_DEC'].mask)].copy()
    best_coord[list(t['IRAS_PSC_RA'].mask)] = 'IRAS_FSC'
    t['BEST_RA'] = ra_best; t['BEST_DEC'] = dec_best; t['BEST_COORD_SOURCE'] = best_coord
    t.write('Ab2015_IRAS_AKARI_combined.xml', format = 'votable', overwrite = True)
    #Also write out a table with only the relevant columns for coordinate matching with AllWISE
    s = Table([t['IRAS_PSC_ID'], t['IRAS_FSC_ID'], t['BEST_RA'], t['BEST_DEC']])
    s['IRAS_PSC_ID'][s['IRAS_PSC_ID'] == ''] = '--'
    s['IRAS_FSC_ID'][s['IRAS_FSC_ID'] == ''] = '--'
    s[:round(len(s)/2)].write('Ab2015_IRAS_AKARI_combined_bestcoord_part1.xml', format = 'votable', overwrite = True)
    s[round(len(s)/2):].write('Ab2015_IRAS_AKARI_combined_bestcoord_part2.xml', format = 'votable', overwrite = True)
    #write_IPAC_VOTable(s, outfile = 'Ab2015_IRAS_AKARI_combined_bestcoord.vo')

def bestcoord():
    """Given the IRASxAKARIxWISEx2MASS match VOTable, store the best coordinate pair
    obtained hierarchically.
    Highest precedence is given to 2MASS position. If this is not available, 
    the IRAS_best_pos_* fields, which contain the best position of 
    (IRAS FSC/PSC, AKARI, WISE) from Abrahamyan et al. 2015, are chosen.
    NOTE: this script overwrites the original xml and csv tables, but the fits
    version remains unchanged."""
    t = Table.read('Ab2015_IRAS_AKARI.xml', format = 'votable')
    ra2m = t['tmass_ra']
    dec2m = t['tmass_dec']
    t['IRAS_AKARI_WISE_2MASS_BEST_RA'] = t['IRAS_best_pos_RA']
    t['IRAS_AKARI_WISE_2MASS_BEST_DEC'] = t['IRAS_best_pos_DEC']
    k = np.where(np.isfinite(t['tmass_ra']))[0]
    if len(k) != 0:
        t[k]['IRAS_AKARI_WISE_2MASS_BEST_RA'] = t[k]['tmass_ra']
        t[k]['IRAS_AKARI_WISE_2MASS_BEST_DEC'] = t[k]['tmass_dec']
    t.write('IRAS_AKARI_WISE_2MASS.xml', overwrite = True, format = 'votable')
    t.write('IRAS_AKARI_WISE_2MASS.csv', overwrite = True, format = 'csv') #deprecated


