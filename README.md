# matching
Finding counterparts to the IRAS PSC/FSC in various large surveys.

The Abrahamyan et al. (2015) catalogue, which has matched the PSC and FSC to the AKARI IRC/FIS and WISE All-Sky surveys. Starting with their table, we extend the matching to large surveys such as AllWISE, 2MASS, Gaia DR2, and many others.

The scripts in this repository are meant to reproduce the entire dataset used for our analysis.

To start, the file doall_Ab2015.py [requirements: numpy, astropy, astroquery] stores the queries in an astropy table, then uses astroquery to run TAP queries from this table. The output from these queries is again in the form of astropy tables, which are combined to generate the complete catalogue.

Once this full catalogue is generated, the rows corresponding to the NESS sample can be extracted (the file NESScatalog.xml contains the IRAS PSC IDs of the NESS targets without the "IRAS " prefix).
