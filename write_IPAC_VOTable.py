import os
import numpy as np

def write_IPAC_VOTable(t, outfile = 'out.vo'):
    #Given an astropy table t, write a VOTable in the IPAC preferred format.
    #This is a terrible kludge because it writes to a temporary file just to manipulate
    #   the VOTable header. Find a better way.
    ncols = len(t.columns)
    hlines = []
    for i in range(ncols):
        line = '<FIELD name="' + t.columns[i].name + '"'
        if ('byte' in t.columns[i].dtype.name) or ('str' in t.columns[i].dtype.name):
            line = line + ' datatype="char" arraysize="*"'
        elif 'int' in t.columns[i].dtype.name:
            line = line + ' datatype="long"'
        elif 'float' in t.columns[i].dtype.name:
            line = line + ' datatype="double"'
        else:
            line = line + ' datatype="' + t.columns[i].dtype.name + '"'
        try:
            line = line + ' unit="' + t.columns[i].unit.name + '"'
        except:
            pass
        line = line + '/>'
        hlines.append(line)
    header = """<?xml version='1.0'?>
    <VOTABLE version="1.3"
    xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
    <RESOURCE>
    <TABLE>"""
    header = header.splitlines()
    header.append('<DESCRIPTION>' + outfile + ' in VO format</DESCRIPTION>')
    header = header + hlines

    #write to a temporary table so it can be read in line by line. Once read, delete the temporary file.
    t.write('temp', format = 'votable', overwrite = True)
    with open('temp') as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    i = content.index('<DATA>')
    outlines = header + content[i:]
    with open(outfile, 'w') as f:
        for line in outlines:
            f.write("%s\n" % line)
    os.remove('temp')

