from __future__ import division

#!/usr/bin/env python
import sys
import re


# use stdin if it's full
if not sys.stdin.isatty():
    input_stream = sys.stdin

# otherwise, read the given filename
else:
    try:
        input_filename = sys.argv[1]
    except IndexError:
        message = 'need filename as first argument if stdin is not full'
        raise IndexError(message)
    else:
        input_stream = open(input_filename, 'rU')

# Initialise sum and number of lines
totsum = 0
totmethsum = 0
CpGobs = 0
CpGtot = 0

for line in input_stream:
    line=line.rstrip("\n")
    
    # Split the lines
    line = line.split("\t")
    chr = line[0]
    totf = int(line[3])
    methf = int(line[4])
    totr = int(line[5])
    methr = int(line[6])
    
    # check if the chromosome is a standard autosome and process if it is
    if (re.search("^chr[0-9]+$", chr)):
        # Now sum the totals for the two strands
        tot = totf + totr
        meth = methf + methr
        CpGtot += 1
        # If this CpG is observed increment counts
        if (tot > 0):
            CpGobs += 1
            totsum += tot
            totmethsum += meth
    # Debugging line
    # print (chr+"\t"+str(totf)+"\t"+str(totr)+"\t"+str(tot)+"\t"+str(meth)+"\t"+str(totsum)+"\t"+str(totmethsum)+"\t"+str(CpGobs)+"\t"+str(CpGtot))

# Now print a report to standard out
print ("Total CpGs processed:\t"+str(CpGtot))
print ("Total CpGs covered:\t"+str(CpGobs))
print ("Total coverage:\t"+str(totsum))
print ("Total methylated coverage:\t"+str(totmethsum))
print ("Mean methylation:\t"+str(totmethsum/totsum))
print ("Mean coverage:\t"+str(totsum/CpGobs))
