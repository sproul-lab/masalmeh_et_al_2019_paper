# Usage
# gzip -dc file.fq.gz | demultiplex_NuGen_RRBS.py [index_table file_suffix]

#!/usr/bin/env python
import sys
import re
# import gzip

# use stdin to allow piping from gzip
if not sys.stdin.isatty():
    input_stream = sys.stdin

# otherwise, give an error
else:
    message = 'Must have input from stdin!'
    raise IndexError(message)

# Initialise the barcodes from file
# Use hash table structure
# Allows checking if the barcode exists and then automatically getting sample ID

# open file and initiate index dictionary
indexFile = sys.argv[1]
f = open(indexFile, 'r')
indexDict= {}

# Capture the file suffix
file_suffix = sys.argv[2]

# read in the index file and parse into the dictionary
# store file names and open for appending
# call using index as key later
for line in f:
    line=line.rstrip("\n")
    index, sample = line.split('\t')
    outFileName = '%s_%s.fq' % (sample, file_suffix)
    indexDict[index] = open(outFileName,'a')

# Close the index file
f.close

# Debug test as hardwired
# indexDict= {}
# indexDict['AGGAGG'] = 'Sample_1'
# indexDict['AAGCCT'] = 'Sample_2'

# Initialise an unrecognised index flag
IndexFlag = False

# Variable to store current index
currIndex=""

# Initialise the counts of the number of reads demultiplexed
readCount = 0
indexCount = 0
noIndexCount = 0

# Initialise a rolling line counter
lineCounter=4

# Start processing the lines of the fastq file
for line in input_stream:

    # Check if this line has a read header (when lineCounter = 4)
    if lineCounter==4:
        # If so, check it is an index we want and capture the index
        # Capture the index
        currIndex=line[-13:-7]

        # Check if it is in the index dictionary and set no index flag accordingly
        # Also increment the read count and reset the lineCounter
        if currIndex in indexDict:
            IndexFlag = True
            readCount += 1
            indexCount += 1
            lineCounter = 0
            # Print a message for every 1000000 reads processed
            if readCount % 1000000 == 0:
                print 'Processed ', readCount, ' reads'

        else:
            IndexFlag = False
            readCount += 1
            noIndexCount += 1
            lineCounter = 0
            # Print a message for every 1000000 reads processed
            if readCount % 1000000 == 0:
                print 'Processed ', readCount, ' reads'

    # Now print the current line to the appropriate file if noIndexFlag is FALSE
    if IndexFlag:
        # If it is an index we want print to appropriate file
        # Call using the dictionary (all should be open)
        indexDict[currIndex].write(line)

        # Debug line
        # print currIndex, indexDict[currIndex], IndexFlag, ': ', line

        # How to do a compressed write
        # Doesn't work easily with dictionary of file handles
        # with gzip.open('%s_%s.fq.gz' % (indexDict[currIndex], file_suffix), 'ab') as fo:
            # fo.write(line)
        # fo.close

    # Now increment the line counter by 1
    lineCounter += 1

# Print a message for finishing the parsing
print '*************************'
print 'Finished!'
print 'Reads parsed: ', readCount
print 'Reads with indexes: ', indexCount
print 'Reads without indexes: ', noIndexCount

# Now close all the file handles
for index in indexDict:
    indexDict[index].close
