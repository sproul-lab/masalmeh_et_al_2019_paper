# Collates various QC output reports from Bismark alignment to make a table for NuGen RRBS only

#!/usr/bin/env python
import sys
import re

# use stdin if it's full
if not sys.stdin.isatty():
    input_stream = sys.stdin

# otherwise, read the given filename
else:
    try:
        input_stream = sys.argv[1]
    except IndexError:
        message = 'need filename as first argument if stdin is not full'
        raise IndexError(message)

# Read through the input line by line
for currPrefix in input_stream:
    currPrefix=currPrefix.rstrip("\n")
    
    # Reset the variables to store the data
    totalReads=0
    readsQC=0
    readsALIGN=0
    dupReads=0
    totalCG=0
    obsCG=0
    methCG=0
    CGcov=0
    totalCHG=0
    methCHG=0
    totalCHH=0
    methCHH=0
    totalLAMBDA=0
    methLAMBDA=0
    
    # Open the trimming report
    # Extract total number of reads
    f = open(currPrefix + "_1_trimming_report.txt", 'r')
    while True:
        line = f.readline()
        if not line: break
        line=line.rstrip("\n")
        m = re.search("^([0-9]+) sequences processed in total$", line)
        if m:
            totalReads=m.group(1)
    f.close

    # Open up the bismark alignment report
    # Get reads passing QC and alignment rate
    f = open(currPrefix + "_bismark_align_report.txt", 'r')
    while True:
        line = f.readline()
        if not line: break
        line=line.rstrip("\n")
        mQC = re.search("^Sequence pairs analysed in total:\t([0-9]+)$", line)
        mALIGN = re.search("^Number of paired-end alignments with a unique best hit:\t([0-9]+)$", line)
        if mQC:
            readsQC=mQC.group(1)
        elif mALIGN:
            readsALIGN=mALIGN.group(1)
    f.close
    
    # Open up the NuGen de-duplication report
    # Get no of duplicates
    # Note have parse a table of numbers to get correct column
    # Use position rather than regular expression to find it, 2 line table
    current_line = 1
    f = open(currPrefix + "_RRBS_read_duplication_report.txt", 'r')
    while True:
      line = f.readline()
      if not line: break
      if current_line == 2:
        dupReads = line.split("\t")[4]
      current_line += 1
    f.close
    
    # Now open up the autosomal CpG report
    # Get number of CpGs in total genome and seen in this alignment
    # Then get mean CpG depth and methylation level
    f = open(currPrefix + "_autosomal_CpG_report.txt", 'r')
    while True:
        line = f.readline()
        if not line: break
        line=line.rstrip("\n")
        mTOTAL = re.search("^Total CpGs processed:\t([0-9]+)$", line)
        mOBS = re.search("^Total CpGs covered:\t([0-9]+)$", line)
        mMETH = re.search("^Mean methylation:\t([0-9]+.[0-9]+)$", line)
        mCOV = re.search("^Mean coverage:\t([0-9]+.[0-9]+)$", line)
        if mTOTAL:
            totalCG=mTOTAL.group(1)
        elif mOBS:
            obsCG=mOBS.group(1)
        elif mMETH:
            methCG=mMETH.group(1)
        elif mCOV:
            CGcov=mCOV.group(1)
    f.close

    # Now open up methylation extractor report
    # Get number of CHG and CHH methylation levels
    f = open(currPrefix + "_bismark_extract_report.txt", 'r')
    while True:
        line = f.readline()
        if not line: break
        line=line.rstrip("\n")
        mTOTAL1 = re.search("^Total C to T conversions in CHG context:\t([0-9]+)$", line)
        mMETH1 = re.search("^Total methylated C's in CHG context:\t([0-9]+)$", line)
        mTOTAL2 = re.search("^Total C to T conversions in CHH context:\t([0-9]+)$", line)
        mMETH2 = re.search("^Total methylated C's in CHH context:\t([0-9]+)$", line)
        if mTOTAL1:
            totalCHG=mTOTAL1.group(1)
        elif mMETH1:
            methCHG=mMETH1.group(1)
        elif mTOTAL2:
            totalCHH=mTOTAL2.group(1)
        elif mMETH2:
            methCHH=mMETH2.group(1)
    f.close

    # Now open up the lambda report
    # Get total CG and meth CG
    f = open(currPrefix + "_lambda_CpG_report.txt", 'r')
    while True:
        line = f.readline()
        if not line: break
        line=line.rstrip("\n")
        mTOTAL = re.search("^Lambda T: ([0-9]+)$", line)
        mMETH = re.search("^Lambda M: ([0-9]+)$", line)
        if mTOTAL:
            totalLAMBDA=mTOTAL.group(1)
        elif mMETH:
            methLAMBDA=mMETH.group(1)
    f.close

    # Now print a report to standard out
    outLIST=[currPrefix, str(totalReads), str(readsQC), str(readsALIGN), str(dupReads), str(totalCG), str(obsCG), str(methCG), str(CGcov), str(totalCHG), str(methCHG), str(totalCHH), str(methCHH), str(totalLAMBDA), str(methLAMBDA)]
    print('\t'.join(map(str,outLIST)))
