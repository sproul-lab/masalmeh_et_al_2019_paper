# Script to parse Bismark BEDgraph Coverage
# Converts to format that is readable by MethylRSeeker

# Load sys module
import sys

# open files needed
f = open(sys.argv[1], 'r')
fo = open(sys.argv[2], 'w')
strand = sys.argv[3]

# Print out file 
print "Processing file: "+sys.argv[1]+" Strand: "+strand

# cycle through file and print out
while True:
    line = f.readline()
    if not line: break
        
    line=line.rstrip("\n")
        
    # Split the lines
    line = line.split("\t")
    chr = line[0]
    pos = int(line[1])
    meth = int(line[4])
    unmeth = int(line[5])
    total = meth + unmeth
    
    # set the start and end of the CpG
    if (strand == "OT"):
        start = pos - 1
        end = pos + 1
    else:
        start = pos - 2
        end = pos

    # Now print out the data
    fo.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+str(total)+"\t"+str(meth)+"\n")

# Close files
f.close()
fo.close()

print "**** Done Processing"
print "Results file: " + sys.argv[2]