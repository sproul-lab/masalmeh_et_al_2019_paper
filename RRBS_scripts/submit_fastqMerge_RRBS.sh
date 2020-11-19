#!/bin/sh

# Options for runnig on a cluster need to be specified
# Have removed the options required for University of Edinburgh and MRC HGU compute cluster

# Usage
# qsub submit_fastqMerge_RRBS.sh -d fileDIR -i sampleID -o outDIR

# Assumes that all files are present in one directory and have two FASTQ files per lane (_1 and_2)
# .fq.gz extension assumed also
# Searches for L00N in name also as this is the format our files were supplied in

# Set an exit method for printing an error
die() {
    printf '%s\n' "$1" >&2;
    exit 1;
}

# Capture command line arguments
while :; do
    case $1 in
        -i)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                sampleID=$2
                shift
            else
                die 'ERROR: "-i" requires a non-empty option argument.'
            fi
            ;;
        -d)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                fileDIR=$2
                shift
            else
                die 'ERROR: "-d" requires a non-empty option argument.'
            fi
            ;;
        -o)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                outDIR=$2
                shift
            else
                die 'ERROR: "-o" requires a non-empty option argument.'
            fi
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)               # Default case: No more options, so break out of the loop.
            break
    esac
    shift
done

# Check required command line options
if [ -z ${sampleID+x} ]; then die 'ERROR: -i is unset, this script needs an input sample name.'; fi
if [ -z ${fileDIR+x} ]; then die 'ERROR: -d is unset, this script needs an input data directory.'; fi
if [ -z ${outDIR+x} ]; then die 'ERROR: -o is unset, this script needs an output directory.'; fi

# Set the working directory
rootDIR=$(pwd)

# Check to see if the output directory exists and if not create it
if [[ ! -d "${outDIR}" ]] ; then
    echo "${outDIR} does not exist, creating...."
    mkdir ${outDIR}
fi

# Add headers to STDOUT and STDERR
echo STDOUT:
echo -e JOB_ID = $JOB_ID"\n"
hostname
echo -e Merge WTCRF RRBS FASTQ from different lanes for a given sample.
echo -e Data Directory: ${fileDIR}
echo -e Sample: ${sampleID}
echo -e Output DIR: ${outDIR}
echo -e Root DIR: ${rootDIR}

echo Start Time:
date

echo STDERR:  >&2
echo -e JOB_ID = $JOB_ID"\n"  >&2
hostname  >&2
echo -e Merge WTCRF RRBS FASTQ from different lanes for a given sample. >&2
echo -e Data Directory: ${fileDIR} >&2
echo -e Sample: ${sampleID} >&2
echo -e Output DIR: ${outDIR} >&2
echo -e Root DIR: ${rootDIR} >&2

# No modules required

################################
# Copy files to TMPDIR on the Node
################################

# Generate a list of files to work on
# Use the R1 files to get them
ls -1 ${fileDIR} | grep -e "^${sampleID}_L00[0-9]_1.fq.gz$" | sed 's/_1.fq.gz//' > $TMPDIR/filePrefixList.txt

echo '************************************'
echo List of file prefixes:
echo '************************************'
cat $TMPDIR/filePrefixList.txt
echo

# Now copy the files over
# Do R1 then R2
cat $TMPDIR/filePrefixList.txt | while IFS= read -r LINE; do cp ${fileDIR}/${LINE}_1.fq.gz $TMPDIR/${LINE}_1.fq.gz; done
cat $TMPDIR/filePrefixList.txt | while IFS= read -r LINE; do cp ${fileDIR}/${LINE}_2.fq.gz $TMPDIR/${LINE}_2.fq.gz; done

# Go to TMPDIR
cd $TMPDIR

echo '************************************'
echo Contents of TMPDIR after copying:
echo '************************************'
ls -1
echo

################################
# Concatenate the files
################################


# Do Read 1
echo "Concatenating R1 files"
cat $(sed 's/$/_1.fq.gz/' filePrefixList.txt) > ${sampleID}_1.fq.gz
# Do Read 2
echo "Concatenating R2 files"
cat $(sed 's/$/_2.fq.gz/' filePrefixList.txt) > ${sampleID}_2.fq.gz

################################
# Copy files back over to Eddie filesystem
################################


echo "Copying merged R1 file back to output directory"
cp ${sampleID}_1.fq.gz ${rootDIR}/${outDIR}/${sampleID}_1.fq.gz
echo "Copying merged R2 file back to output directory"
cp ${sampleID}_2.fq.gz ${rootDIR}/${outDIR}/${sampleID}_2.fq.gz

echo '************************************'
echo Contents of TMPDIR at end of script
echo '************************************'
ls -1
echo

echo '************************************'
echo Finished:
echo '************************************'
date

####################################
# Tidy up, copy StdErr and StdOut
####################################

# Return to root
cd ${rootDIR}

# Move STDOUT and STDERR to files based on SRR ID for easier diagnosis of problems
mv ${JOB_ID}_out ${outDIR}/${sampleID}_${JOB_ID}_fastqMerge_RRBS_stdout
mv ${JOB_ID}_err ${outDIR}/${sampleID}_${JOB_ID}_fastqMerge_RRBS_stderr
