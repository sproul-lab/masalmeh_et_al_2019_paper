#!/bin/sh

# Options for running on a cluster need to be specified
# Have removed the options required for University of Edinburgh and MRC HGU compute cluster

# Usage
# qsub submit_RRBS_demultiplex_NuGen.sh -d fileDIR -i fileID -m indexFILE -o outDIR [--scriptDIR ~]

# assumes you have paired end files

# Set root directory of scripts by default to ~
# This is usually passed by the executable to launch the qsub job
scriptDIR=~;

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
                fileID=$2
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
        -m)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                indexFILE=$2
                shift
            else
                die 'ERROR: "-m" requires a non-empty option argument.'
            fi
            ;;
        --scriptDIR)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                scriptDIR=$2
                shift
            else
                die 'ERROR: "--scriptDIR" requires a non-empty option argument.'
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
if [ -z ${fileID+x} ]; then die 'ERROR: -i is unset, this script needs an input file name.'; fi
if [ -z ${fileDIR+x} ]; then die 'ERROR: -d is unset, this script needs an input data directory.'; fi
if [ -z ${outDIR+x} ]; then die 'ERROR: -o is unset, this script needs an output directory.'; fi
if [ -z ${indexFILE+x} ]; then die 'ERROR: -m is unset, this script needs an index map file.'; fi

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
echo -e Demultiplex WTCRF NextSeq RRBS.
echo -e Data Directory: ${fileDIR}
echo -e Sample: ${fileID}
echo -e Index File: ${indexFILE}
echo -e Output DIR: ${outDIR}
echo -e Root DIR: ${rootDIR}

echo Start Time:
date

echo STDERR:  >&2
echo -e JOB_ID = $JOB_ID"\n"  >&2
hostname  >&2
echo -e Demultiplex WTCRF NextSeq RRBS. >&2
echo -e Data Directory: ${fileDIR} >&2
echo -e Sample: ${fileID} >&2
echo -e Index File: ${indexFILE} >&2
echo -e Output DIR: ${outDIR} >&2
echo -e Root DIR: ${rootDIR} >&2

# Configure modules
# These are specific versions used on University of Edinburgh MRC HGU cluster
# Will need to be adapted to run on other systems
module add python/2.7.10
module add R/3.3.0

# Copy FASTQ files to node
# Check to see if the file has extension 'fq.gz' or 'fastq.gz'
# Here we have hard wired PE
# Note: in both cases file is copied as 'fq.gz'
if [[ -f "${fileDIR}/${fileID}_1.fastq.gz" ]] ;
then
    echo "File is ${fileID}_1.fastq.gz"
    cp ${fileDIR}/${fileID}_1.fastq.gz $TMPDIR/${fileID}_1.fq.gz
    cp ${fileDIR}/${fileID}_2.fastq.gz $TMPDIR/${fileID}_2.fq.gz
elif [[ -f "${fileDIR}/${fileID}_1.fq.gz" ]] ;
then
    echo "Files is ${fileID}_1.fq.gz"
    cp ${fileDIR}/${fileID}_1.fq.gz $TMPDIR/${fileID}_1.fq.gz
    cp ${fileDIR}/${fileID}_2.fq.gz $TMPDIR/${fileID}_2.fq.gz
else
    die "ERROR: Can't find ${fileDIR}/${fileID} as fastq.gz or fq.gz";
fi

# Copy over the index file
cp ${indexFILE} $TMPDIR/indexFile.txt

# Go to TMPDIR and start analysis
cd $TMPDIR

echo '************************************'
echo Start index histogram:
echo '************************************'
date

# Note: this is only generated for the read 1 file as should be identical (always have been to date 15/3/19)

# First get the barcodes with AWK and analyse their frequency
gzip -dc ${fileID}_1.fq.gz  | awk '/^@/{a=length($0);print substr($0,a-11,6)}' | sort | uniq -c | sort -n > ${fileID}_IndexCounts.txt

# Get a histogram from the index counts
Rscript --vanilla ${scriptDIR}/eddie3_RRBSprocess/get_index_histogram_RRBS_NuGen.R indexFile.txt ${fileID}_IndexCounts.txt

# Output should be: {fileID}_IndexPlot.png

echo '************************************'
echo Start demultiplexing:
echo '************************************'
date

# Get the suffix to attach to the files
# Get rid of bits we don't need (eg _001) and change R1/2 to 1/2
outSuffix=$(echo ${fileID} | sed -e 's/^Undetermined_S0_//')

# Demultiplex the files
# R1
gzip -dc ${fileID}_1.fq.gz | python ${scriptDIR}/eddie3_RRBSprocess/demultiplex_RRBS_NuGen.py indexFile.txt ${outSuffix}_1
# R2
gzip -dc ${fileID}_2.fq.gz | python ${scriptDIR}/eddie3_RRBSprocess/demultiplex_RRBS_NuGen.py indexFile.txt ${outSuffix}_2

# Output should be set of files with: indexSampleID_${outSuffix}.fq
# indexSampleID comes from the ${indexFile} second column of sample names

# Compress the output FASTQ files
gzip *_${outSuffix}_1.fq
gzip *_${outSuffix}_2.fq

# Copy back the results
# Histogram plot of index frequency
cp ${fileID}_IndexPlot.png ${rootDIR}/${outDIR}/${fileID}_IndexPlot.png

# The output fastq files
# Do as loop from the indexfile
cut -f 2 indexFile.txt | while IFS= read -r LINE; do cp ${LINE}_${outSuffix}_1.fq.gz ${rootDIR}/${outDIR}/${LINE}_${outSuffix}_1.fq.gz; done
cut -f 2 indexFile.txt | while IFS= read -r LINE; do cp ${LINE}_${outSuffix}_2.fq.gz ${rootDIR}/${outDIR}/${LINE}_${outSuffix}_2.fq.gz; done

echo '************************************'
echo Error checking file list
echo '************************************'
ls -1
echo

echo '************************************'
echo Finished:
echo '************************************'
date


# Return to root
cd ${rootDIR}

# Move STDOUT and STDERR to files based on SRR ID for easier diagnosis of problems
mv ${JOB_ID}_out ${outDIR}/${fileID}_${JOB_ID}_RRBS_demultiplex_NuGen_stdout
mv ${JOB_ID}_err ${outDIR}/${fileID}_${JOB_ID}_RRBS_demultiplex_NuGen_stderr
