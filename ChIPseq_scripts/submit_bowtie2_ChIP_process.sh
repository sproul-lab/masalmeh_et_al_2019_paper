#!/bin/sh

# Options for running on a cluster need to be specified
# Have removed the options required for University of Edinburgh and MRC HGU compute cluster


# Usage
# qsub submit_bowtie2_ChIP_process.sh -d fileDIR -i sampleID -o outDIR -g indexDIR [--PE --noDeDup --noClean --indexROOT]

# Designed to work on SE phred33 files by default
# Will take a FASTQ file and align it against a given genome index
# Following alignment it will strip multimapping reads and PCR duplicates before sorting and indexing the output BAM file by default
# This cleanup can be turned off with flags if needed

# Note: no read clipping and phred 33 hardwired and can't be changed at present
# For PE is hard wired at maximum fragment size of 1000

# Arguments, these are all required
# -i sampleID, input file ID, note without '_1' or '_2' or fastq gz extension
# -d fileDIR, data directory, where FASTQ files are stored
# -o outDIR, directory to write output files to
# -g indexDIR, name of directory containing the index for alignment

# Optional command line parameters
# --PE, use PE FASTQ files
# --noDeDup, don't remove PCR duplicates
# --noClean, don't remove low quality mapping reads/fragments
# Will look inside the --indexROOT directory for genome indexes
# Was set as default for our cluster
# Will need to be modified for other systems
# So assumes index has same name as the indexDIR

# Note only produces a STDOUT (setting qsub option -j yes)

# Need to set the index root to use on other systems
# Can be set differently if needed (also usually passed from executable)
indexROOT=yourINDEXROOT;

# Set defaults, SE and for cleaning up BAM (remove PCR duplicates and low MAPQ scores)
readTYPE=SE
DeDupFLAG=TRUE
CleanFLAG=TRUE

# Set an exit method for printing an error
die() {
    printf '%s\n' "$1" >&2;
    exit 1;
}

# Add header to standard IN and ERR for debugging
echo STDOUT:
echo -e JOB_ID = $JOB_ID"\n"
hostname
echo -e Script call: submit_bowtie2_ChIP_process.sh "$@"

# Now capture arguments
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
        -g)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                indexDIR=$2
                shift
            else
                die 'ERROR: "-g" requires a non-empty option argument.'
            fi
            ;;
        --indexROOT)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                indexROOT=$2
                shift
            else
                die 'ERROR: "--indexROOT" requires a non-empty option argument.'
            fi
            ;;
        --noDeDup)   # Note: flag so no need to check for option argument
            DeDupFLAG=FALSE
            ;;
        --noClean)   # Note: flag so no need to check for option argument
            CleanFLAG=FALSE
            ;;
        --PE)   # Note: flag so no need to check for option argument
            readTYPE=PE
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
if [ -z ${indexDIR+x} ]; then die 'ERROR: -g is unset, this script needs a genome index.'; fi

# Set the working directory
rootDIR=$(pwd)

# Check to see if the output directory exists and if not create it
if [[ ! -d "${outDIR}" ]] ; then
    echo "${outDIR} does not exist, creating...."
    mkdir ${outDIR}
fi

# Continue headers for STDOUT
echo -e Bowtie2 alignment.
echo -e sample ID: ${sampleID}
echo -e input DIR: ${fileDIR}
echo -e Output DIR: ${outDIR}
echo -e Genome Index: ${indexROOT}/${indexDIR}/${indexDIR}
echo -e Read type: ${readTYPE}
echo -e PCR de-duplicate: ${DeDupFLAG}
echo -e Remove low MAPQ: ${CleanFLAG}

echo Start Time:
date

# Configure modules
# These are specific versions used on University of Edinburgh MRC HGU cluster
# Will need to be adapted to run on other systems
# Note CutAdapt is a module under python on our system
# Runs after adding python
module add python/2.7.10
module add TrimGalore/0.4.1
module add bowtie/2.3.1

# Have to add ncurses lib to run new samtools
module add ncurses/6.0
module add samtools/1.6

# Sambamba to de-duplicate
module add sambamba/0.5.9

# Copy FASTQ files to node
# Check to see if the file has extension 'fq.gz' or 'fastq.gz'
# Also needs to check if we have SE or PE and copy accordingly
# Note: in all cases file(s) are copied as 'fq.gz'
if [[ -f "${fileDIR}/${sampleID}_1.fastq.gz" && -f "${fileDIR}/${sampleID}_2.fastq.gz" && "${readTYPE}" = 'PE' ]] ;
then
    echo "Files are ${sampleID}_1.fastq.gz and ${sampleID}_2.fastq.gz"
    cp ${fileDIR}/${sampleID}_1.fastq.gz $TMPDIR/${sampleID}_1.fq.gz
    cp ${fileDIR}/${sampleID}_2.fastq.gz $TMPDIR/${sampleID}_2.fq.gz
elif [[ -f "${fileDIR}/${sampleID}_1.fq.gz" && -f "${fileDIR}/${sampleID}_2.fq.gz" && "${readTYPE}" = 'PE' ]] ;
then
    echo "Files are ${sampleID}_1.fq.gz and ${sampleID}_2.fq.gz"
    cp ${fileDIR}/${sampleID}_1.fq.gz $TMPDIR/${sampleID}_1.fq.gz
    cp ${fileDIR}/${sampleID}_2.fq.gz $TMPDIR/${sampleID}_2.fq.gz
elif [[  -f "${fileDIR}/${sampleID}.fastq.gz" && "${readTYPE}" = 'SE' ]] ;
then
    echo "File is ${sampleID}.fastq.gz"
    cp ${fileDIR}/${sampleID}.fastq.gz $TMPDIR/${sampleID}.fq.gz
elif [[  -f "${fileDIR}/${sampleID}.fq.gz" && "${readTYPE}" = 'SE' ]] ;
then
    echo "File is ${sampleID}.fq.gz"
    cp ${fileDIR}/${sampleID}.fq.gz $TMPDIR/${sampleID}.fq.gz
elif [[  -f "${fileDIR}/${sampleID}_1.fastq.gz" && "${readTYPE}" = 'SE' ]] ;
then
    echo "File is ${sampleID}_1.fastq.gz"
    cp ${fileDIR}/${sampleID}_1.fastq.gz $TMPDIR/${sampleID}.fq.gz
# SE and "sampleID_1.fq.gz"
elif [[  -f "${fileDIR}/${sampleID}_1.fq.gz" && "${readTYPE}" = 'SE' ]] ;
then
    echo "File is ${sampleID}_1.fq.gz"
    cp ${fileDIR}/${sampleID}_1.fq.gz $TMPDIR/${sampleID}.fq.gz
else
    die "ERROR: Can't find ${fileDIR}/${sampleID}[1/2] as fastq.gz or fq.gz";
fi

# Go to TMPDIR and start analysis
cd $TMPDIR

echo '************************************'
echo Start trim_galore:
echo '************************************'
date

# Run trim_galore
# Note - no read clipping set at present
# Check if we have single end or paired end
if [ "${readTYPE}" = 'SE' ] ; then
    # Single End version
    trim_galore --phred33 ${sampleID}.fq.gz
else
    # Paired end version
    trim_galore --phred33 --paired ${sampleID}_1.fq.gz ${sampleID}_2.fq.gz
    mv ${sampleID}_1_val_1.fq.gz ${sampleID}_trimmed_1.fq.gz
    mv ${sampleID}_2_val_2.fq.gz ${sampleID}_trimmed_2.fq.gz
fi

# Prior to alignment, get rid of original FASTQ to save space
if [ "${readTYPE}" = 'SE' ] ; then
    rm ${sampleID}.fq.gz
else
    rm ${sampleID}_1.fq.gz ${sampleID}_2.fq.gz
fi

echo '************************************'
echo Start bowtie2 alignment:
echo '************************************'
date

# run alignment
# Pipe this directly to a sorted BAM file
# Previously output was SAM which we then converted and sorted
# Different versions for SE and PE
# PE use --no-mixed and --no-discordant to make sure we have a clean BAM file
# Also set maximum expected fragment size to 1000bp. Fragments above this size will be discarded
# Hard-wired at present
if [ "${readTYPE}" = 'SE' ] ; then
    # Single End version
    bowtie2 --threads 3 --phred33 -N 1 -L 20 --no-unal -x ${indexROOT}/${indexDIR}/${indexDIR} -U ${sampleID}_trimmed.fq.gz | samtools view -b - | samtools sort -@ 4 - > ${sampleID}.bam
else
    # Paired end version
    bowtie2 --threads 3 --phred33 -N 1 -L 20 --no-unal --no-mixed --no-discordant -X 1000 -x ${indexROOT}/${indexDIR}/${indexDIR} -1 ${sampleID}_trimmed_1.fq.gz -2 ${sampleID}_trimmed_2.fq.gz | samtools view -b - | samtools sort -@ 4 - > ${sampleID}.bam
fi

# Following alignment, get rid of trimmed FASTQ to save space
if [ "${readTYPE}" = 'SE' ] ; then
    rm ${sampleID}_trimmed.fq.gz
else
    rm ${sampleID}_trimmed_1.fq.gz ${sampleID}_trimmed_2.fq.gz
fi

echo '************************************'
echo Start duplicate removal:
echo '************************************'
date

# Run SAMBAMBA to remove duplicates if specified
# No need to sort as bam is sorted

# Check if we want to do this
if [ "${DeDupFLAG}" = 'TRUE' ] ; then
    # Do the de-duplication (default)
    sambamba markdup --nthreads 4 -r ${sampleID}.bam ${sampleID}_DeDUP.bam
else
    # --noDeDUP is specified then we just copy the bam file to keep things consistent downstream
    cp ${sampleID}.bam ${sampleID}_DeDUP.bam
fi

echo '************************************'
echo Start removal of low mapping quality reads:
echo '************************************'
date

# Run samtools to remove reads with poor mapping quality (eg multi-aligning) if specified

# Check if we want to do this
if [ "${CleanFLAG}" = 'TRUE' ] ; then
    # Do the removal of low MAPQ reads (default)
    samtools view -bq 10 ${sampleID}_DeDUP.bam > ${sampleID}_DeDUP_excl.bam
else
    # --noClean is specified then we just copy the bam file to keep things consistent downstream
    cp ${sampleID}_DeDUP.bam ${sampleID}_DeDUP_excl.bam
fi

# Generate an index
samtools index -b ${sampleID}_DeDUP_excl.bam

echo '************************************'
echo Get number of reads in the file after each stage:
echo '************************************'
date

# Note this may reveal no change if --noDeDUP and/or --noClean are set
# Write to report files
# SE version:
if [ "${readTYPE}" = 'SE' ] ; then
  echo Number of reads after alignment:
  sambamba view -c --nthreads 4 ${sampleID}.bam
  sambamba view -c --nthreads 4 ${sampleID}.bam > total_align_count.txt

  echo Number of reads after deduplication and before low mapping quality removal:
  sambamba view -c --nthreads 4 ${sampleID}_DeDUP.bam
  sambamba view -c --nthreads 4 ${sampleID}_DeDUP.bam > total_dedup_count.txt

  echo Number of reads after deduplication and low mapping quality removal:
  sambamba view -c --nthreads 4 ${sampleID}_DeDUP_excl.bam
  sambamba view -c --nthreads 4 ${sampleID}_DeDUP_excl.bam > total_dedup_clean_count.txt
echo
# PE version - divide the counts by 2
else
  echo Number of reads after alignment:
  sambamba view -c --nthreads 4 ${sampleID}.bam | awk '{print $1/2}'
  sambamba view -c --nthreads 4 ${sampleID}.bam | awk '{print $1/2}' > total_align_count.txt

  echo Number of reads after deduplication and before low mapping quality removal:
  sambamba view -c --nthreads 4 ${sampleID}_DeDUP.bam | awk '{print $1/2}'
  sambamba view -c --nthreads 4 ${sampleID}_DeDUP.bam | awk '{print $1/2}' > total_dedup_count.txt

  echo Number of reads after deduplication and low mapping quality removal:
  sambamba view -c --nthreads 4 ${sampleID}_DeDUP_excl.bam | awk '{print $1/2}'
  sambamba view -c --nthreads 4 ${sampleID}_DeDUP_excl.bam | awk '{print $1/2}' > total_dedup_clean_count.txt
  echo
fi

echo '************************************'
echo Copy output back
echo '************************************'
date

# Copy output back to eddie
cp ${sampleID}_DeDUP_excl.bam ${rootDIR}/${outDIR}/${sampleID}.bam
cp ${sampleID}_DeDUP_excl.bam.bai ${rootDIR}/${outDIR}/${sampleID}.bam.bai

# Also copy back the counts
cp total_align_count.txt ${rootDIR}/${outDIR}/${sampleID}_alignCount.txt
cp total_dedup_count.txt ${rootDIR}/${outDIR}/${sampleID}_dedupCount.txt
cp total_dedup_clean_count.txt ${rootDIR}/${outDIR}/${sampleID}_dedupCleanCount.txt

# Debugging ls
echo '************************************'
echo End file list on TMPDIR:
echo '************************************'
ls -1

# Return to root
cd ${rootDIR}
# Add end time stamp to STDOUT and STDERR
echo '************************************'
echo End Time:
echo '************************************'
date


echo '************************************'
echo output files:
echo '************************************'
echo ${sampleID}.bam
echo ${sampleID}.bam.bai
echo ${sampleID}_alignCount.txt
echo ${sampleID}_dedupCount.txt
echo ${sampleID}_finalCount.txt


# Move STDOUT file based on sample ID for easier diagnosis of problems
# Note should only produce a STDOUT (setting qsub option -j yes)
mv ${JOB_ID}_out ${outDIR}/${sampleID}_${JOB_ID}_bowtie2_process_stdout
