#!/bin/sh

# Options for running on a cluster need to be specified
# Have removed the options required for University of Edinburgh and MRC HGU compute cluster

# Usage
# qsub submit_RRBS_process.sh -d fileDIR -i sampleID -o outDIR -g indexDIR [--scriptDIR ~ --indexROOT]

# Only designed to work on PE phred33 files only at present which is default that we recieved.
# Includes NuGen diversity base trimming and deduplication based on index UMI

# Will look inside the --indexROOT directory for genome indexes
# Was set as default for our cluster
# Will need to be modified for other systems
# Will check for CG reference file inside this in folder: CG_locations (this can't be modified)
# CG reference file is a BED file of all CG locations in the genome

# Note: This version requested lots of resources to ensure that jobs don't die
# Also uses parallel Bismark (2 cores align, 3 cores extract - note they already use multiple cores per process, so expect ~10 used for align and 9 for extract)

# Clipping is hard-wired to our best settings

# Arguments, these are all required
# -i sampleID, input file ID, note without '_1' or '_2' or fastq gz extension
# -d fileDIR, data directory, where FASTQ files are stored
# -o outDIR, directory to write output files to
# -g indexDIR, name of directory containing the index for alignment

# Optional command line parameter to set script root DIR
# --scriptDIR - set root DIR of script
# Required to locate other scripts needed for processing

# Optional command line parameter to set path to folder containing indexes:
# indexROOT

# Need to set the index root to use on other systems
# Can be set differently if needed (also usually passed from executable)
indexROOT=yourINDEXROOT;

# Set root directory of scripts by default to ~
# This is usually passed by the executable to launch the qsub job
scriptDIR=~;

# Set an exit method for printing an error
die() {
    printf '%s\n' "$1" >&2;
    exit 1;
}

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
        --scriptDIR)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                scriptDIR=$2
                shift
            else
                die 'ERROR: "--scriptDIR" requires a non-empty option argument.'
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
if [ -z ${indexDIR+x} ]; then die 'ERROR: -g is unset, this script needs a bisulfite converted index.'; fi

# Set the working directory
rootDIR=$(pwd)

# Check to see if the output directory exists and if not create it
if [[ ! -d "${outDIR}" ]] ; then
    echo "${outDIR} does not exist, creating...."
    mkdir ${outDIR}
fi

# Have to find reference CG file:
# In ${indexROOT}/${indexDIR}/CG_locations/
# Note this current script will cause a problem if more than one file is present
# Arbitrarely pick the first one in the ls to avoid this problem a bit
# Name is non standard but should always begin with cg.loc.
fileCGRef=$(ls ${indexROOT}/${indexDIR}/CG_locations/cg.loc.* | head -n 1)

# Add headers to STDOUT and STDERR
echo STDOUT:
echo -e JOB_ID = $JOB_ID"\n"
hostname
echo -e RRBS PE processing alingment to BEDgraph by Strand.
echo -e FASTQ1: ${fileDIR}/${sampleID}_1.fq.gz
echo -e FASTQ2: ${fileDIR}/${sampleID}_2.fq.gz
echo -e Output DIR: ${outDIR}
echo -e Genome Index: ${indexROOT}/${indexDIR}
echo -e CG reference file: ${indexROOT}/${indexDIR}/CG_locations/${fileCGRef}

echo Start Time:
date

# Configure modules
# These are specific versions used on University of Edinburgh MRC HGU cluster
# Will need to be adapted to run on other systems
# Note CutAdapt is a module under python on our system
# Runs after adding python
module add bismark/0.18.1
module add python/2.7.10
module add TrimGalore/0.4.1
module add bowtie/2.3.1

# Have to add ncurses and htslib to run new samtools
# Needed by bismark
module add htslib/1.6
module add ncurses/6.0
module add samtools/1.6

# Add R to use Mbias script at end
module add R/3.6.0

# Copy FASTQ files to node
cp ${fileDIR}/${sampleID}_1.fq.gz $TMPDIR/${sampleID}_1.fq.gz
cp ${fileDIR}/${sampleID}_2.fq.gz $TMPDIR/${sampleID}_2.fq.gz

# Go to TMPDIR and start analysis
cd $TMPDIR

echo '************************************'
echo Start trim_galore:
echo '************************************'
date

# Run trim_galore
trim_galore  --phred33 --paired --adapter AGATCGGAAGAGC --adapter2 AAATCAAAAAAAC ${sampleID}_1.fq.gz ${sampleID}_2.fq.gz

echo '************************************'
echo Start Diversity Trimming:
echo '************************************'
date

# Run the diversity trimming script
python ${scriptDIR}/eddie3_RRBSprocess/NuGen_trimRRBSdiversityAdaptCustomers_v1_260916.py -1 "${sampleID}_1_val_1.fq.gz" -2 "${sampleID}_2_val_2.fq.gz"

echo '************************************'
echo Check file list: Post Diversity trimming
echo '************************************'
ls -1
echo


echo '************************************'
echo Start bismark alignment:
echo '************************************'
date

# run alignment
bismark --multicore 2 --phred33-quals -N 0 -L 20 ${indexROOT}/${indexDIR} -1 ${sampleID}_1_val_1.fq_trimmed.fq.gz -2 ${sampleID}_2_val_2.fq_trimmed.fq.gz


echo '************************************'
echo Check file list: Post Bismark alignment
echo '************************************'
ls -1
echo

echo '************************************'
echo Start alignment de-duplication:
echo '************************************'
date

# Re-name output for easier processing downstream
mv ${sampleID}_1_val_1.fq_trimmed_bismark_bt2_pe.bam ${sampleID}_bismark.bam

# Sort the alignment based on position
samtools sort ${sampleID}_bismark.bam > ${sampleID}_bismark.sort.bam

# De-duplicate
# Note: set with default settings for UMI co-ordinates 6bp from end and 6bp long
python ${scriptDIR}/eddie3_RRBSprocess/NuGen_nudup_RRBS_RemoveDups_v1_090217.py -o ${sampleID}_RMDUP --paired-end --rmdup-only ${sampleID}_bismark.sort.bam

# Re-sort on read name to keep correct for Bismark
samtools sort -n ${sampleID}_RMDUP.sorted.dedup.bam > ${sampleID}_RMDUP_NAMESORT_bismark.bam

# Note report file is:
# ${sampleID}_RMDUP_dup_log.txt

echo '************************************'
echo Check file list: Post De-duplication
echo '************************************'
ls -1
echo


echo '************************************'
echo Start bismark methylation extractor:
echo '************************************'
date

# Do bismark methylation extraction
# As SAMtools is running, output will be bam
# Input file name: ${sampleID}_1_val_1_bismark_bt2_pe.bam
bismark_methylation_extractor --multicore 3 -p --no_overlap --no_header ${sampleID}_RMDUP_NAMESORT_bismark.bam

# Tidy up output
# Get rid of all non-CpG output
rm CHG_* CHH_*

# Parse Bismark methylation extractor outputs to BEDgraph
# File IDs for parsing to BEDgraph
# CpG_OT_${sampleID}_RMDUP_NAMESORT_bismark.txt
# CpG_OB_${sampleID}_RMDUP_NAMESORT_bismark.txt

echo '************************************'
echo Error checking file list: Post Bismark Methylation Extractor
echo '************************************'
ls -1
echo

echo '************************************'
echo Start Bismark BEDgraph conversions:
echo '************************************'

date

bismark2bedGraph -o CpG_OT_${sampleID}_BEDgraph.txt --buffer_size 50% --ample_memory CpG_OT_${sampleID}_RMDUP_NAMESORT_bismark.txt

echo '************************************'
echo Done Bismark BEDgraph conversion OT:
echo '************************************'
date

bismark2bedGraph -o CpG_OB_${sampleID}_BEDgraph.txt --buffer_size 50% --ample_memory CpG_OB_${sampleID}_RMDUP_NAMESORT_bismark.txt

echo '************************************'
echo Done Bismark BEDgraph conversion OB:
echo '************************************'
date

# Now we want to parse the output to give CpG interval and report methylated + total
# Output files (want coverage ones)
# Bismark to BEDgraph now does: CpG_OT_test_BEDgraph.txt.gz and CpG_OT_test_BEDgraph.txt.bismark.cov.gz
# CpG_OT_${sampleID}_BEDgraph.txt
# CpG_OB_${sampleID}_BEDgraph.txt
# CpG_OT_${sampleID}_BEDgraph.txt.bismark.cov
# CpG_OB_${sampleID}_BEDgraph.txt.bismark.cov
echo
ls
echo

echo '************************************'
echo Start parsing BEDgraph outputs:
echo '************************************'
date

# Need to decompress first
echo
echo Start GZIP time: $(date)
gzip -dc CpG_OT_${sampleID}_BEDgraph.txt.gz.bismark.cov.gz > CpG_OT_${sampleID}_BEDgraph.txt.bismark.cov
gzip -dc CpG_OB_${sampleID}_BEDgraph.txt.gz.bismark.cov.gz > CpG_OB_${sampleID}_BEDgraph.txt.bismark.cov
echo End GZIP time: $(date)
echo

python ${scriptDIR}/eddie3_RRBSprocess/prepare_strand_match_BEDgraphCov.py CpG_OT_${sampleID}_BEDgraph.txt.bismark.cov CpG_OT_${sampleID}_BEDgraph_parsed.txt OT
python ${scriptDIR}/eddie3_RRBSprocess/prepare_strand_match_BEDgraphCov.py CpG_OB_${sampleID}_BEDgraph.txt.bismark.cov CpG_OB_${sampleID}_BEDgraph_parsed.txt OB

echo '************************************'
echo Done python parsing:
echo '************************************'
date

# Now match up files to all CpGs in genome
# Match against reference generated by C script
# Use AWK script that prints CpGs that are found in data with values and 0 for T and M if not
# Files to be parsed:
# CpG_OT_${sampleID}_BEDgraph_parsed.txt
# CpG_OB_${sampleID}_BEDgraph_parsed.txt

# Now match data to the reference CG file
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR{a[$1,$2,$3]=$4 FS $5;next} {if (($1,$2,$3) in a) {print $0,a[$1,$2,$3]} else {print $0,0,0} }' CpG_OT_${sampleID}_BEDgraph_parsed.txt ${fileCGRef} > CpG_OT_${sampleID}_BEDgraph_allCpG.txt

awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR{a[$1,$2,$3]=$4 FS $5;next} {if (($1,$2,$3) in a) {print $0,a[$1,$2,$3]} else {print $0,0,0} }' CpG_OB_${sampleID}_BEDgraph_parsed.txt ${fileCGRef} > CpG_OB_${sampleID}_BEDgraph_allCpG.txt

echo '************************************'
echo Done matching against reference CpG:
echo '************************************'
date

# Now sort the resulting files
# Do by chr and start
# Files to be sorted
# CpG_OT_${sampleID}_BEDgraph_allCpG.txt
# CpG_OB_${sampleID}_BEDgraph_allCpG.txt
sort -k1,1 -k2,2n CpG_OT_${sampleID}_BEDgraph_allCpG.txt > CpG_OT_${sampleID}_BEDgraph_allCpG_sorted.txt
sort -k1,1 -k2,2n CpG_OB_${sampleID}_BEDgraph_allCpG.txt > CpG_OB_${sampleID}_BEDgraph_allCpG_sorted.txt

echo '************************************'
echo Done sorting files:
echo '************************************'
date

echo '************************************'
echo Check number of lines in outputs:
echo '************************************'
echo No. lines OT:
wc -l CpG_OT_${sampleID}_BEDgraph_allCpG_sorted.txt
echo No. lines OB:
wc -l CpG_OB_${sampleID}_BEDgraph_allCpG_sorted.txt
echo ''

# Combine to get one file per input with forward and reverse strand in one
# Note this will duplicate the chr, start, end colums
paste CpG_OT_${sampleID}_BEDgraph_allCpG_sorted.txt CpG_OB_${sampleID}_BEDgraph_allCpG_sorted.txt > CpG_${sampleID}_BEDgraph_allCpG_sorted.txt
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,$4,$5,$9,$10}' CpG_${sampleID}_BEDgraph_allCpG_sorted.txt > CpG_${sampleID}_BEDgraph_allCpG_sorted_clean.txt

echo '************************************'
echo Done combining OT and OB:
echo '************************************'
date

# Compress the output
gzip CpG_${sampleID}_BEDgraph_allCpG_sorted_clean.txt

# Now copy back the data BEDgraph file
cp CpG_${sampleID}_BEDgraph_allCpG_sorted_clean.txt.gz ${rootDIR}/${outDIR}/CpG_${sampleID}_BEDgraph.txt.gz

# Copy back the Trim Galore reports
cp ${sampleID}_1.fq.gz_trimming_report.txt ${rootDIR}/${outDIR}/${sampleID}_1_trimming_report.txt
cp ${sampleID}_2.fq.gz_trimming_report.txt ${rootDIR}/${outDIR}/${sampleID}_2_trimming_report.txt

# Copy back the Bismark alignment reports
cp ${sampleID}_1_val_1.fq_trimmed_bismark_bt2_PE_report.txt ${rootDIR}/${outDIR}/${sampleID}_bismark_align_report.txt

# De-duplication report
cp ${sampleID}_RMDUP_dup_log.txt  ${rootDIR}/${outDIR}/${sampleID}_RRBS_read_duplication_report.txt

# Do the M-bias plots
# Have to change M-bias output file name first
mv ${sampleID}_RMDUP_NAMESORT_bismark.M-bias.txt ${sampleID}_bismark_M-bias.txt
Rscript --vanilla ${scriptDIR}/eddie3_RRBSprocess/draw_RRBS_Mbias_plot_PE.R ${sampleID}

# Copy Methylation extraction output back to eddie space along with the M-bias plot
cp ${sampleID}_RMDUP_NAMESORT_bismark_splitting_report.txt ${rootDIR}/${outDIR}/${sampleID}_bismark_extract_report.txt
cp ${sampleID}_bismark_M-bias.txt ${rootDIR}/${outDIR}/${sampleID}_bismark_Mbias.txt
cp ${sampleID}_Mbias_plots.pdf ${rootDIR}/${outDIR}/${sampleID}_Mbias_plots.pdf

# Generate an autosomal CpG report for coverage, etc
# This part doesn't care if it is SE or PE
gzip -dc CpG_${sampleID}_BEDgraph_allCpG_sorted_clean.txt.gz | python ${scriptDIR}/eddie3_RRBSprocess/get_BEDgraph_CpG_report.py > ${sampleID}_autosomal_CpG_report.txt
cp ${sampleID}_autosomal_CpG_report.txt ${rootDIR}/${outDIR}/${sampleID}_autosomal_CpG_report.txt

# Now get a report of methylation on Lambda
# Note: will not throw up an error as long as chrLambda was included in reference
# This part doesn't care if it is SE or PE
gzip -dc CpG_${sampleID}_BEDgraph_allCpG_sorted_clean.txt.gz | grep -e 'chrLambda_J0245' | awk '{ t += $4 + $6; m += $5 + $7; } END { print "Lambda T:", t; print "Lambda M:",  m; }' > ${sampleID}_lambda_CpG_report.txt
cp  ${sampleID}_lambda_CpG_report.txt  ${rootDIR}/${outDIR}/${sampleID}_lambda_CpG_report.txt

echo '************************************'
echo Final check contents of TMPDIR at end of script:
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
echo Results:
echo CpG_${sampleID}_BEDgraph.txt.gz
echo Reports:
echo ${sampleID}_[1/2]_trimming_report.txt
echo ${sampleID}_bismark_align_report.txt
echo ${sampleID}_RRBS_read_duplication_report.txt
echo ${sampleID}_bismark_extract_report.txt
echo ${sampleID}_bismark_Mbias.txt
echo ${sampleID}_Mbias_plots.pdf
echo ${sampleID}_autosomal_CpG_report.txt
echo ${sampleID}_lambda_CpG_report.txt

# Note main output BEDgraph file has following columns:
# chr, start, end, T.for, M.for, T.rev, M.rev

# Move STDOUT and STDERR to file based on sample ID for easier diagnosis of problems
# Single file as -j is set
mv ${JOB_ID}_out ${outDIR}/${sampleID}_${JOB_ID}_RRBS_process_stdout
