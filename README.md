## Custom scripts used for Masalmeh et al 2019
#### De novo DNA methyltransferase activity in colorectal cancer is directed towards H3K36me3 marked CpG islands
Roza H. Ali Masalmeh, Francesca Taglini, Cristina Rubio-Ramon, Kamila I. Musialik, Jonathan Higham, Hazel Davidson-Smith, Ioannis Kafetzopoulos, Kamila P. Pawlicka, Hannah M. Finan, Richard Clark, Jimi Wills, Andrew J. Finch, Lee Murphy, Duncan Sproul
https://doi.org/10.1101/676346

Created by Duncan Sproul

Includes 3 sets of scripts used for different parts of the paper:
- *bsPCR_scripts:* R scripts used to process bsPCR data output from BISMA and BiQ <br />
- *RRBS_scripts* Scripts used to process RRBS on the University of Edinburgh MRC HGU computer cluster <br />
- *ChIPseq_scripts* Script used to process ChIP-seq on the University of Edinburgh MRC HGU computer cluster <br />

Note: The RRBS and ChIP-seq scripts were set up to run on the University of Edinburgh MRC HGU cluster. References to programs and directories specific to this cluster have been removed and these will need to be modified to be set up for another system.

### bsPCR_scripts
These R scripts both work on the HTML output of BISMA and BiQ. They produce a simple lollipop diagram from the clones as well as a nicely formatted table.

#### draw_bsPCR_lollipops_BISMA.R <br />
Usage: *Rscript --vanilla draw_bsPCR_lollipops_BISMA.R sampleID* <br />

Will process look for sampleID.html output file from BISMA and output a PDF image of the clones and CpGs where (*sampleID_plot.pdf*): <br />
Unmethylated = White <br />
Methylated = Black <br />
No data = Grey <br />
Will also output a tab-delimited text table of the data where (*sampleID_table.txt*): <br />
Unmethylated = 0 <br />
Methylated = 1 <br />
No data = NA <br />

The script requires R to be present in the path. Ie if you type R on the command line, it loads R. <br />

*Input file names:* <br />
By default expects to find a file sampleID.html in the current directory. <br />

Example: <br />
`Rscript --vanilla draw_bsPCR_lollipops_BISMA.R SFRP1_HCT116_WT` <br />

Will produce *SFRP1_HCT116_WT_plot.pdf* and *SFRP1_HCT116_WT_table.txt*. <br />

Note: this script does not check for any problems so error messages may not be particularly useful. <br />

#### draw_bsPCR_lollipops_BIQ.R <br />
Usage: *Rscript --vanilla draw_bsPCR_lollipops_BIQA.R sampleID* <br />

Will process look for sampleID.html output file from BISMA and output a PDF image of the clones and CpGs where (*sampleID_plot.pdf*): <br />
Unmethylated = White <br />
Methylated = Black <br />
No data = Grey <br />
Will also output a tab-delimited text table of the data where (*sampleID_table.txt*): <br />
Unmethylated = 0 <br />
Methylated = 1 <br />
No data = NA <br />

The script requires R to be present in the path. Ie if you type R on the command line, it loads R. <br />

*Input file names:* <br />
By default expects to find a file sampleID.html in the current directory. <br />

Example: <br />
`Rscript --vanilla draw_bsPCR_lollipops_BIQ.R SFRP1_HCT116_WT` <br />

Will produce *SFRP1_HCT116_WT_plot.pdf* and *SFRP1_HCT116_WT_table.txt*. <br />

Note: this script does not check for any problems so error messages may not be particularly useful. <br />

### RRBS_scripts
These scripts process raw multiplexed FASTQ files from Illumina sequenced NuGen RRBS into a processed BEDgraph like file giving the coverage on the forward and reverse strand for every CpG in the genome. <br />

The submitted scripts should be run in the following order: <br />
- *submit_RRBS_demultiplex_NuGen.sh:* Demultiplexes the raw FASTQ files based on the index read in the read header. <br />
- *submit_fastqMerge_RRBS.sh:* Combines multiple FASTQ files for a single sample. <br />
- *submit_RRBS_process.sh:* Aligns and processes the FASTQ file using Bismark. <br />

All other scripts are called and required from within these scripts. <br />

#### submit_RRBS_demultiplex_NuGen.sh <br />
Usage: *qsub submit_RRBS_demultiplex_NuGen.sh -d fileDIR -i fileID -m indexFILE -o outDIR [--scriptDIR ~* <br />

*-i* fileID, input file ID, note without '_1' or '_2' or fastq gz extension <br />
*-d* fileDIR, data directory, where FASTQ files are stored <br />
*-o* outDIR, directory to write output files to <br />
*-m* A tab-delimited text file giving the index sequence and sample IDs to look for <br />
*--scriptDIR* directory to find the required scripts in, defaults to ~.<br />
Assumes you have paired end files and will process both. <br />
Will only look for exact matches to the index sequences. <br />

#### submit_fastqMerge_RRBS.sh <br />
Usage: *qsub submit_fastqMerge_RRBS.sh -d fileDIR -i sampleID -o outDIR* <br />

*-i* sampleID, sample ID, note without '_1' or '_2' or fastq gz extension <br />
*-d* fileDIR, data directory, where FASTQ files are stored <br />
*-o* outDIR, directory to write output files to <br />

Assumes that all files are present in one directory and have two FASTQ files per lane (_1 and_2). <br />
.fq.gz extension assumed also. <br />
Searches for L00N in name also as this is the format our files were supplied in. <br />

#### submit_RRBS_process.sh <br />
Usage: *qsub submit_RRBS_process.sh -d fileDIR -i sampleID -o outDIR -g indexDIR [--scriptDIR ~ --indexROOT]* <br />

*-i* sampleID, input file ID, note without '_1' or '_2' or fastq gz extension  <br />
*-d* fileDIR, data directory, where FASTQ files are stored  <br />
*-o* outDIR, directory to write output files to  <br />
*-g* indexDIR, name of directory containing the index for alignment  <br />

Optional command line parameter to set script root DIR <br />
*--scriptDIR* - set root DIR of script <br />
Required to locate other scripts needed for processing <br />

Optional command line parameter to set path to folder containing indexes: <br />
*--indexROOT*

Only designed to work on PE phred33 files only at present which is default that we recieved <br />
Includes NuGen diversity base trimming and deduplication based on index UMI <br />

Will look inside the --indexROOT directory for genome indexes <br />
Was set as default for our cluster <br />
Will need to be modified for other systems <br />
Will check for CG reference file inside this in folder: CG_locations (this can't be modified) <br />
CG reference file is a BED file of all CG locations in the genome <br />

Note: This version requested lots of resources to ensure that jobs don't die <br />
Also uses parallel Bismark (2 cores align, 3 cores extract - note they already use multiple cores per process, so expect ~10 used for align and 9 for extract) <br />

Clipping is hard-wired to our best settings <br />

Output file is a modified BEDgraph format giving coverage at every CG in the genome on both strands (set by the CG index file): <br />
chr, start, end, T.for, M.for, T.rev, M.rev <br />
Where T is total coverage and M is methylated coverage <br />

### ChIPseq_scripts
These scripts process raw FASTQ files from Illumina sequenced ChIP-seq to produce aligned, deduplicated BAM files which remove multimapping reads. <br />

#### submit_bowtie2_ChIP_process.sh <br />
Usage: *qsub submit_bowtie2_ChIP_process.sh -d fileDIR -i sampleID -o outDIR -g indexDIR [--PE --noDeDup --noClean --indexROOT]* <br />

*-i* sampleID, input file ID, note without '_1' or '_2' or fastq gz extension <br />
*-d* fileDIR, data directory, where FASTQ files are stored <br />
*-o* outDIR, directory to write output files to <br />
*-g* indexDIR, name of directory containing the index for alignment <br />

Optional command line parameters <br />
*--PE*, use PE FASTQ files <br />
*--noDeDup*, don't remove PCR duplicates <br />
*--noClean*, don't remove low quality mapping reads/fragments  <br />
Will look inside the *--indexROOT* directory for genome indexes  <br />
Was set as default for our cluster  <br />
Will need to be modified for other systems  <br />
Assumes index has same name as the indexDIR  <br />

Designed to work on SE phred33 files by default  <br />
Will take a FASTQ file and align it against a given genome index  <br />
Following alignment it will strip multimapping reads and PCR duplicates before sorting and indexing the output BAM file by default  <br />
This cleanup can be turned off with flags if needed  <br />

Note: no read clipping and phred 33 hardwired and can't be changed at present  <br />
For PE is hard wired at maximum fragment size of 1000  <br />
