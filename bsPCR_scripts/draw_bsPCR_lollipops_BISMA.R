# Script to read BISMA parsed bsPCR output
# Generate a neat plot of the alleles methylation status
# Also output a nicely formatted table of the data

#!/usr/bin/env Rscript

# Collect the file ID from the command line
# Also the output ID
args <- commandArgs(trailingOnly=TRUE);
sample.id <- args[1];

# Set the input file ID
file.in <- paste(sample.id, "html", sep=".");

# Read in the file
data.in <- readLines(file.in);

# Get the line number for the beginning of the table we want
table.head.line.no <- grep("Methylation data for these CpG dinucleotides", data.in);

# Get the CpG numbers and parse
cpg.no <- strsplit(data.in[table.head.line.no + 1], split="\t")[[1]];
cpg.no <- as.numeric(sub("^CpG_", "", cpg.no, perl=TRUE)[2:length(cpg.no)]);

# Now we need to get the lines with the table of information
# grep for [1]SeqID\t
meth.status.table.line.nos <- grep('^\\[[0-9]+\\].+\t', data.in, perl=TRUE)[grep('^\\[[0-9]+\\].+\t', data.in, perl=TRUE) > table.head.line.no];

# Split the table
meth.status.table <- sapply(data.in[meth.status.table.line.nos], strsplit, split="\t");

# Re-arrange the table
# Suppress warnings here as x will be converted to NA
options(warn=-1);
meth.status.table <- data.frame(t(sapply(meth.status.table, function(x) {as.vector(x[2:length(x)], mode="numeric");})), stringsAsFactors=FALSE);
options(warn=0);

# Clean up the rownames and colnames
rownames(meth.status.table) <- sub("^\\[[0-9]+\\]", "", sub("\t.+$", "", rownames(meth.status.table), perl=TRUE), perl=TRUE);
colnames(meth.status.table) <- paste("CpG", cpg.no, sep="_");

# Write out the data table for analyses
write.table(meth.status.table, paste(sample.id,"_table.txt", sep="_"), sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

# Now we want to draw the plot
# Set up the table with NA being 3, (1=Unmeth and 2=Meth)
meth.status.table.plot <- meth.status.table + 1;
meth.status.table.plot[is.na(meth.status.table)] <- 3;

no.cpgs <- ncol(meth.status.table.plot);
no.clones <- nrow(meth.status.table.plot);

# Draw the plot
# Initiate
pdf(paste(sample.id,"_plot.pdf", sep="_"), useDingbats=FALSE);
par(mar=c(0,0,0,0));
plot(1,1, xlim=c(1,no.cpgs), ylim=c(1,no.clones), type="n", axes=FALSE, xlab="", ylab="", main="");
# Add the lines
invisible(sapply(c(1:no.clones), function(y, x=no.cpgs) {lines(x=c(1,x), y=c(y,y))}));
# Add the points
invisible(sapply(c(1:no.clones), function(y, x=no.cpgs, table=meth.status.table.plot) {points(x=c(1:no.cpgs), y=rep(y, times=no.cpgs),
  pch=21, cex=2, bg=c("white","black", "grey")[as.vector(table[y,], mode="numeric")])}));
invisible(dev.off());

q(save = "no");
