#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Usage
# Rscript --vanilla get_index_histogram_v1_060217.R [index_file hist_file]

# Get the file names to read
index_file <- args[1]
hist_file <- args[2]

# Read in the index file
index.file.in <- read.table(index_file, sep="\t",
col.names=c("index","sample"), stringsAsFactors=FALSE)

# Read in the histogram
data.in <- read.table(hist_file,
col.names=c("count","index"), stringsAsFactors=FALSE)

# Get frequency
data.in$freq <- data.in$count/sum(data.in$count)

# Split into those that are in the index file and not
data.index <- data.in[data.in$index %in% index.file.in$index,]
data.other <- data.in[!(data.in$index %in% index.file.in$index),]

# Now split remainder by frequency >= 0.01 and < 0.01
data.other.highCount <- data.other[data.other$freq >= 0.01,]
data.other.lowCount <- data.other[data.other$freq < 0.01,]

# Now create a combined table
data.other.lowCount.sum <- data.frame(sum(data.other.lowCount$count),
"Any Others", sum(data.other.lowCount$freq), stringsAsFactors=FALSE)
colnames(data.other.lowCount.sum) <- c("count","index","freq")
data.hist <- rbind(data.index, data.other.highCount, data.other.lowCount.sum)
data.hist <- data.hist[order(data.hist$freq),]

# Now create an output graph
png(filename = paste(sub("Counts.txt","",hist_file), "Plot.png", sep=""),
width = 600, height=100+(50*nrow(data.hist)))

# Formatting to allow the names
longest.name <- max(nchar(data.hist$index))
par(mar=c(5.1,4+(longest.name/2.5),4.1,2.1))

# Now set the colours
hist.cols <- rep('grey', nrow(data.hist))
hist.cols[data.hist$index %in% index.file.in$index] <- "#377EB8"

# Sum up how much the sample indexes explain
sample.freq.sum <- sum(data.index$freq)

# Do the plot
barplot(
data.hist$freq*100,
names.arg=data.hist$index,
horiz=TRUE,
las=1,
col=hist.cols,
xlab="Percentage of reads",
main=paste("Index Frequencies: ", sub("_IndexCounts.txt","",hist_file), ".fq.gz", sep=""))
box()

mtext(paste("Samples explain ", round(sample.freq.sum*100, digits = 2),"% of Reads", sep=""))

dev.off()
