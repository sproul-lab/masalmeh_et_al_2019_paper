# Script to read in single end M-bias output from Bismark
# Outputs a PDF of the M-bias data in each nucleotide context

#!/usr/bin/env Rscript

# Collect the sample ID from the command line
args <- commandArgs(trailingOnly=TRUE);
sample.id <- args[1];

# Function to read in read in the data
read.mbias.file <- function(x)
{
    data.in <- readLines(x);
    # Get the line that each header is on
    line.no.cpg.r1 <- grep("CpG context \\(R1\\)", data.in);
    line.no.cpg.r2 <- grep("CpG context \\(R2\\)", data.in);
    line.no.chg.r1 <- grep("CHG context \\(R1\\)", data.in);
    line.no.chg.r2 <- grep("CHG context \\(R2\\)", data.in);
    line.no.chh.r1 <- grep("CHH context \\(R1\\)", data.in);
    line.no.chh.r2 <- grep("CHH context \\(R2\\)", data.in);
    # Get the line numbers of the spaces
    line.no.space <- grep("^[A-Za-z0-9]|=", data.in, invert=TRUE);

    # Now get the data
    # CpG Read 1
    data.cpg.r1 <- do.call(rbind, sapply(data.in[(line.no.cpg.r1+3):(line.no.space[1]-1)], strsplit, split="\t"));
    class(data.cpg.r1) <- "numeric";
    rownames(data.cpg.r1) <- NULL;
    colnames(data.cpg.r1) <- strsplit(data.in[(line.no.cpg.r1+2)], split="\t")[[1]];
    data.cpg.r1 <- data.frame(data.cpg.r1);
    
    # CpG Read 2
    data.cpg.r2 <- do.call(rbind, sapply(data.in[(line.no.cpg.r2+3):(line.no.space[4]-1)], strsplit, split="\t"));
    class(data.cpg.r2) <- "numeric";
    rownames(data.cpg.r2) <- NULL;
    colnames(data.cpg.r2) <- strsplit(data.in[(line.no.cpg.r2+2)], split="\t")[[1]];
    data.cpg.r2 <- data.frame(data.cpg.r2);
    
    # CHG Read 1
    data.chg.r1 <- do.call(rbind, sapply(data.in[(line.no.chg.r1+3):(line.no.space[2]-1)], strsplit, split="\t"));
    class(data.chg.r1) <- "numeric";
    rownames(data.chg.r1) <- NULL;
    colnames(data.chg.r1) <- strsplit(data.in[(line.no.chg.r1+2)], split="\t")[[1]];
    data.chg.r1 <- data.frame(data.chg.r1);
    
    # CHG Read 2
    data.chg.r2 <- do.call(rbind, sapply(data.in[(line.no.chg.r2+3):(line.no.space[5]-1)], strsplit, split="\t"));
    class(data.chg.r2) <- "numeric";
    rownames(data.chg.r2) <- NULL;
    colnames(data.chg.r2) <- strsplit(data.in[(line.no.chg.r2+2)], split="\t")[[1]];
    data.chg.r2 <- data.frame(data.chg.r2);
    
    # CHH Read 1
    data.chh.r1 <- do.call(rbind, sapply(data.in[(line.no.chh.r1+3):(line.no.space[3]-1)], strsplit, split="\t"));
    class(data.chh.r1) <- "numeric";
    rownames(data.chh.r1) <- NULL;
    colnames(data.chh.r1) <- strsplit(data.in[(line.no.chh.r1+2)], split="\t")[[1]];
    data.chh.r1 <- data.frame(data.chh.r1);
    
    # CHG Read 2
    data.chh.r2 <- do.call(rbind, sapply(data.in[(line.no.chh.r2+3):(line.no.space[6]-1)], strsplit, split="\t"));
    class(data.chh.r2) <- "numeric";
    rownames(data.chh.r2) <- NULL;
    colnames(data.chh.r2) <- strsplit(data.in[(line.no.chh.r2+2)], split="\t")[[1]];
    data.chh.r2 <- data.frame(data.chh.r2);
    
    list(cpg.r1=data.cpg.r1, cpg.r2=data.cpg.r2, chg.r1=data.chg.r1,
    chg.r2=data.chg.r2, chh.r1=data.chh.r1, chh.r2=data.chh.r2);

}


# Now read in the M-bias data for the sample
data.in <- read.mbias.file(paste(sample.id, "_bismark_M-bias.txt", sep=""));

# Now get out the different parts
data.cpg.r1 <- data.in[["cpg.r1"]][,c(2:3,5)];
data.cpg.r2 <- data.in[["cpg.r2"]][,c(2:3,5)];
data.chg.r1 <- data.in[["chg.r1"]][,c(2:3,5)];
data.chg.r2 <- data.in[["chg.r2"]][,c(2:3,5)];
data.chh.r1 <- data.in[["chh.r1"]][,c(2:3,5)];
data.chh.r2 <- data.in[["chh.r2"]][,c(2:3,5)];

# Now output the graphs
pdf(paste(sample.id,"Mbias_plots.pdf", sep="_"));

# CpG Context
data.cpg.r1$per.meth <- data.cpg.r1$count.methylated/data.cpg.r1$coverage;
data.cpg.r1$per.meth[is.na(data.cpg.r1$per.meth)] <- 0;

data.cpg.r2$per.meth <- data.cpg.r2$count.methylated/data.cpg.r2$coverage;
data.cpg.r2$per.meth[is.na(data.cpg.r2$per.meth)] <- 0;

plot(data.cpg.r1$per.meth,
type="l", col="blue", ylim=c(0,1),
xlab="Position", ylab="Methylation",
main=paste("M-bias CpG Context - ", sample.id, sep=""));

points(data.cpg.r2$per.meth, type="l", col="blue", lty="dashed")

# CHG Context
data.chg.r1$per.meth <- data.chg.r1$count.methylated/data.chg.r1$coverage;
data.chg.r1$per.meth[is.na(data.chg.r1$per.meth)] <- 0;

data.chg.r2$per.meth <- data.chg.r2$count.methylated/data.chg.r2$coverage;
data.chg.r2$per.meth[is.na(data.chg.r2$per.meth)] <- 0;

plot(data.chg.r1$per.meth,
type="l", col="blue", ylim=c(0,1),
xlab="Position", ylab="Methylation",
main=paste("M-bias CHG Context - ", sample.id, sep=""));

points(data.chg.r2$per.meth, type="l", col="blue", lty="dashed")

# CHH Context
data.chh.r1$per.meth <- data.chh.r1$count.methylated/data.chh.r1$coverage;
data.chh.r1$per.meth[is.na(data.chh.r1$per.meth)] <- 0;

data.chh.r2$per.meth <- data.chh.r2$count.methylated/data.chh.r2$coverage;
data.chh.r2$per.meth[is.na(data.chh.r2$per.meth)] <- 0;

plot(data.chh.r1$per.meth,
type="l", col="blue", ylim=c(0,1),
xlab="Position", ylab="Methylation",
main=paste("M-bias CHH Context - ", sample.id, sep=""));

points(data.chh.r2$per.meth, type="l", col="blue", lty="dashed")

dev.off();


q(save = "no");
