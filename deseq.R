#! /usr/bin/env Rscript

###############################################################################
#
# DESCRIPTION
#
# Abstract: This R script takes in gene counts table (output from python script
#           format_counts_table.py) and performs differential expression
#           analysis using DESeq2. The outputs are a results table with log2
#           fold changes, p values and adjusted p values and MA plot as tab-
#           separated and PDF file, respectively.
#
# Author: Akshay Paropkari
#
# -----------------------------------------------------------------------------
#
# USAGE and EXAMPLE
#
# Rscript --vanilla raw_gene_counts_file metadata_file results_file output_MA_file
# Rscript --vanilla gene_counts.txt metadata.txt deseq2_lfc.txt MA_plot.pdf
#
#-----------------------------------------------------------------------------
#
# DESeq2: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#
###############################################################################

library("DESeq2", warn.conflicts = FALSE, quietly = T, verbose = F)
library("readxl", warn.conflicts = FALSE, quietly = T, verbose = F)
library("IHW", warn.conflicts = FALSE, quietly = T, verbose = F)

# Get all arguments
args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied 'gene_raw_counts.txt'")
}

# Get metadata as input
coldata <- as.data.frame(read_excel(args[2]))
row.names(coldata) <- coldata$SampleID
print(rownames(coldata))

# Read in counts data
raw.counts <- as.matrix(read.delim(file = args[1], row.names = 1))

# Reorder samples (columns) to match samples (rows) in metadata file
raw.counts <- raw.counts[, rownames(coldata)]

# MUST BE TRUE
all(rownames(coldata) %in% colnames(raw.counts))
all(rownames(coldata) == colnames(raw.counts))

# Create a DGEList object
dds <- DESeqDataSetFromMatrix(countData = raw.counts, colData = coldata,
                              design = ~ Condition)
dds$Condition <- factor(dds$Condition, levels = c("WT", "Mutant"))
print(resultsNames(dds))

# Prefiltering out genes that don't have count of at least 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run differential expression and get helpful analysis messages
dds <- DESeq(dds)
resIHW <- results(dds, filterFun=ihw, alpha = 0.05)
print(summary(resIHW))
print(sum(resIHW$padj < 0.05, na.rm=TRUE))
print(metadata(resIHW)$ihwResult)

# Obtain  MA plot
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
write.table(as.data.frame(resLFC[order(resLFC$pvalue),]), file=args[3],
            quote = F, sep = "\t")
pdf(file = args[4])
plotMA(resLFC, ylim=c(-2,2))
abline(h = c(-1, 1), col="dodgerblue", lwd=2)
dev.off()
