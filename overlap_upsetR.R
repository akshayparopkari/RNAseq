#########################################################################################
# DESCRIPTION
# Using output from DESeq2, run and plot UpSet diagram showing significantly and
# differentially expressed genes in multiple conditions.
#
# AUTHOR
# Akshay Paropkari
#
# VERSION
# 0.1.0
#########################################################################################

# library imports
library(reshape2, warn.conflicts = FALSE, quietly = T, verbose = F)
library(tools, warn.conflicts = FALSE, quietly = T, verbose = F)
library(UpSetR, warn.conflicts = FALSE, quietly = T, verbose = F)

# Get all arguments
args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please specify 'up' or 'down' for significantly and differentially up-regulated "
       "or down-regulated gene overlap plot, respectively.")
} else {
  if (tolower(args[1]) == "up"){
    expt <- "upregulated"
  }
  if (tolower(args[1]) == "down"){
    expt <- "downregulated"
  }
}

# Note that all deseq2 tab-separated output files must end with
# "[mutant]_deseq2_output.tsv" and be available inside <INPUT DIRECTORY>
# E.g. for Can1 mutant, filename can by can1_deseq2_output.tsv

input.folder <- <INPUT_DIRECTORY> # CHANGE THIS TO YOUR INPUT FOLDER LOCATION
setwd(input.folder)
deseq.res.files <- list.files(path = ".", pattern = "*_deseq2_output.tsv",
                              full.names = T)

# initiate list to save all significant ORFs
final.list <- list()

# iterate through all files in input folder, and calculate significantly
# upregulated ORFs
for (f in deseq.res.files) {
  message(paste0("PROCESSING ", f))
  if (length(unlist(strsplit(x = f, split = "/"))) == 2) {
    mutant <- unlist(strsplit(x = f, split = "/"))[2]
  } else {
    mutant <- unlist(strsplit(x = f, split = "/"))[1]
  }
  mutant <- unlist(strsplit(mutant, split = "_"))[1]
  message(paste0("FORMATTING ", toTitleCase(mutant),"-KO DATA SET"))
  deseq.data <- read.csv(file = f, sep = "\t", header = T, row.names = 1)

  # consider entries with ORF ids only
  deseq.data <- deseq.data[grep(pattern = "^orf*", x = rownames(deseq.data)), ]

  # remove ORFs with NA values in any column
  deseq.data <- deseq.data[complete.cases(deseq.data), ]

  # collect significant upregulated genes in mutant

  # significant and upregulated ORFs
  if (tolower(args[1]) == "up"){
    message(paste0("COLLECTING SIGNIFICANTLY UPREGULATED ORFs IN ",
                 toTitleCase(mutant),"-KO DATA SET"))
    orf.vec <- row.names(deseq.data)[which(deseq.data$padj < 0.05 & deseq.data$log2FoldChange > 0)]
    }

  # significant and downregulated ORFs
  if tolower(args[1]) == "down"{
    message(paste0("COLLECTING SIGNIFICANTLY DOWNREGULATED ORFs IN ", toTitleCase(mutant),
                 "-KO DATA SET"))
    orf.vec <- row.names(deseq.data)[which(deseq.data$padj < 0.05 & deseq.data$log2FoldChange < 0)]
    }

  # significant ORFs
#  if tolower(args[1] == NULL){
#    message(paste0("COLLECTING SIGNIFICANT AND DIFFERENTIALLY EXPRESSED ORFs IN ",
#                   toTitleCase(mutant),"-KO DATA SET"))
#    expr <- "differentially_expressed"
#    orf.vec <- row.names(deseq.data)[which(deseq.data$padj < 0.05)]
#  }

  # save ORF list into a list
  message(paste0("SAVING ORF LIST FOR ", toTitleCase(mutant), "-KO DATA SET"))
  name <- paste0(toTitleCase(mutant), "-KO")
  final.list[[name]] <- orf.vec
  message(paste0("=========================================================="))
}

# save ORF list to file
fnh <- paste0(input.folder, "/", expr, "_orf_list.tsv")
message(paste0("SAVING ORF LIST DATA TO FILE AT ", fnh))
final.list.df <- melt(data = final.list, value.name = "orf.ids")
final.list.df <- dcast(data = final.list.df, formula = orf.ids~L1,
                       value.var="orf.ids")
final.list.df$orf.ids <- NULL
final.list.df$row_sum <- apply(X = final.list.df, MARGIN = 1,
                                    FUN = function(x) length(which(!is.na(x))))
write.table(x = final.list.df, file = fnh, quote = F, sep = "\t", row.names = F)

# plot upset plot and save as  JPEG file on the Desktop
# NO SPACES in file name or output folder
image.fnh <- paste0(input.folder, "/", expr, "_orf_upset_plot.jpeg")
message(paste0("SAVING UPSET PLOT TO FILE AT ", image.fnh))
jpeg(filename = image.fnh, width = 21, height = 9, units = "in", res=300)
upset(data = fromList(final.list), nsets = 100, nintersects = NA,
      order.by = "freq", point.size = 4,
      mainbar.y.label = "Intersection set size",
      sets.x.label = "Set size",
      text.scale = c(2.5, 0, 2, 1.3, 2, 1.75))
dev.off()
