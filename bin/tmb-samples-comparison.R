#!/usr/bin/env Rscript
library("ggplot2")

args <- commandArgs(TRUE)
output_Rdata <- args[1]
tmb_tsv <- args[2]
anno_tsv <- args[3] 

# values used in the ANNOVAR output for nastring
NA_vals <- c('.')  

tmb_df_og <- read.delim(file = tmb_tsv, header = TRUE, sep = '\t')
anno_df <- read.delim(file = anno_tsv, header = TRUE, sep = '\t', na.strings = NA_vals)

# keep only LoFreq entries for this analysis
tmb_df <- tmb_df_og[which(tmb_df_og[["Caller"]] == "LoFreq"), ]


# figure out which samples have the following mutations (expected to be mutually exclusive)
# BRAF, NRAS, NF1
control_samples <- c("SC-SERACARE", "NC-HAPMAP")
all_samples <- levels(tmb_df[["SampleID"]])[which(! levels(tmb_df[["SampleID"]]) %in% control_samples)]

variant_groups <- list(
    BRAF = unique(as.character(anno_df[which(anno_df[["Gene.refGene"]] == "BRAF" & anno_df[["Sample"]] %in% all_samples ), ][["Sample"]])),
    NRAS = unique(as.character(anno_df[which(anno_df[["Gene.refGene"]] == "NRAS" & anno_df[["Sample"]] %in% all_samples ), ][["Sample"]])),
    NF1 = unique(as.character(anno_df[which(anno_df[["Gene.refGene"]] == "NF1" & anno_df[["Sample"]] %in% all_samples ), ][["Sample"]])),
    other = unique(as.character(anno_df[which(! anno_df[["Gene.refGene"]] %in% c("NF1", "NRAS", "BRAF") & anno_df[["Sample"]] %in% all_samples ), ][["Sample"]]))
)


# ggplot(data = tmb_df, aes(x = SampleID, y = TMB)) +
#     geom_point()

# hold off on this until later...
save.image(output_Rdata)