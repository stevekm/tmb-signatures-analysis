#!/usr/bin/env Rscript
# Script for running deconstructSigs to produce genomic signatures
# need at least 55 variants per sample !! 
# https://cancer.sanger.ac.uk/cosmic/signatures

library("BSgenome.Hsapiens.UCSC.hg19")
library("deconstructSigs")
signature_type <- 'signatures.cosmic'
tri.counts.method <- 'exome'

args <- commandArgs(TRUE)
cohortID <- args[1]
output_Rdata <- args[2]
signatures_plot_pdf <- args[3]
signatures_pieplot_pdf <- args[4]
signatures_weights_tsv <- args[5]

input_vcfs <- args[6:length(args)]

cohort_label <- sprintf("Cohort.%s", cohortID)
vcf_colnames <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample")

vcf_df <- do.call('rbind', lapply(input_vcfs, function(x){
    df <- read.delim(file = x, header = FALSE, sep = '\t',
                     comment.char = '#', col.names = vcf_colnames, check.names = FALSE)
    return(df)
}))

# subset cols
vcf_df <- vcf_df[, c("CHROM", "POS", "REF", "ALT")]
vcf_df[["SampleID"]] <- cohort_label

# keep only entries with chroms in the reference data
genome_seqnames <- seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
vcf_df <- vcf_df[which(as.character(vcf_df[["CHROM"]]) %in% genome_seqnames), ]

# remove duplicates
vcf_df <- vcf_df[which(!duplicated(vcf_df)), ]

# convert to signatures format
sigs.input <- mut.to.sigs.input(mut.ref = vcf_df,
                                 sample.id = "SampleID",
                                 chr = "CHROM",
                                 pos = "POS",
                                 ref = "REF",
                                 alt = "ALT")

signatures <- whichSignatures(tumor.ref = sigs.input,
                               signatures.ref = signatures.cosmic, # signature_type
                               sample.id = cohort_label,
                               contexts.needed = TRUE,
                               tri.counts.method = tri.counts.method)

signature_weights <- signatures[["weights"]]
signature_weights[["SampleID"]] <- cohort_label
signature_weights[["SignatureType"]] <- signature_type
write.table(x = signature_weights, file = signatures_weights_tsv, sep = '\t', row.names = FALSE, col.names = TRUE)

pdf(file = signatures_plot_pdf)
print(plotSignatures(signatures, sub = signature_type))
dev.off()

pdf(file = signatures_pieplot_pdf)
print(makePie(signatures, sub = signature_type))
dev.off()


save.image(file = output_Rdata)