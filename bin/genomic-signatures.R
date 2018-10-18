#!/usr/bin/env Rscript
# Script for running deconstructSigs to produce genomic signatures
# https://github.com/raerose01/deconstructSigs
# need at least 55 variants per sample !! 
# https://cancer.sanger.ac.uk/cosmic/signatures

library("BSgenome.Hsapiens.UCSC.hg19")
library("deconstructSigs")
signature_type <- 'signatures.cosmic'
tri.counts.method <- 'exome'

args <- commandArgs(TRUE)
sampleID <- args[1]
input_vcf <- args[2]
output_Rdata <- args[3]
signatures_Rds <- args[4]
signatures_plot_Rds <- args[5]
signatures_pie_plot_Rds <- args[6]
signatures_plot_pdf <- args[7]
signatures_pie_plot_pdf <- args[8]
signatures_weights_tsv <- args[9]
signatures_input_Rds <- args[10]

# load .vcf
vcf_colnames <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sampleID)
variants <- read.delim(file = input_vcf, header = FALSE, sep = '\t',
                       comment.char = '#', col.names = vcf_colnames, check.names = FALSE)

# add sample ID column
variants[["SampleID"]] <- sampleID

# keep only entries with chroms in the reference data
genome_seqnames <- seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
variants <- variants[which(as.character(variants[["CHROM"]]) %in% genome_seqnames), ]

# convert to signatures format
sigs.input <- mut.to.sigs.input(mut.ref = variants,
                                sample.id = "SampleID",
                                chr = "CHROM",
                                pos = "POS",
                                ref = "REF",
                                alt = "ALT")

# save signatures input
saveRDS(object = sigs.input, file = signatures_input_Rds, compress = TRUE)

# make the signatures
signatures <- whichSignatures(tumor.ref = sigs.input,
                              signatures.ref = signatures.cosmic, # signature_type
                              sample.id = sampleID,
                              contexts.needed = TRUE,
                              tri.counts.method = tri.counts.method)

signature_weights <- signatures[["weights"]]
signature_weights[["SampleID"]] <- sampleID
signature_weights[["SignatureType"]] <- signature_type
write.table(x = signature_weights, file = signatures_weights_tsv, sep = '\t', row.names = FALSE, col.names = TRUE)


# save signatures
saveRDS(object = signatures, file = signatures_Rds, compress = TRUE)


# make plots
# https://stackoverflow.com/a/29583945/5359531
pdf(file = signatures_plot_pdf)
dev.control(displaylist="enable") 
print(plotSignatures(signatures, sub = signature_type))
saveRDS(object = recordPlot(), file = signatures_plot_Rds)
dev.off()

pdf(file = signatures_pie_plot_pdf)
dev.control(displaylist="enable") 
print(makePie(signatures, sub = signature_type))
saveRDS(object = recordPlot(), file = signatures_pie_plot_Rds)
dev.off()

save.image(file = output_Rdata)