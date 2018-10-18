params.outputDir = "output"
params.inputDir = "samples"

params.ANNOVAR_BUILD_VERSION ="hg19"
params.ANNOVAR_PROTOCOL = "refGene,cosmic70,dbnsfp33a,fathmm,1000g2015aug_all,exac03,snp138,snp129"
params.ANNOVAR_OPERATION ="g,f,f,f,f,f,f,f"

def outputDir = "${params.outputDir}"

Channel.fromPath('targets.bed').into { targets_bed1; targets_bed2; targets_bed3; targets_bed4 }
Channel.fromPath("${params.refDir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa").into { genome_fa; genome_fa2 }
Channel.fromPath("${params.refDir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai").into { genome_fai; genome_fai2 }
Channel.fromPath("${params.refDir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict").into { genome_dict; genome_dict2 }

Channel.fromPath("${params.ANNOVAR_DB_DIR}").set { annovar_db_dir }

Channel.fromFilePairs("${params.inputDir}/**{.dd.ra.rc.bam,.dd.ra.rc.bam.bai}").map{ sampleID, items ->
    def bamfile = items[0]
    def baifile = items[1]
    return([sampleID, bamfile, baifile])
}
.set { sample_bams }

Channel.fromPath("${params.inputDir}/**.HaplotypeCaller.filtered.vcf")
.map { path ->
    def dirpath = new File("${path}").parent
    def dirname = new File("${dirpath}").name
    def sampleID = "${dirname}"
    def caller = "HaplotypeCaller"
    return([sampleID, caller, path])
}.set { hc_vcfs }

Channel.fromPath("${params.inputDir}/**.LoFreq.filtered.vcf")
.map { path ->
    def dirpath = new File("${path}").parent
    def dirname = new File("${dirpath}").name
    def sampleID = "${dirname}"
    def caller = "LoFreq"
    return([sampleID, caller, path])
}.set { lofreq_vcfs }


process merge_targets {
    publishDir "${outputDir}", mode: "copy", overwrite: true

    input:
    file(targets) from targets_bed2

    output:
    file("${merged}") into (targets_merged, targets_merged2, targets_merged3)

    script:
    merged = "merged_targets.bed"
    """
    # merge all regions, regardless of strand
    bedtools merge -i "${targets}" > "${merged}"
    """
}


process count_regions {
    publishDir "${outputDir}", mode: "copy", overwrite: true

    input:
    file(bed) from targets_bed1.concat(targets_merged)

    output:
    file("${output_file}") into region_counts

    script:
    output_file = "${bed}.region-count.txt"
    """
    wc -l "${bed}" > "${output_file}"
    """
}


process count_distance {
    publishDir "${outputDir}", mode: "copy", overwrite: true

    input:
    file(bed) from targets_bed3.concat(targets_merged2)

    output:
    file("${output_file}") into bed_distances

    script:
    output_file = "${bed}.distance.txt"
    """
    bed-size.py "${bed}" > "${output_file}"
    """
}

sample_bams.combine(genome_fa)
.combine(genome_fai)
.combine(genome_dict)
.combine(targets_bed4)
.tap { sample_inputs_ref }

process calculate_callable_loci {
    tag "${prefix}"
    publishDir "${outputDir}/${sampleID}", mode: "copy", overwrite: true

    input:
    set val(sampleID), file(bamfile), file(baifile), file(genomeFa), file(genomeFai), file(genomeDict), file(targetsBed) from sample_inputs_ref

    output:
    set val("${sampleID}"), file("${output_summary}"), file("${output_bed_pass}") into called_loci
    file("${output_summary}")
    file("${output_bed}")
    file("${output_bed_pass}")


    script:
    prefix = "${sampleID}"
    output_summary = "${prefix}.CallableLoci.summary.txt"
    output_bed = "${prefix}.CallableLoci.bed"
    output_bed_pass = "${prefix}.CallableLoci.pass.bed"
    output_bed_pass_size = "${prefix}.CallableLoci.pass.distance.txt"
    """
    # minDepth 500 = NGS580 call threshold
    # minMappingQuality, minBaseQuality 20 = NGS580 DepthOfCoverage threshold

    gatk.sh \
    -T CallableLoci \
    -R "${genomeFa}" \
    -I "${bamfile}" \
    -summary "${output_summary}" \
    --minMappingQuality 20 \
    --minBaseQuality 20 \
    --minDepth 500 \
    --intervals "${targetsBed}" \
    -o "${output_bed}"

    grep -E 'CALLABLE|PASS' "${output_bed}" > "${output_bed_pass}" || touch "${output_bed_pass}" # exit code 1 if no matches

    # bed-size.py "${output_bed_pass}" > "${output_bed_pass_size}" # not in the GATK3.8 .simg
    """
}

process callable_loci_table {
    tag "${prefix}"
    publishDir "${outputDir}/${sampleID}", mode: "copy", overwrite: true

    input:
    set val(sampleID), file(summary), file(bed) from called_loci

    output:
    file("${output_summary}") into loci_tables
    set val(sampleID), file("${output_summary}"), file("${output_txt}") into loci_tables2

    script:
    prefix = "${sampleID}"
    output_summary = "${prefix}.CallableLoci.summary.tsv"
    output_txt = "${prefix}.CallableLoci.txt"
    tmpFile = "tmp"
    """
    callable-loci-table.py "${summary}" "${tmpFile}"
    paste-col.py -i "${tmpFile}" --header "Sample" -v "${sampleID}" > "${output_summary}"
    grep 'CALLABLE' "${tmpFile}" | cut -f2 > "${output_txt}"
    """
}

process collect_loci_tables {
    publishDir "${outputDir}", overwrite: true, mode: 'copy'

    input:
    file('t*') from loci_tables.collect()

    output:
    file("${output_table}")

    script:
    output_table = "all_callable_loci.tsv"
    """
    concat-tables.py * > "${output_table}"
    """
}

process filter_vcf{
    tag "${prefix}"
    publishDir "${outputDir}/${sampleID}", mode: "copy", overwrite: true

    input:
    set val(sampleID), val(caller), file(vcf), file(genomeFa), file(genomeFai), file(genomeDict) from lofreq_vcfs.concat(hc_vcfs).combine(genome_fa2).combine(genome_fai2).combine(genome_dict2)

    output:
    set val(sampleID), val(caller), file("${output_filtered_vcf}"), file("${output_filtered_tsv}") into filtered_vcfs_tsvs
    set val(sampleID), val(caller), file("${output_filtered_vcf}") into filtered_vcfs
    file("${output_tsv}")

    script:
    prefix = "${sampleID}.${caller}"
    output_filtered_vcf = "${prefix}.variants.re-filtered.vcf"
    output_filtered_tsv = "${prefix}.variants.re-filtered.tsv"
    output_tsv = "${prefix}.variants.tsv"
    if( caller == 'HaplotypeCaller' )
        """
        # report if
        # alternate allele freq (allele depth / depth) greater than 0.5
        # more than 5 variant call supporting reads
        # quality reads present (reported depth >0)
        gatk.sh -T SelectVariants \
        -R "${genomeFa}" \
        -V "${vcf}" \
        --sample_name "${sampleID}" \
        -select "vc.getGenotype('${sampleID}').getAD().1 / vc.getGenotype('${sampleID}').getDP() > 0.50" \
        -select "vc.getGenotype('${sampleID}').getAD().1 > 5" \
        -select "vc.getGenotype('${sampleID}').getDP() > 0" \
        > "${output_filtered_vcf}"

        gatk.sh -T VariantsToTable \
        -R "${genomeFa}" \
        -V "${output_filtered_vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN \
        -GF AD -GF DP \
        -o "${output_filtered_tsv}"

        gatk.sh -T VariantsToTable \
        -R "${genomeFa}" \
        -V "${vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN \
        -GF AD -GF DP \
        -o "${output_tsv}"
        """
    else if( caller == 'LoFreq' )
        """
        # do not report if frequency is less than 5%
        gatk.sh -T SelectVariants \
        -R "${genomeFa}" \
        -V "${vcf}" \
        -select "AF > 0.05"  \
        > "${output_filtered_vcf}"

        gatk.sh -T VariantsToTable \
        -R "${genomeFa}" \
        -V "${output_filtered_vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F DP -F AF -F SB -F INDEL -F CONSVAR -F HRUN \
        -o "${output_filtered_tsv}"

        gatk.sh -T VariantsToTable \
        -R "${genomeFa}" \
        -V "${vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F DP -F AF -F SB -F INDEL -F CONSVAR -F HRUN \
        -o "${output_tsv}"
        """
    else
        error "Invalid caller: ${caller}"
}

def deconstructSigs_variant_min = 55 // min num variants for signatures library
filtered_vcfs.tap { filtered_vcfs2 }
filtered_vcfs2.filter{ sampleID, caller, vcf ->
    // make sure there are enough variants in the VCF to proceed!
    // need at least ~55 variants
    def line_count = 0
    def num_variants = 0
    def enough_variants = false // bad by default

    // count number of variants in the vcf file
    vcf.withReader { reader ->
        while (line = reader.readLine()) {
            if (!line.startsWith("#")) num_variants++
            if (num_variants > deconstructSigs_variant_min) {
                enough_variants = true
                break
                }
            line_count++
        }
    }
    return(enough_variants)
}
.tap { filtered_vcfs_numchecked; filtered_vcfs_numchecked2 }

filtered_vcfs_numchecked2.map{ sampleID, caller, vcf ->
    return([caller, vcf])
}.groupTuple().set{ filtered_vcfs_numchecked_groups }

process genomic_signatures {
    tag "${prefix}"
    publishDir "${outputDir}/${sampleID}", mode: "copy", overwrite: true

    input:
    set val(sampleID), val(caller), file(vcf) from filtered_vcfs_numchecked

    output:
    file("${output_Rdata}")
    file("${signatures_Rds}")
    file("${signatures_plot_Rds}")
    file("${signatures_pieplot_Rds}")
    file("${signatures_plot_pdf}")
    file("${signatures_pieplot_pdf}")
    file("${signatures_weights_tsv}") into genomic_signature_weights
    set file("${signatures_plot_pdf}"), file("${signatures_pieplot_pdf}") into genomic_signature_plots
    // file("${signatures_input_Rds}") into genomic_signatures_inputs

    script:
    prefix = "${sampleID}.${caller}"
    output_Rdata = "${prefix}.signatures.Rdata"
    signatures_Rds = "${prefix}.signatures.Rds"
    signatures_plot_Rds = "${prefix}.signatures.plot.Rds"
    signatures_pieplot_Rds = "${prefix}.signatures.pieplot.Rds"
    signatures_plot_pdf = "${prefix}.signatures.plot.pdf"
    signatures_pieplot_pdf = "${prefix}.signatures.pieplot.pdf"
    signatures_weights_tsv = "${prefix}.signatures.weights.tsv"
    signatures_input_Rds = "${prefix}.signatures.input.Rds"
    """
    genomic-signatures.R "${prefix}" "${vcf}" "${output_Rdata}" "${signatures_Rds}" "${signatures_plot_Rds}" "${signatures_pieplot_Rds}" "${signatures_plot_pdf}" "${signatures_pieplot_pdf}" "${signatures_weights_tsv}" "${signatures_input_Rds}"
    """
}


genomic_signature_plots.flatten().toSortedList({ a, b ->
    // sort by file basename
    a.name <=> b.name
    }).set{ sorted_genomic_signature_plots }

process collect_signature_plots {
    publishDir "${outputDir}", mode: "copy", overwrite: true

    input:
    file("p*.pdf") from sorted_genomic_signature_plots

    output:
    file("${output_pdf}") into all_samples_signature_plot

    script:
    output_pdf = "all_genomic_signatures.pdf"
    """
    # 'natural sort' numeric order...
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None -sOutputFile="${output_pdf}" \$(ls -v *.pdf)
    """
}

process genomic_signatures_cohort {
    tag "${prefix}"
    publishDir "${outputDir}", mode: "copy", overwrite: true

    input:
    set val(caller), file('*.vcf') from filtered_vcfs_numchecked_groups

    output:
    file("${signatures_weights_tsv}") into genomic_signature_cohort_weights
    set file("${signatures_plot_pdf}"), file("${signatures_pieplot_pdf}") into genomic_signatures_cohort_plots
    file("${output_Rdata}")

    script:
    prefix = "${caller}"
    output_Rdata = "${prefix}.all.signatures.Rdata"
    signatures_plot_pdf = "${prefix}.signatures.plot.pdf"
    signatures_pieplot_pdf = "${prefix}.signatures.pieplot.pdf"
    signatures_weights_tsv = "${prefix}.signatures.weights.tsv"
    """
    genomic-signatures-cohort.R "${prefix}" "${output_Rdata}" "${signatures_plot_pdf}" "${signatures_pieplot_pdf}" "${signatures_weights_tsv}" *
    """
}
genomic_signatures_cohort_plots.flatten().set { genomic_signatures_cohort_plots_flattened }

process collect_signature_weights {
    publishDir "${outputDir}", mode: "copy", overwrite: true

    input:
    file('t*') from genomic_signature_weights.concat(genomic_signature_cohort_weights).collect()

    output:
    file("${output_table}") into all_signatures_weights

    script:
    output_table = "all_genomic_signature_weights.tsv"
    """
    concat-tables.py * > "${output_table}"
    """
}


process annotate_vcf {
    tag "${prefix}"
    publishDir "${outputDir}/${sampleID}", mode: "copy", overwrite: true

    input:
    set val(sampleID), val(caller), file(vcf), file(tsv), file(annovar_db) from filtered_vcfs_tsvs.combine(annovar_db_dir)

    output:
    set val(sampleID), val(caller), file("${avinput_file}"), file("${annovar_output_txt}") into (annotations, annotations2)

    script:
    prefix = "${sampleID}.${caller}"
    avinput_file = "${prefix}.avinput"
    annovar_output_txt = "${prefix}.${params.ANNOVAR_BUILD_VERSION}_multianno.txt"
    // annovar_output_vcf = "${prefix}.${params.ANNOVAR_BUILD_VERSION}_multianno.vcf"
    """
    convert2annovar.pl \
    -includeinfo \
    -format vcf4 \
    "${vcf}" > \
    "${avinput_file}"

    table_annovar.pl "${avinput_file}" "${annovar_db}" \
    --buildver "${params.ANNOVAR_BUILD_VERSION}" \
    --remove \
    --protocol "${params.ANNOVAR_PROTOCOL}" \
    --operation "${params.ANNOVAR_OPERATION}" \
    --nastring . \
    --onetranscript \
    --outfile "${prefix}"
    """
}

process update_annotation_tables {
    tag "${prefix}"
    publishDir "${outputDir}/${sampleID}", mode: "copy", overwrite: true

    input:
    set val(sampleID), val(caller), file(avinput), file(annovar_txt) from annotations2

    output:
    file("${output_table}") into updated_annotation_tables

    script:
    prefix = "${sampleID}.${caller}"
    output_table = "${prefix}.annotations.tsv"
    """
    cat "${annovar_txt}" | \
    paste-col.py --header "Sample" -v "${sampleID}" | \
    paste-col.py --header "Caller" -v "${caller}" > "${output_table}"
    """
}

process collect_annotation_tables {
    publishDir "${outputDir}", overwrite: true, mode: 'copy'

    input:
    file('t*') from updated_annotation_tables.collect()

    output:
    file("${output_table}") into all_annotations

    script:
    output_table = "all_annotations.tsv"
    """
    concat-tables.py --quoteOutput * > "${output_table}"
    """
}

// only keep files with at least 1 variant
annotations.filter { sampleID_anno, caller, avinput, anno_txt ->
    def count = anno_txt.readLines().size()
    count > 1
}.set { annotations_filtered }

process tmb_filter_variants {
    tag "${prefix}"
    publishDir "${outputDir}/${sampleID}", mode: "copy", overwrite: true

    input:
    set val(sampleID), val(caller), file(avinput), file(anno_txt) from annotations_filtered

    output:
    file("${output_Rdata}")
    file("${output_variants}") into filtered_variants
    set val(sampleID), val(caller), file("${output_variants}") into filtered_variants2

    script:
    prefix = "${sampleID}.${caller}"
    output_Rdata = "${prefix}.Rdata"
    output_variants = "${prefix}.annotations.filtered.tsv"
    tmpFile = "tmp"
    """
    tmb-variant-filter.R "${output_Rdata}" "${tmpFile}" "${anno_txt}"
    cat "${tmpFile}" | \
    paste-col.py --header "Sample" -v "${sampleID}" | \
    paste-col.py --header "Caller" -v "${caller}" > "${output_variants}"
    """
}

process collect_filtered_annotation_tables {
    publishDir "${outputDir}", overwrite: true, mode: 'copy'

    input:
    file('t*') from filtered_variants.collect()

    output:
    file("${output_table}") into (output_filtered_annotations, output_filtered_annotations2)

    script:
    output_table = "all_filtered_annotations.tsv"
    """
    concat-tables.py --quoteOutput * > "${output_table}"
    """
}



loci_tables2.map{ sampleID, summary_table, callable_txt ->
    def callable_loci = callable_txt.readLines()[0]
    return([sampleID, summary_table, "${callable_loci}"])
}.combine(filtered_variants2)
.filter{ sampleID_loci, loci_summary_table, callable_loci, sampleID_anno, caller, filtered_annotations ->
    sampleID_anno == sampleID_loci
}
.map{ sampleID_loci, loci_summary_table, callable_loci, sampleID_anno, caller, filtered_annotations ->
    def num_variants = filtered_annotations.readLines().size() - 1
    return([sampleID_loci, caller, callable_loci, num_variants])
}
.set { loci_annotations }

process calculate_tmb {
    tag "${prefix}"
    publishDir "${outputDir}/${sampleID}", mode: "copy", overwrite: true

    input:
    set val(sampleID), val(caller), val(loci), val(variants) from loci_annotations

    output:
    file("${output_tmb}") into tmbs

    script:
    prefix = "${sampleID}.${caller}"
    output_tmb = "${prefix}.tmb.tsv"
    """
    # tmb="\$(( ${variants} / ( ${loci} / 1000000 ) ))"
    tmb=\$( python -c 'print(float(${variants}) / float(${loci}) * 1000000 )' )
    printf 'SampleID\tCaller\tnBases\tnVariants\tTMB\n' > "${output_tmb}"
    printf "${sampleID}\t${caller}\t${loci}\t${variants}\t\${tmb}\n" >> "${output_tmb}"
    """
}
tmbs.collectFile(name: 'all_tmbs.tsv', keepHeader: true, storeDir: "${outputDir}").into { all_tmb; all_tmb2 }

process tmb_samples_comparison {
    publishDir "${outputDir}", overwrite: true, mode: 'copy'

    input:
    set file(tmbs), file(variants) from all_tmb2.combine(output_filtered_annotations2)

    output:
    file("${output_Rdata}")

    script:
    output_Rdata = "tmb-samples-compare.Rdata"
    """
    tmb-samples-comparison.R "${output_Rdata}" "${tmbs}" "${variants}"
    """
}

process collect_output {
    publishDir "${outputDir}", overwrite: true, mode: 'copy'

    input:
    file("*") from all_tmb.concat(output_filtered_annotations, targets_merged3, region_counts, bed_distances, all_annotations, all_samples_signature_plot, all_signatures_weights, genomic_signatures_cohort_plots_flattened).collect()

    output:
    file("${output_zip}")

    script:
    output_zip = "tmb_analysis_results.zip"
    """
    zip "${output_zip}" *
    """
}
