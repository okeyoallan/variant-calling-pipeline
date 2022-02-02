// Enable DSL 2 syntax

nextflow.enable.dsl = 2

// Process 1a Quality checking. Tool: fastqc

process QUALITY_CHECK {
    publishDir path: "${params.outdir}/Qc", mode: 'copy'
    tag "Quality Checking"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_logs"

    script:
    """
    mkdir ${sample_id}_logs
    fastqc -o ${sample_id}_logs -f fastq -q ${reads}
    """
}

// Process 1b Multi_QC raw reads. Tool: Multiqc

process MULTIQC {
        publishDir path: "${params.outdir}/MultiQc", mode: 'copy'
        tag "Raw Multiqc report"

        input:
        file(qualitycheck_out)

        output:
        file('multiqc_report.html')

        script:

        """
        multiqc .
        """
}

// Process 1c Trimming raw reads. Tool: Trimmomatic

process TRIMMOMATIC {
        publishDir path: "${params.outdir}/Trimming", mode: 'copy'
        tag "Trimming raw reads"

        input:
        tuple val(sample_id), path(reads)
        path adapter

        output:
        tuple path(fq_1_paired), path(fq_2_paired)

        script:
        fq_1_paired = sample_id + '_R1.paired.fastq'
        fq_1_unpaired = sample_id + '_R1.unpaired.fastq'
        fq_2_paired = sample_id + '_R2.paired.fastq'
        fq_2_unpaired = sample_id + '_R2.unpaired.fastq'

        """
        trimmomatic \
        PE -phred33 \
        ${reads[0]} \
        ${reads[1]} \
        $fq_1_paired \
        $fq_1_unpaired \
        $fq_2_paired \
        $fq_2_unpaired \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:${adapter}:2:40:15
        """
}

// Process 1d Post Fastqc. Tool: Fastqc

process POST_FASTQC {
        publishDir path: "${params.outdir}/PostFastQc", mode: 'copy'
        tag " Post fastQc"

        input:
        tuple path(read_R1), path(read_R2)

        output:
        path "${sample_id}_log"

        script:
        sample_id = ( read_R1 =~ /(.+)_R1.paired.fastq/ )[0][1]

        """
        mkdir ${sample_id}_log
        fastqc -o ${sample_id}_log -f fastq -q ${read_R1} ${read_R2}
        """
}

// Process 1e Post MultiQc. Tool: Multiqc

process MULTIQC_P {
        publishDir path: "${params.outdir}/PostMultiQC", mode: 'copy'
        tag "MultiQc report post trimming"

        input:
        file(postfastqc_out)

        output:
        file('multiqc_report.html')

        script:
        """
        multiqc .
        """
}

// Process 2 indexing ref_genome & alignment. Tool: bwa index/mem

process ALIGNMENT{
        publishDir path: "${params.outdir}/aligned", mode:'copy'
        tag "Alignment"

        input:
        path ref_ch
        tuple path(read_R1), path(read_R2)

        output:
        path "${sample_id}.sam", emit: aligned_sam

        script:
        sample_id = ( read_R1 =~ /(.+)_R1.paired.fastq/ )[0][1]

        template 'align.sh'
}

// Process 3 Merge samfiles. Tool: Samtools merge

process MERGE_SAM {
        publishDir path : "${params.outdir}/merged", mode: 'copy'
        tag "Merging Samfiles"

        input:
        path align_sam

        output:
        path "merged.sam", emit: merge_sam

        script:
        merged = "merged.sam"

        """
        samtools merge ${merged} ${align_sam}
        """
}

// Process 4  convert to bam format, sort and index. Tool Samtools

process CONVERT_TO_BAM {
      publishDir path : "${params.outdir}/conversion", mode: 'copy'
      tag "Conversion"

      input:
      path merge_sam

      output:
      path "merged.bam", emit: merge_bam
      path "merged_sorted.bam", emit: sort_bam
      path "merged.bam.bai", emit: index_bam

      script:
      bam = "merged.bam"
      sort = "merged_sorted.bam"
      index = "merged.bam.bai"

      """
      samtools view -Sb ${merge_sam} > ${bam}

      samtools sort -O bam -o ${sort} ${bam}

      samtools index ${sort} > ${index}
      """
}

// Process 5 marking and removing duplicates. Tool: gatk MarkDuplicates

process REMOVE_DUPLICATES {
      publishDir path: "${params.outdir}/Dedups", mode:'copy'
      tag "Removing Duplicates"

      input:
      path sort_bam

      output:
      path "marked_dups.bam", emit: marked_dups
      path "marked_dups_metrics.txt", emit: met_dedups

      script:
      dedups = "marked_dups.bam"
      metdedups = "marked_dups_metrics.txt"

      """


      gatk MarkDuplicates --INPUT ${sort_bam} --OUTPUT ${dedups} \
      --METRICS_FILE ${metdedups} --REMOVE_DUPLICATES true

      samtools index ${dedups}

      """
}

// Process 6 sequence dictionary. Tool: gatk CreateSequenceDictionary

process CREATE_SEQ_DICTIONARY {

        publishDir path: "${params.outdir}"
        tag " Creating Sequence Dictionary"

        input:
        path ref_chann

        output:
        path "*.dict", emit: genome_dict
        path "${fai}", emit: genome_fai



        script:
        fai = "${ref_chann}.fai"

        """
        gatk CreateSequenceDictionary -R ${ref_chann}

        samtools faidx ${ref_chann} > ${fai}
        """
}


// Process 7 variant call first round uncalibrated

process VARIANT_CALL_1 {
     publishDir path: "${params.outdir}/raw_variants", mode: 'copy'
     tag "Variant call_1"

     input:
     path marked_dup
     path ref_chh
     path genome_dict
     path genome_fai

     output:
     path "raw_variants.vcf", emit: raw_variants


     when:
     !params.knownsites

     script:
     raw_variants = "raw_variants.vcf"
     """
     gatk HaplotypeCaller -I ${marked_dup} -O ${raw_variants} -R ${ref_chh}
     """
}


// variant filter for the raw variants
process VARIANT_FILTER_1 {
      publishDir path: "${params.outdir}/filteredrawvariants", mode: 'copy'
      tag "Filtering raw variants"

      input:
      path variants_raw

      output:
      path "SRR_filteredb_raw_variants.vcf", emit: filteredrawvariants

      script:
      raw_var_filt = "SRR_filteredb_raw_variants.vcf"

      """
      gatk VariantFiltration -V ${variants_raw} -filter  "QD < 6.0" \
      --filter-name "QD6" -filter "MQ < 40.0" \
      --filter-name "MQ40" -filter "FS > 50.0" \
      --filter-name "FS50" -filter "SOR > 3.75" \
      --filter-name "SOR3.75" -filter "MQRankSum < -5.0" \
      --filter-name "MQRankSum5" -filter "ReadPosRankSum < -8.0" \
      --filter-name "ReadPosRankSum-8" -O variants_raw_tagged.vcf

      gatk SelectVariants -V variants_raw_tagged.vcf -O ${raw_var_filt} -select 'vc.isNotFiltered()'

      """
}

// Process Indexing and BaseRecalibration of the raw variants Tool: gatk

process BASERECALIBRATION_1{
        publishDir path: "${params.outdir}/Recal_1", mode: 'copy'
        tag "BaseRecalibration and Indexing known-features"

        input:
        path raw_filt
        path marked_dp
        path ref_ch
        path genome_dict
        path genome_fai

        output:
        path "recal.table", emit: recal_table_1
        path "recal.bam", emit: recal_bam_1

        script:
        Rtable_1 = "recal.table"
        Recal_1 = "recal.bam"

        """
        gatk IndexFeatureFile -F  ${raw_filt}

        gatk BaseRecalibrator -I ${marked_dp} --known-sites ${raw_filt} -R ${ref_ch} -O ${Rtable_1}

        gatk ApplyBQSR -R ${ref_ch} -I ${marked_dp} -bqsr ${Rtable_1} -O ${Recal_1}

        """
}


// calling variants a second round to generate the known sites for baserecalibration

process VARIANT_CALL_2 {
     publishDir path: "${params.outdir}/raw_2_variants", mode: 'copy'
     tag "Variant call_2"
     input:
     path recal_1_t
     path ref_chnn
     path genome_dict
     path genome_fai
     output:
     path "raw_2_variants.vcf", emit: raw_2_variants
     script:
     raw_2_variants = "raw_2_variants.vcf"

     """
     gatk HaplotypeCaller -I ${recal_1_t} -O ${raw_2_variants} -R ${ref_chnn}
     """
}



// filtering variants from the second variant call

process VARIANT_FILTER_2 {
      publishDir path: "${params.outdir}/filtered_2_rawvariants", mode: 'copy'
      tag "filtering raw variants_2"

      input:
      path variants_2_raw

      output:
      path "variants_2_raw.vcf", emit: filteredraw_2_variants

      script:
      raw_2_var_filt = "variants_2_raw.vcf"

      """
      gatk VariantFiltration -V ${variants_2_raw} -filter "QD < 6.0" \
      --filter-name "QD6" -filter "MQ < 40.0" \
      --filter-name "MQ40" -filter "FS > 50.0" \
      --filter-name "FS50" -filter "SOR > 3.75" \
      --filter-name "SOR3.75" -filter "MQRankSum < -5.0" \
      --filter-name "MQRankSum5" -filter "ReadPosRankSum < -8.0" \
      --filter-name "ReadPosRankSum-8" -O variants_2_raw.vcf

      gatk SelectVariants -V variants_2_raw.vcf -O ${raw_2_var_filt} -select 'vc.isNotFiltered()'

      """
}

process BASERECALIBRATION_2{
        publishDir path: "${params.outdir}/Recal", mode: 'copy'
        tag "BaseRecalibration and Indexing known-features"

        input:
        path knownv
        path marked_dp
        path ref_ch
        path genome_dict
        path genome_fai


        output:
        path "recal.table", emit: recal_table
        path "recal.bam", emit: recal_bam

        script:
        Rtable = "recal.table"
        Recal = "recal.bam"


        """
        gatk IndexFeatureFile -F  ${knownv}

        gatk BaseRecalibrator -I ${marked_dp} --known-sites ${knownv} -R ${ref_ch} -O ${Rtable}

        gatk ApplyBQSR -R ${ref_ch} -I ${marked_dp} -bqsr ${Rtable} -O ${Recal}

        """
}



// Process 7 Indexing and BaseRecalibration. Tool: gatk
process BASERECALIBRATION{
        publishDir path: "${params.outdir}/Recal", mode: 'copy'
        tag "BaseRecalibration and Indexing known-features"

        input:
        path knownch
        path marked_dp
        path ref_ch
        path genome_dict
        path genome_fai


        output:
        path "recal.table", emit: recal_table
        path "recal.bam", emit: recal_bam

        script:
        Rtable = "recal.table"
        Recal = "recal.bam"


        if (params.knownsites)

        """
        gatk IndexFeatureFile -F  ${knownch}

        gatk BaseRecalibrator -I ${marked_dp} --known-sites ${knownch} -R ${ref_ch} -O ${Rtable}

        gatk ApplyBQSR -R ${ref_ch} -I ${marked_dp} -bqsr ${Rtable} -O ${Recal}

        """

}

// Process 8 variant call. Tool: gatk HaplotypeCaller

process VARIANT_CALL {
     publishDir path: "${params.outdir}/variants", mode: 'copy'
     tag "Variant call"

     input:
     path recal_b
     path ref_c
     path genome_dict
     path genome_fai

     output:
     path "variantsGHC.vcf", emit: variants_ghc

     script:
     variantsG = "variantsGHC.vcf"

     """
     gatk HaplotypeCaller -I ${recal_b} -O ${variantsG} -R ${ref_c}

     """
}

// Process 9 filtering variants. Tool: GATK filter

process VARIANT_FILTER {
      publishDir path: "${params.outdir}/filteredvariants", mode: 'copy'
      tag "filtering variants"

      input:
      path variants

      output:
      path "SRR_filteredb_variants.vcf", emit: filteredvariants

      script:
      var_filt = "SRR_filteredb_variants.vcf"

      """
      gatk VariantFiltration -V ${variants} -filter "QD < 6.0" \
      --filter-name "QD6" -filter "MQ < 40.0" \
      --filter-name "MQ40" -filter "FS > 50.0" \
      --filter-name "FS50" -filter "SOR > 3.75" \
      --filter-name "SOR3.75" -filter "MQRankSum < -5.0" \
      --filter-name "MQRankSum5" -filter "ReadPosRankSum < -8.0" \
      --filter-name "ReadPosRankSum-8" -O variants_tagged.vcf

      gatk SelectVariants -V variants_tagged.vcf -O ${var_filt} -select 'vc.isNotFiltered()'


      """

}

// Process 10 Variant decomposing and normalizing. Tool: vt

process NORMALIZATION {
        publishDir path: "${params.outdir}/decomposed", mode: 'copy'
        tag "Decomposing variants"

        input:
        path filtvar
        path ref_gen
        path genome_dict
        path genome_fai

        output:
        path "SRR_normalized_var.vcf", emit: norm_var

        script:
        norm = "SRR_normalized_var.vcf"
        decomp = "SRR_decomposed_var.vcf"

        """
        vt decompose -s -o ${decomp} ${filtvar}

        vt normalize -r ${ref_gen} -o ${norm} ${decomp}

        """
}

// Process 11 variant annotation. Tool: snpEff

process ANNOTATION {
      publishDir path: "${params.outdir}/annotated", mode: 'copy'
      tag "Variant annotation"

      input:
      path normvar
      path snpeff_data

      output:
      path "SRR_anno_variants.vcf", emit: annotatevar
      path "snpEff_genes.txt", emit: snpEff_genes
      path "snpEff_summary.html", emit: snpEff_html

      script:
      annota_var = "SRR_anno_variants.vcf"
      if (!params.snpeff_data)

      template 'snpeffdb.sh'


     else

     """
     snpEff eff $params.snpeff_db ${normvar}  > ${annota_var}

     """


}
