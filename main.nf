// Enable DSL 2 syntax

nextflow.enable.dsl = 2



println """
===========================================
VARIANT CALLING PIPELINE USING GATK4 BEST PRACTICES
============================================
"""


// import modules here

include { QUALITY_CHECK; MULTIQC; TRIMMOMATIC; POST_FASTQC; MULTIQC_P; ALIGNMENT; MERGE_SAM; CONVERT_TO_BAM; REMOVE_DUPLICATES; CREATE_SEQ_DICTIONARY; VARIANT_CALL_1; VARIANT_FILTER_1; BASERECALIBRATION_1; VARIANT_CALL_2; VARIANT_FILTER_2; BASERECALIBRATION_2; BASERECALIBRATION; VARIANT_CALL; VARIANT_FILTER; NORMALIZATION; ANNOTATION } from "./Modules/gatkHC.nf"

// set input channels
Channel.fromFilePairs( params.total_reads, checkExists:true )
        .set { read_pairs_ch }

Channel.fromPath ( params.reference, checkIfExists:true )
        .set { reference_ch }

if (params.knownsites)
Channel.fromPath (params.knownsites)
        .set { known_ch }

Channel.fromPath ( params.adapter, checkIfExists:true )
       .set { adapter_ch }

Channel.fromPath ( params.snpeff_data, checkIfExists:true )
       .set { snpeff_ch }



// Run the workflow
workflow {
// step 1a Quality Checking
QUALITY_CHECK(read_pairs_ch)

// step 1b MultiQc_raw reads
MULTIQC(QUALITY_CHECK.out.collect())

// step 1c Trimming reads
TRIMMOMATIC(read_pairs_ch, adapter_ch.collect())

// step 1d Post trimming fastqc
POST_FASTQC(TRIMMOMATIC.out)

// step 1e Post MultiQC
MULTIQC_P(POST_FASTQC.out.collect())

// step 2 Alignment
ALIGNMENT(reference_ch.collect(), TRIMMOMATIC.out)

// step 3 Merging aligned sam files
MERGE_SAM(ALIGNMENT.out.collect())

// step 4 Bam Conversion
CONVERT_TO_BAM(MERGE_SAM.out)

// step 5 remove duplicates
REMOVE_DUPLICATES(CONVERT_TO_BAM.out.sort_bam)

// step 6 creating sequence dictionary
CREATE_SEQ_DICTIONARY(reference_ch)

// step first round uncalibrated to generate knownsites
VARIANT_CALL_1( REMOVE_DUPLICATES.out.marked_dups, reference_ch, CREATE_SEQ_DICTIONARY.out)

// filtering raw variants
VARIANT_FILTER_1(VARIANT_CALL_1.out)

// Baserecalibration first round for filtered raw variants
BASERECALIBRATION_1(VARIANT_FILTER_1.out, REMOVE_DUPLICATES.out.marked_dups, reference_ch, CREATE_SEQ_DICTIONARY.out)

// Second round of variant calling to generate knownsites vcf
VARIANT_CALL_2( BASERECALIBRATION_1.out.recal_bam_1, reference_ch, CREATE_SEQ_DICTIONARY.out)

// filtering second round of raw variants
VARIANT_FILTER_2(VARIANT_CALL_2.out)

// step 7 Baserecalibration
BASERECALIBRATION_2(VARIANT_FILTER_2.out, REMOVE_DUPLICATES.out.marked_dups, reference_ch, CREATE_SEQ_DICTIONARY.out)



// step 7 Baserecalibration when one has their knownsite
if (params.knownsites)
BASERECALIBRATION( known_ch, REMOVE_DUPLICATES.out.marked_dups, reference_ch, CREATE_SEQ_DICTIONARY.out)

// step 8 variant calling
if (params.knownsites)
VARIANT_CALL(BASERECALIBRATION.out.recal_bam,  reference_ch, CREATE_SEQ_DICTIONARY.out)
if (!params.knownsites)
VARIANT_CALL(BASERECALIBRATION_2.out.recal_bam,  reference_ch, CREATE_SEQ_DICTIONARY.out)



// step 9 variant filter
VARIANT_FILTER(VARIANT_CALL.out)

// step 10 Normalization
NORMALIZATION(VARIANT_FILTER.out, reference_ch, CREATE_SEQ_DICTIONARY.out)

// step 11 annotation
ANNOTATION(NORMALIZATION.out, snpeff_ch)

}
