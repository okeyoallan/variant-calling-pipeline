// Enable DSL 2 syntax

nextflow.enable.dsl = 2

// import modules here

include { QUALITY_CHECK; MULTIQC; TRIMMOMATIC; POST_FASTQC; MULTIQC_P; ALIGNMENT; MERGE_SAM; CONVERT_TO_BAM; REMOVE_DUPLICATES; CREATE_SEQ_DICTIONARY; BASERECALIBRATION; VARIANT_CALL; VARIANT_FILTER; DECOMPOSITION; ANNOTATION } from "./Modules/gatkHC.nf"

// set input channels

/*  Channel.fromFilePairs( params.reads, checkExists:true )
        .set { read_pairs_ch }

Channel.fromPath ( params.genome, checkIfExists:true )
        .set { reference_ch }

Channel.fromPath (params.variants, checkIfExists:true )
	.set { known_ch }

Channel.fromPath ( params.adapter, checkIfExists:true )
       .set { adapter_ch }
*/
  read_pairs_ch = Channel.fromFilePairs ( params.reads, checkExists:true )
  reference_ch = Channel.fromPath ( params.genome, checkIfExists:true )
  known_ch  = Channel.fromPath ( params.variants, checkIfExists:true )
  adapter_ch = Channel.fromPath ( params.adapter, checkIfExists:true )
  
   
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

// step 7 Baserecalibration
BASERECALIBRATION(known_ch, REMOVE_DUPLICATES.out.marked_dups, reference_ch, CREATE_SEQ_DICTIONARY.out)

// step 8 variant calling
VARIANT_CALL(BASERECALIBRATION.out.recal_bam, reference_ch, CREATE_SEQ_DICTIONARY.out)

// step 9 variant filter
VARIANT_FILTER(VARIANT_CALL.out)

// step 10 decomposition
DECOMPOSITION(VARIANT_FILTER.out, reference_ch)

// step 11 annotation
ANNOTATION(DECOMPOSITION.out)

}
