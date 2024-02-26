//
// Prepare indices
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX                           } from '../../modules/nf-core/samtools/index/main'

workflow PREPARE_INDICES {
    take:
    input           // channel: [meta, alignment (BAM), []]
    
    main:

    ch_versions = Channel.empty()
    ch_out      = Channel.empty()

    // Determine if INDEX provided
    input.branch{
        is_indexed:  it[0].index == true
        to_index:    it[0].index == false
    }.set{samtools_input}

    // Remove empty INDEX [] from channel
    input_to_index = samtools_input.to_index.map{ it -> [it[0], it[1]] }

    // INDEX BAM/CRAM only if not provided
    SAMTOOLS_INDEX(input_to_index)
    ch_versions    = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_index_files = Channel.empty().mix(SAMTOOLS_INDEX.out.bai, SAMTOOLS_INDEX.out.crai)
    

    // Combine channels
    ch_new = input_to_index.join(ch_index_files)
    ch_out = samtools_input.is_indexed.mix(ch_new)

    // Gather versions of all tools used
    emit:
    genome_bam_bai = ch_out
    versions       = ch_versions
}