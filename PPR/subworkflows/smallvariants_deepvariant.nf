include { VG_VIEW } from '../modules/vg_view.nf'
include { VG_INDEX_XG } from '../modules/vg_index_xg.nf'
include { VG_SURJECT } from '../modules/vg_surject.nf'
include { BAM_PREPROCESSING } from '../modules/bam_preprocessing.nf'
include { DEEPVARIANT } from '../modules/deepvariant.nf'

workflow SMALLVARIANTS_DEEPVARIANT {
    take:
    gam_ch      // tuple val(sample_name), path(gam_file)
    gbz_ch      // tuple val(sample_name), path(gbz_file)

    main:
    // Convert GBZ to VG format
    VG_VIEW(gbz_ch)

    // Create XG index
    VG_INDEX_XG(VG_VIEW.out.vg)

    // Surject GAM to BAM
    VG_SURJECT(gam_ch, VG_INDEX_XG.out.xg)

    // BAM preprocessing (sort, header fix, index)
    BAM_PREPROCESSING(VG_SURJECT.out.bam)

    // DeepVariant calling
    DEEPVARIANT(BAM_PREPROCESSING.out.bam)

    emit:
    vcf = DEEPVARIANT.out.vcf
    gvcf = DEEPVARIANT.out.gvcf
    bam = BAM_PREPROCESSING.out.bam
}
