#!/usr/bin/env nextflow

/*
 * Pipeline for reproducing the figures and tables presented in Linder et al 2022 implemented with Nextflow
 * 
 * Author:
 * - Robert Linder <rlinder02@gmail.com>
 */

nextflow.enable.dsl = 2


/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.haployptes = "$projectDir/assets/Jan2022.allhaps.txt.gz"
params.founders = "$projectDir/assets/Founder_names.txt"
params.snp_count = "$projectDir/assets/SNP.count.Dec28.txt.gz"
params.snp_freqs = "$projectDir/assets/SNP.freq.Dec28.txt.gz"
params.offsets = "$projectDir/assets/newoffsets.txt"
params.repeats = "$projectDir/assets/repeatpos2.txt"
params.helper = "$projectDir/assets/helper_functions.R"
params.metadata = "$projectDir/assets/sexual_paper_sample_metadata_revised.txt"
params.outdir = "$projectDir/results"

// log.info """\
//  Linder et al 2022 - N F   P I P E L I N E
//  ===================================
//  genome: ${params.genome}
//  annotations  : ${params.annotations} 
//  reads        : ${params.reads}
//  outdir       : ${params.outdir}
//  """

include {SNP_REFORMAT} from './modules/snp_reformat/main'
include {COVERAGE} from './modules/genomewide_coverage/main'
include {HAP_REFORMAT} from './modules/haplotype_reformat/main'
include {HAP_SPLIT} from './modules/splitting_haplotypes/main'
include {HAP_HET_DIFFS} from './modules/haplotype_heterozygosity_differences/main'
include {HETEROZYGOSITY} from './modules/heterozygosity_analysis/main'
include {PREP_COLLAPSE} from './modules/Preparing_collapsing_files/main'
include {COLLAPSE} from './modules/collapse_haplotypes/main'
include {TRANSFORM_HAPS} from './modules/transforming_haplotype_freqs/main'
include {TBLE_55_REPS} from './modules/Output_55_replicates_using/main'
include {CHISQ} from './modules/genomewide_chisq_test/main'
include {AVERAGE_HAPS} from './modules/averaged_haplotype_frequencies/main'
include {SNP_DIFFS} from './modules/absolute_snp_differences/main'
include {TOP_23_HAPS} from './modules/top_23_haps/main'
include {CORRELATION_PREP} from './modules/correlations/main'
include {CORRELATION_CALC} from './modules/calculating_correlations/main'
include {CORRELATIONS_COMBINE} from './modules/combine_correlations/main'

/* 
* main script flow
*/
	
workflow {
    SNP_REFORMAT(params.snp_count, params.snp_freqs, params.offsets)
    snp_cov = SNP_REFORMAT.out.map {  it.listFiles()[1]  }
    COVERAGE(snp_cov,  params.offsets, params.repeats, params.helper)
    HAP_REFORMAT(params.haployptes, params.founders, params.offsets, params.helper)
    low_cov = COVERAGE.out[1].map { it.listFiles()}.flatten().filter( ~/.*low_cov_samples_DT.txt$/ )
    HAP_SPLIT(HAP_REFORMAT.out, low_cov, params.founders, params.helper)
    treatment_key = HAP_SPLIT.out[1]
    HAP_HET_DIFFS(HAP_SPLIT.out[2], params.founders, params.helper)
    treatments = HAP_SPLIT.out[0].flatten().filter( ~/.*SEE01_unique_samples.txt$/ )
    HETEROZYGOSITY(HAP_HET_DIFFS.out, params.founders, treatments, treatment_key, params.helper)
    reps_using = HETEROZYGOSITY.out.map { it.listFiles() }.flatten().filter( ~/.*SEE01_Reps_Using_Hets.txt$/ )
    PREP_COLLAPSE(reps_using, params.founders, HAP_REFORMAT.out, treatments, treatment_key, params.helper)
    chemical_ch = PREP_COLLAPSE.out[0].splitCsv( header: false, sep: '\t' ).map { row -> row }
    COLLAPSE(PREP_COLLAPSE.out[1].first(), chemical_ch, treatment_key.first(), params.founders, params.helper)
    TRANSFORM_HAPS(COLLAPSE.out, params.founders, params.helper)
    TBLE_55_REPS(TRANSFORM_HAPS.out.collect(), params.helper) 
    CHISQ(TRANSFORM_HAPS.out, params.founders, params.helper)
    AVERAGE_HAPS(COLLAPSE.out)
    avg_snp_freqs = SNP_REFORMAT.out.map {  it.listFiles()[0]  }
    SNP_DIFFS(avg_snp_freqs, params.founders, params.helper)
    TOP_23_HAPS(AVERAGE_HAPS.out.collect(), params.helper)
    CORRELATION_PREP(HAP_HET_DIFFS.out, params.metadata, params.helper)
    split_ch = CORRELATION_PREP.out[0].splitCsv( header: false, sep: '\t' ).map { row -> row }
    CORRELATION_CALC(split_ch, CORRELATION_PREP.out[1].first(), CORRELATION_PREP.out[2].first())
    CORRELATIONS_COMBINE(CORRELATION_CALC.out.collect())
}

/* 
 * completion handler
 */

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Find the results here: --> $params.outdir\n" : "Oops .. something went wrong" )
}