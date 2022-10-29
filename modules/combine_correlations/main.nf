params.outdir = 'results'

process CORRELATIONS_COMBINE {
	tag "Combine all Spearman correlation tables generated previously into a single table"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path correlations

	output:
	path "complete_chems_all_reps_spearman_cors.txt"

	script:

	"""
	analysis_step13c_combining_haplotype_correlations.R "${correlations}" 
	"""
}