params.outdir = 'results'

process HAP_SPLIT {
	tag "Split haplotypes from individual samples into separate files"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path haplotypes
	path low_cov
	path founders
	path helper

	output:
	path "SEE01_unique*"
	path "treatment_key.txt"
	path "Tables/Individual_tx_tables"

	script:

	"""
	mkdir -p Tables/Individual_tx_tables
	preprocess_step4-Splitting_haplotype_files.R "${haplotypes}" "${low_cov}" "${founders}" "${helper}"
	"""
}