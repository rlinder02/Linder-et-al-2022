params.outdir = 'results'

process COLLAPSE {
	tag "Collapse founder haplotypes"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	file haplotypes
	val chemical
	file treatment_key
	path founders
	path helper

	output:
	path "Tables/All_reps_tx_tables/*_DT.txt"

	script:

	"""
	mkdir -p Tables/All_reps_tx_tables
	preprocess_step7b_collapsing_haplotypes_hpc.R "${haplotypes}" "${chemical}" "${treatment_key}" "${founders}" "${helper}"
	"""
}