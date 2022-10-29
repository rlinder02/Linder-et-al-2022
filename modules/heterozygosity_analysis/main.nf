params.outdir = 'results'

process HETEROZYGOSITY {
	tag "Assemble a look-up table of replicates using in this study"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path hap_diffs
	path founders
	path treatments
	path treatment_key
	path helper

	output:
	path "Tables/Summary_tables"

	script:

	"""
	mkdir -p Tables/Summary_tables
	preprocess_step6-heterozygosity_filtering_and_further_analyses.R "${hap_diffs}" "${founders}" "${treatments}" "${treatment_key}" "${helper}"
	"""
}