params.outdir = 'results'

process PREP_COLLAPSE {
	tag "Assemble tables needed to collapse founder haplotypes"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path reps_using
	path founders
	path haplotypes
	path treatments
	path treatment_key
	path helper

	output:
	path "SEE01_chems_revised.txt"
	path "SEE01_reps_using_haps.txt"

	script:

	"""
	preprocess_step7a_normalizing_collapsed_founders.R "${reps_using}" "${founders}" "${haplotypes}" "${treatments}" "${treatment_key}" "${helper}"
	"""
}