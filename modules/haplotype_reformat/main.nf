params.outdir = 'results'

process HAP_REFORMAT {
	tag "Reformat haplotype data"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path haplotypes
	path founders
	path offsets
	path helper

	output:
	path "*.allhaps.restructured.*"

	script:

	"""
	preprocess_step3-reformatting_haplotype_call_data.R "${haplotypes}" "${founders}" "${offsets}" "${helper}"
	"""
}