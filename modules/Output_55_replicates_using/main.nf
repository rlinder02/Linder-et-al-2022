params.outdir = 'results'

process TBLE_55_REPS {
	tag "Assemble a table of information about the 55 replicates used in this study"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path transformed_haplotypes
	path helper

	output:
	path "Individual_haplotype_differences_55_selected_populations.txt"

	script:

	"""
	analysis_step8b-saving_to_disk_55_replicates_using_info.R "${transformed_haplotypes}" "${helper}"
	"""
}