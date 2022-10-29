params.outdir = 'results'

process TRANSFORM_HAPS {
	tag "Transforming founder haplotypes"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path coll_haplotypes
	path founders
	path helper

	output:
	file "Tables/Ind_tx_ind_haps_sum_of_squared_diffs_tables/*_DT.txt"

	script:

	"""
	mkdir -p Tables/Ind_tx_ind_haps_sum_of_squared_diffs_tables
	analysis_step8a-saving_to_disk_ind_treatments_vs_base_differences.R "${coll_haplotypes}" "${founders}" "${helper}"
	"""
}