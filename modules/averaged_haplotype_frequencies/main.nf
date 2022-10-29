params.outdir = 'results'

process AVERAGE_HAPS {
	tag "Averaging evolved replicate haplotype frequencies"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path coll_haplotypes

	output:
	file "Tables/Averaged_and_sd_tx_tables/*_DT.txt"

	script:

	"""
	mkdir -p Tables/Averaged_and_sd_tx_tables
	analysis_step10_calc_avgs_difs_collapsed_haps_hpc.R "${coll_haplotypes}"
	"""
}