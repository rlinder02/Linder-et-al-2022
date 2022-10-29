params.outdir = 'results'

process HAP_HET_DIFFS {
	//conda 'bioconductor-preprocesscore'
	tag "Calculate haplotype frequency and heterozygosity differences"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path haplotypes
	path founders
	path helper

	output:
	path "Tables/Ind_tx_ind_haps_het_diffs_tables"

	script:

	"""
	mkdir -p Tables/Ind_tx_ind_haps_het_diffs_tables
	preprocess_step5-saving_to_disk_ind_heterozygosity_differences.R "${haplotypes}" "${founders}" "${helper}"
	"""
}