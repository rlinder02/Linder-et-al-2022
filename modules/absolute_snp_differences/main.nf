params.outdir = 'results'

process SNP_DIFFS {
	tag "Calculate absolute SNP frequency differences between evolved populations and the base population"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path avg_snp_freqs
	path founders
	path helper

	output:
	path "Tables/Avg_snps_abs_diffs_dfs"

	script:

	"""
	mkdir -p Tables/Avg_snps_abs_diffs_dfs
	analysis_step11_calc_abs_snp_difs.R "${avg_snp_freqs}" "${founders}" "${helper}"
	"""
}