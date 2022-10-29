params.outdir = 'results'

process CORRELATION_PREP {
	tag "Prep for computing genome-wide Spearman correlations of haplotype frequencies"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path hap_diffs
	path metadata
	path helper

	output:
	path "SEE01_gposIdx.txt"
	path "Reformatted_SEE01_hap_freqs_v3.txt"
	path "Unique_chemPos_v2.txt"

	script:

	"""
	analysis_step13a_calculating_haplotype_correlations.R "${hap_diffs}" "${metadata}" "${helper}"
	"""
}