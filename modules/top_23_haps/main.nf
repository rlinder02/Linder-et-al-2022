params.outdir = 'results'

process TOP_23_HAPS {
	tag "Finding the top 23 most over-represented founder/synthetic haplotypes"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path avg_haplotypes
	path helper

	output:
	path "topHapCombosAllChems.txt"
	
	script:

	"""
	analysis_step12-top23_haplotypes.R "${avg_haplotypes}" "${helper}"
	"""
}