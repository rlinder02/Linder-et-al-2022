params.outdir = 'results'

process CORRELATION_CALC {
	tag "Computing genome-wide Spearman correlations of haplotype frequencies"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	val gposition 
	file haplotypes
	file chemPos

	output:
	path "Tables/Correlation_tables/*_cors.txt"

	script:

	"""
	mkdir -p Tables/Correlation_tables
	analysis_step13b_calculating_haplotype_correlations.R "${gposition}" "${haplotypes}" "${chemPos}"
	"""
}