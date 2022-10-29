params.outdir = 'results'

process CHISQ {
	tag "Run a chi-squared test genomewide to look for significant haplotype frequency differences from the base popoulation"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path transformed_haplotypes
	path founders
	path helper

	output:
	path "Tables/Hap_adjusted_LOD_tables/*_DT.txt"
	path "Tables/Repeatability_tables/*_DT.txt"
	path "Tables/Pleiotropy_tables/*.txt"

	script:

	"""
	mkdir -p Tables/Hap_adjusted_LOD_tables
	mkdir -p Tables/Repeatability_tables
	mkdir -p Tables/Pleiotropy_tables
	analysis_step9-saving_to_disk_chisq_LOD_scores.R "${transformed_haplotypes}" "${founders}" "${helper}"
	"""
}