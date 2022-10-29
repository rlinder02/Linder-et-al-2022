params.outdir = 'results'

process COVERAGE {
	tag "Calculate average genomewide coverage"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path snp_cov
	path offsets
	path repeats
	path helper 

	output:
	path "Tables/Genomewide_plotting_coverages_2kb_window"
	path "Tables/Genomewide_avg_coverages"
	path "Tables/Treatment_vs_base_coverage_tables_2kb"

	script:

	"""
	mkdir -p Tables/Genomewide_plotting_coverages_2kb_window
	mkdir -p Tables/Genomewide_avg_coverages
	mkdir -p Tables/Treatment_vs_base_coverage_tables_2kb
	preprocess_step2-coverage_genomewide_calculations.R "${snp_cov}" "${offsets}" "${repeats}" "${helper}"
	"""
}