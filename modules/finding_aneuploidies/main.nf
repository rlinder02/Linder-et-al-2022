params.outdir = 'results'

process ANEUPLOIDIES {
	tag "Find aneuploidies"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path norm_covs
	path low_cov
	path offsets
	path helper

	output:
	path "Plots/Ind_normalized_coverage_mutations_plots"
	file "Tables/Treatment_vs_base_coverage_tables_2kb/*_hmm_DT.txt"
	file "Tables/Aneuploid_windows/*_aneuploid_windows_DT.txt"

	script:
	"""
	mkdir -p Plots/Ind_normalized_coverage_mutations_plots
	mkdir -p Tables/Treatment_vs_base_coverage_tables_2kb
	mkdir -p Tables/Aneuploid_windows
	analysis_step14-Aneuploidies_calling_HMM.R "${norm_covs}" "${low_cov}" "${offsets}" "${helper}"
	"""
}