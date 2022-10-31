params.outdir = 'results'

process ANEUPLOIDY_TABLE {
	tag "Assemble a table of aneuploidies detected amongst the 55 replicates used in this study"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path aneuploidies
	path reps_using
	path treatments
	path treatment_key
	path low_cov
	path helper
	path offsets

	output:
	path "Aneuploidy_table_v4.txt"

	script:

	"""
	analysis_step14b-Aneuploidies_Summary_table.R "${aneuploidies}" "${reps_using}" "${treatments}" "${treatment_key}" "${low_cov}" "${helper}" "${offsets}"
	"""
}