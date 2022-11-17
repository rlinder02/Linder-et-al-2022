params.outdir = 'results'

process MUT_TABLE {
	tag "Assemble the final table of annotated SNVs"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path mutations
	path helper
	path treatments
	path treatment_key
	path low_cov
	path reps_using

	output:
	path "SNVsAnnotated_TFs.txt"

	script:

	"""
	analysis_step15d_annotate_mutations_assemble_final_SNV_table.R "${mutations}" "${helper}" "${treatments}" "${treatment_key}" "${low_cov}" "${reps_using}"
	"""
}