params.outdir = 'results'

process MUT_FURTHER_ANNOTATE {
	tag "Further annotate the list of potential SNVs"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path ann_mutations
	path gff
	path helper
	path old_muts
	path new_muts
	path offsets

	output:
	path "All_mutations_annotated_less_stringent2.txt"

	script:

	"""
	analysis_step15b_annotate_mutations_post_SNPeff.R "${ann_mutations}" "${gff}" "${helper}" "${old_muts}" "${new_muts}" "${offsets}"
	"""
}