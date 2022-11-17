params.outdir = 'results'

process MUT_GO {
	tag "Annotate mutations with GO terms"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path mutations
	path helper

	output:
	path "SNVsAnnotatedAbr.txt"

	script:

	"""
	analysis_step15c_annotate_mutations_KEGG_GO.R "${mutations}" "${helper}"
	"""
}