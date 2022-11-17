params.outdir = 'results'

process MUT_PREPROCESS {
	tag "Preprocess mutation list for input into SNPeff"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path old_muts
	path new_muts
	path offsets
	path helper

	output:
	path "SNP_mutations_sexual_restructured.vcf"

	script:

	"""
	analysis_step15a_annotate_mutations_preprocess.R "${old_muts}" "${new_muts}" "${offsets}" "${helper}"
	"""
}