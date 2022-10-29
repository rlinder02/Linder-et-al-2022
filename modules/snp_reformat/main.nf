params.outdir = 'results'

process SNP_REFORMAT {
	tag "Reformat SNP tables for downstream analyses"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	//debug true

	input:
	path snp_count
	path snp_freqs
	path offsets 

	output:
	path "reformatted_snps"

	script:

	"""
	mkdir reformatted_snps
	preprocess_step1-reformat_snp_files.R "${snp_count}" "${snp_freqs}" "${offsets}"
	mv SNPtable.* reformatted_snps
	"""
}