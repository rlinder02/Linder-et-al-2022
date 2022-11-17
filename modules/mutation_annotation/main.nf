params.outdir = 'results'

process MUT_ANNOTATE {
	tag "Annotate mutations using SNPeff"
	//label 'big_mem'
	publishDir params.outdir, mode:'copy'
	debug true

	input:
	path mutations

	output:
	file "mutation_list/*_filtered2.vcf"

	shell:

	'''
	mkdir mutation_list
	mv !{mutations} mutation_list
	cp -rp /usr/local/share/snpeff-5.1-2/* mutation_list
	cp -rp /usr/local/share/snpsift-5.1-0/* mutation_list
	cd mutation_list
	java -Xmx8g -jar snpEff.jar Saccharomyces_cerevisiae !{mutations} > SNP_mutations_ann.vcf
	cat SNP_mutations_ann.vcf | ./scripts/vcfEffOnePerLine.pl > SNP_mutations_ann_split.vcf
	java -jar SnpSift.jar filter "!((ANN[*].EFFECT = 'upstream_gene_variant') & (ANN[*].DISTANCE > 1000))" SNP_mutations_ann_split.vcf > SNP_mutations_ann_filtered.vcf
	java -jar SnpSift.jar filter "!((ANN[*].EFFECT = 'downstream_gene_variant') & (ANN[*].DISTANCE > 200))" SNP_mutations_ann_filtered.vcf > SNP_mutations_ann_filtered2.vcf
	'''
}