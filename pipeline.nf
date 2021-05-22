#!/usr/bin/env nextflow

// Script parameters
db        = "$baseDir/dataset/Cervix/"
input     = "$baseDir/input.tsv"

// Tool paths
input_to_table   = "$baseDir/input_to_table.R"
norm_DE          = "$baseDir/DE.R"

// Channel creation (optional???)
design_ch        = channel.create()
count_ch         = channel.create()

process input_to_table {
	//conda 'bioconda::r-r.utils=1.1.2'
	conda 'r-r.utils'
	
	input:
	//input.tsv
	
	//htseq dir
	
	output:
	//design (geen file?)
	file design into design_ch

	//htseq count table
	file counts_table into count_ch

	script:
	//Rscript ./input_to_table.R '/home/gebruiker/input.tsv' '/home/gebruiker/dataset/Cervix/'
	"""
	Rscript $input_to_table $baseDir/input.tsv $db
	"""
}

process normalisation_DE {
	//conda 'r.utils' 'edger' 'ggplot2'
	
	input:
	//design 
	file design from design_ch
	
	//counts_table
	file counts_table from count_ch
	
	//output:
	//results
	//rdata als normalisatie apart

	script:
	//Rscript ./normalization.R '/home/gebruiker/design.csv' '/home/gebruiker/dataset/Cervix/'
	"""
	Rscript $norm_DE $input $db
	"""
}
