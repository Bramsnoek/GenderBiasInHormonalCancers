#!/usr/bin/env nextflow

// Input parameters
params.db	= "$baseDir/dataset/Cervix/"
params.input	= "$baseDir/input.tsv"
params.site	= "Pancreas"

// Tool paths
input_to_table	= "$baseDir/scripts/merge_files.R"
norm_de	= "$baseDir/scripts/DE.R"

// Format params
site		= "$params.site"

// Initiate channel
dataset_design = Channel.fromPath( 'input.tsv' )

// Pre-proccessing the Htseq counts
process merged_files {
	conda 'conda-forge::r-r.utils=2.7.0'
	
	output:
	file "merged.csv" into csv_channel
	
	script:
	println("Running pipeline for:\t${site}")
	"""
	Rscript ${input_to_table} "${baseDir}/dataset/${site}/" "merged.csv"
	"""
}

// Normalisation and Analysis of the Htseq counts 
process normalization {
	conda 'conda-forge::r-r.utils bioconda::bioconductor-edger r::r-ggplot2 r::r r::r-dplyr conda-forge::r-readr conda-forge::r-cli'
	publishDir "${baseDir}/output/", mode: 'copy'
	
	input:
	file merged_set from csv_channel
	file design from dataset_design
	
	output:
	set file("${site}_lm_results.csv"), file("${site}_lm_results_05.csv") into output
	
	script:
	"""
	Rscript ${norm_de} ${design} ${merged_set} "${site}_lm_results.csv" "${site}_lm_results_05.csv"
	"""
	
}



