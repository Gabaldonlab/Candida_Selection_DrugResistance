## General information

This repository contains the code and software environment definitions used to generate the datasets, results, tables and figures shown our paper 'Genome-wide signatures of recent selection and drug resistance across Candida opportunistic pathogens'.

## Code

The folder 'custom_scripts' includes various R and python scripts used to the generate variant calls, selection scores, processed GWAS datasets and figures of the paper. This is custom code, highly related to our datasets and analyses, so that you cannot run it directly on new sequences. However, you can check it to understand how we generated our datasets and results. This folder includes the following scripts:

- 'get_data.py' was used to generate large datasets: variant calling,  tree generation, assignment of resistance, clade inference, clonal cluster inference, obtention of gene functional annotations, SNP distance measurements, ancestral state reconstruction of variants and running of GWAS. This script actually runs some of the others (described below).

- 'integrate_varcalls.py' was used to integrate the variants from different runs.

- 'get_metadata.py' was used to screen the literature for metadata.

- 'get_tree_from_snps_df_resamplingHeteroSNPs.py' and 'get_consensus_tree_with_branchLengths.R' were used to get the resampled trees for diploids.

- 'run_calculate_pairwise_varDifferences_one_sample.py' was used to calculate pairwise SNP distances between samples.

- 'run_function_get_filtering_stats_df_one_species_drug_gwas_method_consistency_btw_pvals.py' was used to calculate how different GWAS parameters/filters work.

- 'run_get_df_diversity_piN_piS_all_steps_function.py' was used to generate the recent selection metrics.

- 'run_jobs_file_in_parallel.py' was used to run ancestral state reconstruction for each variant.

- 'get_figures_tables.py' was used to generate all figures and tables of the paper.

- 'pipeline_gwas_validation.py' was run to perform the validation of GWAS results. It is a straightforward script to understand how we did most of the steps (sample collection, variant calling, tree generation and GWAS).

- 'Cmine_functions.py' is a module with the python functions used by most of the other scripts.

The folder 'manually_curated_data' includes manually-generated tables based on the literature. These were used in 'get_metadata.py' to generate the strain metadata.

The folder 'ancestral_GWAS_pipeline' includes the scripts that constitute the convergence GWAS pipeline, which we applied to each species and drug. This is a standalone pipeline that can be used on any input dataset, so that it is valuable for further projects. This is the function of each script:

- 'get_ASR_mutations.py' generates the Ancestral State Reconstruction (ASR) for a set of variants.

- 'run_pastml.py' runs pastml on a given tree and variant / phenotype to get ASR results.

- 'get_resampled_phenotypesASR.py' generates resampled phenotypes.

- 'get_GWAS_jobs.py' takes the results of 'get_ASR_mutations.py' and 'get_resampled_phenotypesASR.py' and generates the raw GWAS results.

- 'ancestral_GWAS_functions.py' is a module with the python functions used by the other scripts.

The folder 'tree_from_SNPs' includes the pipeline to get trees from SNP data. Check 'tree_from_SNPs/README.md' for more information on how to install and run it.

The folder 'supplementary_tables' contains the supplementary tables in csv format (sepparated by tabs).

## Software environments

We used a combination of conda environments and manually installed software to run all the scripts above. The folder 'conda_envs' includes the .yml files that define these environments:

- Candida_mine_env: base environment, used by most custom scripts.

- ancestral_GWAS_env: environment has the dependencies to run the ASR and the GWAS pipeline.

- augustus_env: environment used to run augustus.

- genometools_env: environment used to run gt to correct the gff provided by augustus.

- InterProScan_env: environment used to run InterproScan.

- phylotools_R_env: environment used to run 'get_consensus_tree_with_branchLengths.R'.

- tree_from_SNPs_env: environment to run the tree reconstruction pipeline.

- sra_tools_env: environment to get the SRA datasets, used in 'pipeline_gwas_validation.py' and 'Cmine_functions.py'.

In addition, we used additional software:

- perSVade (v0.6) to generate all variant calls, available [here](https://github.com/Gabaldonlab/perSVade). To exactly reproduce our results, you may use the one-liner script of v0.6. However, if you want to use perSVade to do variant calling we recommend using the latest version. Note that to install perSVade you need to create a various conda environments, used by various custom_scripts (in our scripts, these environments are called 'perSVade_env_*').

- Pathway Tools API (v25.0) to work with Metacyc annotations. We installed it interactively following [these guidelines](https://bioinformatics.ai.sri.com/ptools/installation-guide/released/index.html).

- Interproscan (v5.52-86.0), available [here](https://interproscan-docs.readthedocs.io/).

