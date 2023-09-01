#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
`
@author: mschikora

This script plots all the figures and tables for the recent selection paper.
"""

#%% ENV

# module imports
import os, sys, scipy
from Bio import SeqIO
import numpy as np
import multiprocessing as multiproc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
import copy as cp
from ete3 import Tree

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    run_in_cluster = False 
    threads = 4   
else:
    run_in_cluster = True    
    ParentDir = "/gpfs/projects/bsc40/mschikora"
    threads = 48
    
# import functions
CurDir = "%s/CandidaMine_data_generation/v1"%ParentDir
sys.path.insert(0, CurDir)
import Cmine_functions as fun

# define dirs
PlotsDir = "%s/Plots_RecentEvolPaper"%CurDir; fun.make_folder(PlotsDir)
PlotsDir_allCmine = "%s/plots"%CurDir
TablesDir = "%s/Tables_RecentEvolPaper"%CurDir; fun.make_folder(TablesDir)
ProcessedDataDir_allCmine = "%s/processed_data"%CurDir # this was generated with ./get_data.py
ProcessedDataDir = "%s/ProcessedData_RecentEvolPaper"%CurDir; fun.make_folder(ProcessedDataDir)
manually_curated_data = "%s/manually_curated_data"%CurDir
DataDir = "%s/data"%CurDir



#%% LOAD PROCESSED DATASETS

# map each species to different files or objects
print("loading objects")
species_to_ref_genome = fun.get_species_to_ref_genome(CurDir)
species_to_gff = fun.get_species_to_gff(CurDir)
species_to_srr_to_sampleID = fun.get_species_to_srr_to_sampleID(CurDir)
species_to_tree = fun.get_species_to_tree(CurDir)

# load the metadata 
print("metadata")
metadata_df = fun.load_object("%s/metadata_df_with_resistance_and_others.py"%ProcessedDataDir_allCmine)

# load the df with pairwise differences
print("pairwise comparisons")
type_comp_to_df_pairwise_diff = {type_comp : fun.get_df_pairwise_diversity(DataDir, ProcessedDataDir, species_to_ref_genome, species_to_gff, metadata_df, species_to_tree, threads=threads, type_comparison=type_comp) for type_comp in ["number_genes_protAltering", "number_vars"]}

# get the repeats df for all species without overlaps and annotated SVs
print("repeats df")
repeats_df = fun.get_repeats_df_all_species(DataDir, ProcessedDataDir, PlotsDir, species_to_ref_genome, threads=threads, replace=False)

# get a df with all the unique filtered SVs and metadata about the mechanisms of SV formation
print("SVs")
unique_SV_CNV_df = fun.get_unique_SV_CNV_df_with_metadata(repeats_df, DataDir, ProcessedDataDir, species_to_ref_genome, threads=threads, replace=False)

# get a df with all the unique filetered small vars and some metadata
print("small vars")
unique_small_vars_df = fun.get_unique_small_vars_df_with_metadata(repeats_df, DataDir, ProcessedDataDir, species_to_ref_genome, threads=threads, replace=False)


# get gene features
gene_features_df = fun.load_object("%s/all_species_gene_features_df.py"%ProcessedDataDir_allCmine)

# load the df with coverage per gene
print("cov per gene")
df_coverage_per_gene = fun.load_object("%s/df_coverage_per_gene_all.py"%ProcessedDataDir_allCmine)

# define general parms
min_pct_overlap_CNV = 25 # the minimum pct to consider that a CNV overlaps repeats
max_fraction_strains_extreme_piS = 0.1 # the maximum number of strains w/ extreme piS
selection_enrichments_max_fraction_genes_pathway_toBeConsidered = 0.25

# get diversity (selection) data for all genes and types
df_diversity_all = fun.get_df_diversity_piN_piS(gene_features_df, unique_SV_CNV_df, unique_small_vars_df, DataDir, "%s/getting_df_diversity_all_piN_piS"%ProcessedDataDir, species_to_gff, df_coverage_per_gene, species_to_srr_to_sampleID, species_to_ref_genome, repeats_df, metadata_df, species_to_tree, PlotsDir, min_pct_overlap_CNV_simpleRepeats=min_pct_overlap_CNV, threads=threads, replace=False)

# get the filtered df
df_diversity_filt = fun.get_df_diversity_filt_piN_piS(df_diversity_all, max_fraction_strains_extreme_piS)


# set the GWAS parms
min_npheno_transitions_GWAS = 5 # min # phenotype transitions to consider a GWAS dataset
type_muts = "only_non_syn" # type of mutations

# get a dict that maps each species and drug to the raw GWAS df file
spp_drug_to_gwas_df_file = fun.get_df_gwas_af_filtering_integrated(DataDir, ProcessedDataDir, threads)

# get the filtering df (I plot twice so that the final is the one used)
for If, (max_n_genes_gwas, max_p_nsignificant_vars, discard_non_bonferroni) in enumerate([(100, 1, False), (100, 0.05, True)]):

    
    filtering_stats_df = fun.get_gwas_filtering_df(spp_drug_to_gwas_df_file, ProcessedDataDir, PlotsDir, gene_features_df, DataDir, threads, species_to_gff, min_npheno_transitions_GWAS, max_n_genes_gwas, max_p_nsignificant_vars)
        
    # Figure GWAS_ChooseFilters (panel with all filters)
    filename = '%s/figure_GWAS_ChooseFilters_allFilters_ngenes<%i_pN<%s_discNonBonf=%s.png'%(PlotsDir, max_n_genes_gwas, max_p_nsignificant_vars, discard_non_bonferroni)
    if fun.file_is_empty(filename): fun.get_figure_GWAS_ChooseFilters_heatmap(filtering_stats_df, filename, "", (10, 11), max_n_genes_gwas, discard_non_bonferroni=discard_non_bonferroni)
    
    # save df all (no filtering)
    if If==0: filtering_stats_df_noMT = cp.deepcopy(filtering_stats_df)

print("final filts GWAS", max_n_genes_gwas, max_p_nsignificant_vars)

# keep only these
valid_spp_drug = set(filtering_stats_df[["species", "drug"]].drop_duplicates().apply(tuple, axis=1))

spp_drug_to_gwas_df_file = {x:f for x,f in spp_drug_to_gwas_df_file.items() if x in valid_spp_drug}

# get the filtered filters, which are best for each spp-drug
filtering_stats_df_filt, designed_GWAS_filters_df  = fun.get_filtering_stats_df_good_filters(filtering_stats_df, TablesDir)

# get the filtered NR gwas results (only protein altering)
df_gwas_filt = fun.get_df_gwas_filtered_and_redundancy_removed_eachVarOnlyInOneGroup(spp_drug_to_gwas_df_file, ProcessedDataDir, PlotsDir, threads, DataDir, designed_GWAS_filters_df, type_muts, species_to_gff, gene_features_df)

# filter and only keep NR results
df_gwas_filt = df_gwas_filt[(df_gwas_filt.type_collapsing.isin({"none", "domains", "genes"})) | ( (~df_gwas_filt.pathway_has_sig_genes) & (df_gwas_filt.pathway_is_NR) )]


# create a df that has the low-confidence GWAS results
gwas_table_df_low_confidence = fun.get_df_gwas_results_low_confidence(spp_drug_to_gwas_df_file, ProcessedDataDir, DataDir, type_muts,  species_to_gff, gene_features_df, threads, TablesDir)

#%% GENERATE TABLES

# Table allStrainData
df_allStrainData = fun.generate_table_allStrainData(metadata_df, species_to_srr_to_sampleID, species_to_tree, TablesDir, DataDir, ProcessedDataDir, threads)

# Table Strain_metadata
df_Strain_metadata = fun.generate_table_Strain_metadata(metadata_df, TablesDir, type_comp_to_df_pairwise_diff)

# Table genes
gene_features_df[["species", "gff_upmost_parent", "gene_name", "Scerevisiae_orthologs", "orthofinder_orthocluster", "description"]].to_excel("%s/gene_features.xlsx"%TablesDir, index=False)

# Table GenesSelection (all genes and nonshared ones). One table for either selected genes or valid genes
selection_dfs_dict = {}
for gene_types in ["under_selection", "all_genes"]: 
    selection_dfs_dict[gene_types] = fun.get_table_GenesSelection(gene_features_df, df_diversity_filt, TablesDir, species_to_gff, DataDir, ProcessedDataDir, gene_types)
    
# Table EnrichmentsSelection
df_enrichment_all = fun.get_table_enrichments_selection(DataDir, gene_features_df, ProcessedDataDir, PlotsDir, df_diversity_filt, species_to_gff, selection_enrichments_max_fraction_genes_pathway_toBeConsidered, TablesDir)

# table with GWAS metadata
df_table_drugs_gwas = fun.get_table_Drugs_for_GWAS(filtering_stats_df, metadata_df, TablesDir, DataDir)

# table with all NR GWAS results
gwas_table_df, gwas_results_df_more_than_1_spp_OGs = fun.get_Table_GroupsGWAS(df_gwas_filt, DataDir, "%s/get_Table_GroupsGWAS_data_nsyn_muts_NR"%ProcessedDataDir, TablesDir, gene_features_df, designed_GWAS_filters_df, type_muts)

#%% GENERATE MERGED TABLES

# Table Selection
fun.get_merged_table_selection(gene_features_df, selection_dfs_dict, df_enrichment_all, TablesDir)

# Table GWAS
fun.get_merged_table_GWAS(gwas_table_df, gwas_results_df_more_than_1_spp_OGs, gwas_table_df_low_confidence, TablesDir)

# Table Strains (only can be run if ./pipeline_gwas_validation.py has been run)
fun.get_merged_table_strains(df_allStrainData, species_to_tree, TablesDir, df_Strain_metadata, df_table_drugs_gwas, CurDir)


#%% ANALYZE SERIAL ISOLATES

# selection analysis
fun.serial_isolates_selection_analysis(metadata_df, ProcessedDataDir, DataDir, PlotsDir, species_to_ref_genome, gene_features_df, species_to_gff, df_coverage_per_gene, selection_dfs_dict["under_selection"][0], selection_dfs_dict["all_genes"][0])

#%% GWAS VALIDATION

# this can only be run after everything else is done, and with previous running of ./pipeline_gwas_validation.py

# keep genes GWAS hits
gwas_table_df_genes = gwas_table_df[gwas_table_df.type_collapsing.isin({"genes", "domains", "none"})]
gwas_results_df_more_than_1_spp_OGs_genes = gwas_results_df_more_than_1_spp_OGs[gwas_results_df_more_than_1_spp_OGs.type_collapsing.isin({"genes", "domains", "none"})]

# add geneID
gwas_table_df_genes["geneID"] = gwas_table_df_genes.apply(fun.get_geneID_gwas_hit, axis=1)
gwas_results_df_more_than_1_spp_OGs_genes["geneID"] = gwas_results_df_more_than_1_spp_OGs_genes.apply(fun.get_geneID_gwas_hit, axis=1)


# get GWAS validation
wrong_spp_drug_combinations = {("Candida_glabrata", "VRC"), ("Candida_auris", "MIF")} # these had <5 pheno transitions, as shown by the previous analysis
df_gwas_validation_all =  fun.get_df_gwas_validation_all(ProcessedDataDir, wrong_spp_drug_combinations, "%s/data_gwas_validation"%CurDir, gwas_table_df_genes)


# plot GWAS validation results

hit_criteria = "best gene hit" # "same grouping" or  "best gene hit"
#fun.plot_gwas_validation_scatterplots(df_gwas_validation_all, ProcessedDataDir, hit_criteria, gwas_table_df_genes, gwas_results_df_more_than_1_spp_OGs_genes, PlotsDir)
fun.plot_gwas_validation_scatterplots_one_pval(df_gwas_validation_all, ProcessedDataDir, hit_criteria, gwas_table_df_genes, gwas_results_df_more_than_1_spp_OGs_genes, PlotsDir, gene_features_df)

    



#%% GENERATE FIGURES

# Figure Dataset_Overview trees with clades 
fun.get_figure_Dataset_Overview_trees_clades_and_donutplots(metadata_df, PlotsDir, ProcessedDataDir, species_to_tree)

# Figure Dataset_Overview species tree
fun.get_figure_Dataset_Overview_species_tree(ProcessedDataDir_allCmine, PlotsDir,  height_rect=300)
fun.get_figure_Dataset_Overview_species_tree(ProcessedDataDir_allCmine, PlotsDir,  height_rect=600)

# Figure DefineBreakpoints
fun.get_figure_DefineBreakpoints(df_allStrainData, PlotsDir, manually_curated_data, metadata_df)

# Figure clade_definition
fun.get_figure_clade_definition(metadata_df, species_to_tree, PlotsDir, ProcessedDataDir)

# Figure Diversity_patterns hmap
ax = fun.get_figure_Diversity_patterns_heatmap(type_comp_to_df_pairwise_diff, "number_vars", PlotsDir, species_to_tree, metadata_df)

# Figure Diversity_patterns boxplots
for type_comparison in type_comp_to_df_pairwise_diff.keys(): fun.get_figure_Diversity_patterns_boxplots(type_comp_to_df_pairwise_diff, type_comparison, PlotsDir, species_to_tree)


# Figure variants_MAF
fun.get_figure_variants_MAF(DataDir, ProcessedDataDir, PlotsDir, species_to_tree)

# Figure mechanisms_SV_formation simple repeats
fun.get_figure_mechanisms_SV_formation_overlap_simple_repeats(unique_SV_CNV_df, unique_small_vars_df, PlotsDir, min_pct_overlap_CNV)
 
# Figure mechanisms_SV_formation actual mechanisms
fun.get_figure_mechanisms_SV_formation_allMechanisms(unique_SV_CNV_df, PlotsDir, repeats_df, min_pct_overlap_CNV, min_fraction_to_qualify_as_repeatSVregion=0.1)


# Figure threshold_for_clonalSamples
fun.get_figure_threshold_for_clonalSamples(metadata_df, PlotsDir, ProcessedDataDir, PlotsDir_allCmine, species_to_tree)

# Figure selection_SNPs: distribution of the fraction of samples tht have an extreme piS
fun.get_selection_SNPs_distribution_fraction_extreme_values(df_diversity_all, "%s/selection_SNPs_fraction_samples_extreme_piN_piS.pdf"%PlotsDir, max_fraction_strains_extreme_piS) 

# Figure selection_SNPs explanation selection score
fun.get_selection_score_scatterplot(df_diversity_filt, "%s/selection_SNPs_selection_score.png"%PlotsDir, ["SNP"], r'$\pi_N > \pi_S$')

# Figure selection_noSNPs explanation selection score
fun.get_selection_score_scatterplot(df_diversity_filt, "%s/selection_noSNPs_selection_score.png"%PlotsDir, ["if_INDEL", "DUP", "DEL"], 'variant')

# Figure selection_SNPs volcano
fun.get_selection_SNPs_selection_score_vs_pval(df_diversity_filt, "%s/selection_SNPs_volcano.png"%PlotsDir)

# Figure selection_noSNPs distribution F value
fun.get_selection_noSNPs_distribution_selecion_score(df_diversity_filt, "%s/selection_noSNP_distribution_score.png"%PlotsDir)

# Figure Selection C. glabrata tree with clades
fun.get_figure_selection_tree_recent_clades(metadata_df, PlotsDir, species_to_tree, 'Candida_glabrata')

# Figure Selection
fun.get_figure_Selection_overlapping_orthogroups(gene_features_df, ProcessedDataDir, PlotsDir, df_diversity_filt, TablesDir, species_to_gff)

# Figure Selection_enrichments
fun.get_figures_Selection_enrichments(df_enrichment_all, selection_enrichments_max_fraction_genes_pathway_toBeConsidered, ProcessedDataDir, PlotsDir, gene_features_df, DataDir)

# Figure GWAS_ChooseFilters (panel with important filters)
fun.get_figure_GWAS_ChooseFilters_heatmap(filtering_stats_df_filt, '%s/figure_GWAS_ChooseFilters_subsetFilters.png'%PlotsDir, "", (10, 11), max_n_genes_gwas, designed_GWAS_filters_df=designed_GWAS_filters_df)

# Figure manhattan plots GWAS
fun.get_figure_GWAS_manhattan(DataDir, ProcessedDataDir, PlotsDir, spp_drug_to_gwas_df_file, threads, species_to_ref_genome, designed_GWAS_filters_df)

# Figure trees GWAS SNP results
df_manhattan = fun.load_object("%s/GWAS_AF_resistance_df_manhattan.py"%(ProcessedDataDir))
fun.get_figure_GWAS_significant_vars_on_trees(DataDir, "%s/data_get_figure_GWAS_significant_vars_on_trees"%ProcessedDataDir, PlotsDir, spp_drug_to_gwas_df_file, threads, species_to_ref_genome, designed_GWAS_filters_df, df_manhattan, metadata_df, gene_features_df)

# Figure GWAS_main distribution resistance across clades
fun.get_figure_GWAS_main_distribution_resistance_across_clades(species_to_tree, metadata_df, spp_drug_to_gwas_df_file, PlotsDir)


# Figure GWAS ASRs phenotypes
fun.plot_GWAS_AFresistance_ASR_phenotypes_tree_withGeneInfo_only_phenotypes(DataDir, PlotsDir, species_to_gff, gene_features_df, metadata_df, designed_GWAS_filters_df)

# Figure with all results
fun.get_figure_GWAS_all_results(PlotsDir, df_gwas_filt)


# Figures with trees affected by GWAS hits
fun.get_figures_trees_with_association_info(DataDir, PlotsDir, species_to_gff, gene_features_df, metadata_df, designed_GWAS_filters_df, df_gwas_filt, spp_drug_to_gwas_df_file)




#%% PRINT NUMBERS 

# number of strains of each type
for t in sorted(set(df_allStrainData.type)): print(t, "%i/%i strains"%(sum(df_allStrainData.type==t), len(df_allStrainData)))

# number of strains with resistance data towards different drugs
for d in sorted([x.split("_")[0] for x in df_allStrainData.keys() if x.endswith("_MIC")]):
    print(d, sum(df_allStrainData[["susceptibility", "intermediate_susceptibility", "resistance"]].apply(lambda r: d in  fun.make_flat_listOflists([(str(x).split(",")) for x in r.values if not pd.isna(x)]), axis=1)))

# number of drugs
print(len([x.split("_")[0] for x in df_allStrainData.keys() if x.endswith("_MIC")]), "total # AF drugs")

# total number of strains
print("total # strains", len(metadata_df))

# print fraction of clinical strains
for spp in fun.sorted_species_byPhylogeny: 
    df_spp = df_allStrainData[df_allStrainData.species_name==spp]
    print("%s: %.2fpct clinical"%(spp, (sum(df_spp.type=="clinical")/len(df_spp))*100  ))

# print fraction of strains with cluster
df_clusters = fun.load_object("%s/get_metadata_df_with_cladeID_clonal_df_plot.py"%ProcessedDataDir)
print(1 - df_clusters[df_clusters["threshold SNPs/kb"]==1].set_index("species")["fraction strains w/out cluster"])

# print number of genes under selection
genes_df_selection = fun.load_object("%s/genes_under_selection.py"%ProcessedDataDir)
print(len(set(genes_df_selection.orthofinder_orthocluster)), "gene families under selection")
for spp in fun.sorted_species_byPhylogeny: print(len(set(genes_df_selection[genes_df_selection.species==spp].orthofinder_orthocluster)), "gene families under selection", spp)

# print number of OGs shared
df = pd.read_excel("%s/GenesSelection.xlsx"%TablesDir)

print("There are %i OGs"%(len(set(df.orthofinder_orthocluster))))

print("There are %i OGs that in some case are shared across diff vars of the same species"%(len(set(df[df.n_types_vars_in_species_orthogroup>=2].orthofinder_orthocluster))))

print("There are %i OGs with selection in different species"%(len(set(df[df.n_species_orthogroup>=2].orthofinder_orthocluster))))

# print total number of OGs
print('total n OGs', len(set(gene_features_df[~pd.isna(gene_features_df.orthofinder_orthocluster)].orthofinder_orthocluster)))

# print number of sites
for taxID, spp in fun.taxID_to_sciName.items():
    taxID_dir = "%s/%s_%i"%(DataDir, spp, taxID)
    if fun.taxID_to_ploidy[taxID]==1: iqtree_files = ["%s/generate_tree_from_SNPs/iqtree_unroted.iqtree"%taxID_dir]
    else: iqtree_files = ["%s/generate_tree_from_SNPs_resamplingHetSNPs/resamplings/resample_%i/iqtree_unroted.iqtree"%(taxID_dir, I) for I in range(1, 101)]
    
    nsites = list(map(fun.get_nsites_iqstree_file, iqtree_files))
    print(spp, "%i-%i sites"%(min(nsites), max(nsites)))
    
# GWAS stats
print("number spp-drug GWAS: %i"%len(spp_drug_to_gwas_df_file))
print("filtering optimization # filters: %i"%(len(filtering_stats_df[(filtering_stats_df.species=="Candida_albicans") & (filtering_stats_df.drug=="FLC")])))

# number of GWAS done
all_types_variants = {'coverageCNV', 'SV', 'SNP', 'INDEL'}

list_types_vars = [("all_vars", all_types_variants),
                ("small_vars", {'SNP', 'INDEL'}), 
                ("SVs", {'SV'}),
                ("coverageCNVs", {'coverageCNV'}),
                ("SVs_and_CNVs", {'coverageCNV', 'SV'}),
                ("small_vars_and_SVs", {'SV', 'SNP', 'INDEL'}),
                ("small_vars_and_CNVs", {'coverageCNV', 'SNP', 'INDEL'})]

nonSyn_types_mutations = {"non_syn_muts", "non_syn_non_truncating_muts", "truncating_muts"}
nonSyn_types_mutations_and_all = {"non_syn_muts", "non_syn_non_truncating_muts", "truncating_muts", "all_muts"}
I=0

for type_vars, type_vars_set in list_types_vars:
    for type_genes in ["all_genes", "only_protein_coding", "only_non_coding"]:
        for type_mutations in ["all_muts", "syn_muts", "non_syn_muts", "non_syn_non_truncating_muts", "truncating_muts"]:
            for type_collapsing in ["none", "genes", "domains", "GO", "MetaCyc", "Reactome"]: 
                if type_collapsing=="none" and (type_vars!="all_vars" or type_genes!="all_genes" or type_mutations!="all_muts"): continue
                if type_collapsing=="genes" and (type_genes!="all_genes" or type_mutations not in nonSyn_types_mutations_and_all): continue
                if type_collapsing in {"domains", "GO", "Reactome", "MetaCyc"} and (type_genes!="only_protein_coding" or type_mutations not in nonSyn_types_mutations): continue
            
                I+=1

print("%i GWAS runs"%I)

# number of orthogroups and pathways in gwas
print(len(set.union(*gwas_table_df[gwas_table_df.orthogroups.apply(lambda x : x.startswith("OG"))].orthogroups.apply(lambda x: set(x.split(","))))), "orthogroups in NR GWAS")

print( len(set(gwas_table_df[gwas_table_df.type_collapsing.isin({"Reactome", "GO", "MetaCyc"})].group_name)), "pathways in NR GWAS")


print("# hits considering >small_vars", sum(gwas_table_df.type_vars!="small_vars")/len(gwas_table_df))

print("# hits considering truncating muts", sum(gwas_table_df.type_mutations=="truncating_muts")/len(gwas_table_df))

for type_col in ["GO", "Reactome", "MetaCyc", "genes", "domains", "none"]:
    
    print("# hits considering %s"%type_col, sum(gwas_table_df.type_collapsing==type_col)/len(gwas_table_df))

print("# hits considering pathways", sum(gwas_table_df.type_collapsing.isin({"GO", "Reactome", "MetaCyc"}))/len(gwas_table_df))


# number orthogroups shared across datasets

print(len(set.union(*gwas_results_df_more_than_1_spp_OGs[gwas_results_df_more_than_1_spp_OGs.orthogroups.apply(lambda x : x.startswith("OG"))].orthogroups.apply(lambda x: set(x.split(","))))), "orthogroups in NR GWAS >1 ds")

print( len(set(gwas_results_df_more_than_1_spp_OGs[gwas_results_df_more_than_1_spp_OGs.type_collapsing.isin({"Reactome", "GO", "MetaCyc"})].group_name)), "pathways in NR GWAS >1 ds")

# print number of hits
for d in ["Candida_auris-ANI", "Candida_glabrata-MIF"]: print("nhits", d, sum(gwas_table_df.species_and_drug==d))

# genes working GWAS
for df, type_filt in [(filtering_stats_df_noMT, "all_filters"), (filtering_stats_df, "MTcorr_filters")]:
    
    int_spp_drug = {x[0]+"-"+x[1] for x in spp_drug_to_gwas_df_file}.difference({"Candida_auris-AMB", "Candida_glabrata-MIF", "Candida_auris-ANI", "Candida_auris-MIF", "Candida_glabrata-VRC", "Candida_auris-VRC", "Candida_glabrata-FLC"})
    df = df[df.spp_and_drug.isin(int_spp_drug)]
    df = df[(df.correct_row_p_nsignificant_vars_and_max_ngenes)]
    
    if type_filt=="MTcorr_filters": df = df[(df.correction_method.isin({"bonferroni"})) | (df.only_maxT_pvals)]

    df["tuple_expected_genes"] = df.genes_expected.apply(sorted).apply(tuple)
    df["tuple_expected_genes_found"] = df.genes_expected_found.apply(sorted).apply(tuple)


    print("\n", type_filt, "\n", df[["spp_and_drug", "tuple_expected_genes", "tuple_expected_genes_found"]].drop_duplicates())
        
# print genes
for spp, gene in [("Candida_auris", "Scer_ERG3")]:
    
    geneID = gene_features_df[(gene_features_df.final_name==gene) & (gene_features_df.species==spp)].gff_upmost_parent.iloc[0]
    print(geneID)
    
    df_s = gwas_table_df_low_confidence[gwas_table_df_low_confidence.species==spp]
    df = df_s[df_s.description.apply(lambda x: geneID in x)].drop_duplicates(subset=["drug", "group_name"]).sort_values(by=["drug", "group_name"])
    
    
    for I,r in df.iterrows(): print(r.drug,r.type_vars, r.type_mutations, r.group_name)
    
# investigate ergosterol genes / pathways
df = gwas_table_df_low_confidence[gwas_table_df_low_confidence.species=='Candida_glabrata']
for I, r in df[df.biological_process_GO.apply(lambda x: "ergosterol" in x)][["drug", "group_name", "description"]].drop_duplicates().iterrows():
    
    print("\n", r.drug, r.description)

# print repeated strains
bs_to_runs = df_allStrainData.groupby("strain").apply(lambda d: set(d.Run))
bs_to_runs_multiple = bs_to_runs[bs_to_runs.apply(len)>1]
print("There are %i/%i strains that are in multiple runs. %.4f"%(len(bs_to_runs_multiple), len(bs_to_runs), len(bs_to_runs_multiple)/len(bs_to_runs)))

# print strains that are equal
df = type_comp_to_df_pairwise_diff["number_vars"]
df = df[df.origin_sample!=df.target_sample]
df_no_dif = df[df["diffent_positions"]==0]
print("There are %i/%i comparisons that are different. %.4f pct"%(len(df_no_dif), len(df), (len(df_no_dif)/len(df))*100))

#%% GRAVEYARD

# check the overlap between the glabrata hits and the genes affected by azole-resistance mutations in Pais 2022
df_GWAS_hits_all_shared = fun.get_tab_as_df_or_empty_df("%s/Tables_RecentEvolPaper/supplementary_tables/TableS3-High-confidence_GWAS_hits_>1_dataset.csv"%CurDir)

glabrata_gene_IDs_gwas = set(df_GWAS_hits_all_shared[(df_GWAS_hits_all_shared.species=="Candida_glabrata") & (df_GWAS_hits_all_shared.drug!="MIF")].apply(fun.get_geneID_gwas_hit, axis=1))
glabrata_gene_IDs_Pais2022 = set(pd.read_excel("%s/manually_curated_data/2022A_Pais-supplementary-table-1.xlsx"%CurDir, "Azole_resistant_exclusive").ORF)
intersecting_genes = glabrata_gene_IDs_gwas.intersection(glabrata_gene_IDs_Pais2022)
print("In Pais 2022 they detect %i genes with mutations exclusively in azole-resistant strains. From our GWAS genes that are in >1 dataset, there are %i/%i genes overlapping this set:\n%s"%(len(glabrata_gene_IDs_Pais2022), len(intersecting_genes), len(glabrata_gene_IDs_gwas), "\n".join(intersecting_genes)))