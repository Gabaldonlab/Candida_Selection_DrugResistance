#!/usr/bin/env python

# runs the get_filtering_stats_df_one_species_drug_gwas_method_consistency_btw_pvals function for a chunk of filters specified via command line

# module imports
import os, sys

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    run_in_cluster = False    

else:
    run_in_cluster = True    
    ParentDir = "/gpfs/projects/bsc40/mschikora"
    
# define dirs
CurDir = "%s/CandidaMine_data_generation/v1"%ParentDir
sys.path.insert(0, CurDir)
import Cmine_functions as fun

# parse the cmd line args
df_gwas_af_s_file, species, drug, gwas_method, filters_df_chunk_file, gene_features_df_s_file, df_gwas_af_no_collapsing_file, df_gwas_af_no_collapsing_resamples_file, filtering_stats_df_sppDrug_chunk_file, threads = sys.argv[1:]
threads = int(threads)

# load dfs
print("loading dfs")
df_gwas_af_s = fun.load_object(df_gwas_af_s_file)
filters_df_chunk = fun.load_object(filters_df_chunk_file)
gene_features_df_s = fun.load_object(gene_features_df_s_file)
df_gwas_af_no_collapsing = fun.load_object(df_gwas_af_no_collapsing_file)
df_gwas_af_no_collapsing_resamples = fun.load_object(df_gwas_af_no_collapsing_resamples_file)

# run function
if fun.file_is_empty(filtering_stats_df_sppDrug_chunk_file):

    # get df
    print("running function")
    filtering_stats_df_sppDrug_chunk = fun.get_filtering_stats_df_one_species_drug_gwas_method_consistency_btw_pvals(0, 0, df_gwas_af_s, species, drug, gwas_method, filters_df_chunk, gene_features_df_s, df_gwas_af_no_collapsing, df_gwas_af_no_collapsing_resamples, threads)

    # checks
    if sorted(filtering_stats_df_sppDrug_chunk.filter_I)!=sorted(filters_df_chunk.filter_I): raise ValueError("there should be one filter as in filter_I")

    print("saving")
    fun.save_object(filtering_stats_df_sppDrug_chunk, filtering_stats_df_sppDrug_chunk_file)

print("success job finished")

