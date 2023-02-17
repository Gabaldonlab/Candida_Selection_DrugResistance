#!/usr/bin/env python

# This script will integrate all the filtered variant calls into a single table, which has also the other samples that have that variant. It should be run on perSVade_env (version 0.9)

##### DEFINE ENVIRONMENT #######

# module imports
import os
import sys
from Bio import SeqIO
import numpy as np
import multiprocessing as multiproc
import pandas as pd

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    run_in_cluster = False    
    threads = 4
else:
    run_in_cluster = True    
    ParentDir = "/gpfs/projects/bsc40/mschikora"
    threads = 48

# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir
sys.path.insert(0, perSVade_dir)

# import functions
print("import functions")
import sv_functions as fun

# define dirs
CurDir = "%s/CandidaMine_data_generation/v1"%ParentDir
DataDir = "%s/data"%CurDir
PlotsDir = "%s/plots"%CurDir

# define vars
taxID_to_sciName = { 5476:"Candida_albicans", # albicans
                     5480:"Candida_parapsilosis", # parapsilosis
                     273372:"Candida_metapsilosis", # metapsilosis --> the length is 236193, which makes sense
                     273371:"Candida_orthopsilosis", # orthopsilosis
                     5482:"Candida_tropicalis", # tropicalis
                     5478:"Candida_glabrata", # glabrata
                     498019:"Candida_auris" # auris
                        }

taxID_to_ploidy = {5476: 2, # albicans
                   5480: 2, # parapsilosis
                   273372: 2, # metapsilosis
                   273371: 2, # orthopsilosis
                   5482: 2, # tropicalis
                   5478: 1, # glabrata
                   498019: 1 # auris
                   }  

# define the number of samples
taxID_to_nSamples = {5476: 642, # albicans (all samples would be 645, three were removed)
                   5480: 51, # parapsilosis
                   273372: 26, # metapsilosis
                   273371: 33, # orthopsilosis
                   5482: 89, # tropicalis
                   5478: 420, # glabrata
                   498019: 754 # auris (all samples would be 763, 9 were removed because the bam files were corrupted)
                   }  

# define SRRs to remove because some of the variant calling steps could not be performed
srrs_to_remove = ['ERR570064', 'ERR331060'] + ['ERR246509', 'ERR246505', 'ERR246510'] + ['ERR321926'] + ['ERR331061', 'ERR570065'] + ["SRR5083834", "SRR1811019", "SRR1106666"] + ["SRR10461151", "SRR7909148", "SRR7909184", "SRR7909183", "SRR7909204", "SRR7909181", "SRR7909233", "SRR7909235", "SRR7909298"] 

# define the interesting taxIDs 
interesting_taxIDs = [5476, 5480, 273371, 5482, 5478, 498019] # all (excluding metapsilosis (273372))

####################################

# go through each taxID
for taxID in interesting_taxIDs:

    # get a folder to store everything
    sciName = taxID_to_sciName[taxID]
    taxID_dir = "%s/%s_%i"%(DataDir, sciName, taxID); fun.make_folder(taxID_dir)
    integrated_calls_outdir = "%s/integrated_varcalls"%taxID_dir
    print(sciName)

    # define vars
    ploidy = taxID_to_ploidy[taxID]
    reference_genome = "%s/genome.fasta"%taxID_dir
    
    if sciName in {"Candida_tropicalis", "Candida_parapsilosis"}: gff = "%s/annotations.gff.added_mtDNA_annotations.gff"%taxID_dir
    else: gff = "%s/annotations.gff"%taxID_dir

    # get coverage df
    print("filtering df")
    min_coverage = 40
    min_pct_covered = 90
    coverage_df = fun.get_tab_as_df_or_empty_df("%s/coverage_all_srrs.tab"%taxID_dir)
    coverage_df_filt = coverage_df[(coverage_df.mean_coverage>=min_coverage) & (coverage_df.pct_covered>=min_pct_covered) & ~(coverage_df.srr.isin(srrs_to_remove))]

    #### define a df with the varcalls ######
    data_dict = {}
    
    # go through each srr
    for Is, srr in enumerate(coverage_df_filt.srr):
        print(Is, srr)
        perSVade_outdir = "%s/varCall_output/%s"%(taxID_dir, srr)

        # define the final file and check that it was generated
        final_file = "%s/perSVade_finished_file.txt"%perSVade_outdir
        #final_file = "%s/smallVars_CNV_output/variant_annotation_ploidy2.tab"%perSVade_outdir # debug (only small vars)
        if fun.file_is_empty(final_file): raise ValueError("%s was not generated"%final_file)

        data_dict[srr] = {"perSVade_outdir":perSVade_outdir, "srr":srr}

    paths_df = pd.DataFrame(data_dict).transpose().sort_values(by="srr")
    paths_df["sampleID"] = list(range(0, len(paths_df)))
    paths_df["sampleID"] = paths_df["sampleID"].apply(str)

    # check that there is an expected number of samples
    #print("There are %i jobs in paths_df"%len(paths_df))
    if len(paths_df)!=taxID_to_nSamples[taxID]: raise ValueError("There should be %i samples in %i. There are %i samples"%(taxID_to_nSamples[taxID], taxID, len(paths_df)))

    # save
    sampleIDmapping = "%s/samples_data.tab"%taxID_dir
    fun.save_df_as_tab(paths_df, sampleIDmapping)

    #########################################


    # get integrated CNV calls
    integrated_CNperWindow_file = fun.get_integrated_CNperWindow_df_severalSamples(paths_df, integrated_calls_outdir, threads=threads)

    # get integrated SV and CNV calls
    fields_varCall_SV_CNV = ["#CHROM", "POS", "ID", "REF", "ALT", "INFO_BPS_TYPE", "INFO_BREAKENDIDs", "INFO_BREAKEND_FILTER", "INFO_BREAKEND_QUAL", "INFO_BREAKEND_coordinates", "INFO_BREAKEND_has_poly16GC", "INFO_BREAKEND_len_inserted_sequence", "INFO_BREAKEND_length_inexactHomology", "INFO_BREAKEND_length_microHomology", "INFO_BREAKEND_overlaps_repeats", "INFO_BREAKEND_real_AF", "INFO_BREAKPOINTID", "INFO_BREAKPOINTIDs", "INFO_END", "INFO_FILTER", "INFO_QUAL", "INFO_QUAL_max", "INFO_QUAL_mean", "INFO_QUAL_min", "INFO_RELCOVERAGE", "INFO_RELCOVERAGE_TO_3", "INFO_RELCOVERAGE_TO_5", "INFO_SVTYPE", "INFO_all_FILTERs", "INFO_any_has_poly16GC", "INFO_any_overlaps_repeats", "INFO_best_FILTER", "INFO_bpIDs", "INFO_has_poly16GC", "INFO_len_inserted_sequence_max", "INFO_len_inserted_sequence_mean", "INFO_len_inserted_sequence_min", "INFO_length_event_max", "INFO_length_event_mean", "INFO_length_event_min", "INFO_length_inexactHomology", "INFO_length_inexactHomology_max", "INFO_length_inexactHomology_mean", "INFO_length_inexactHomology_min", "INFO_length_microHomology", "INFO_length_microHomology_max", "INFO_length_microHomology_mean", "INFO_length_microHomology_min", "INFO_median_coverage_corrected", "INFO_median_relative_CN_AneuFinder", "INFO_median_relative_CN_HMMcopy", "INFO_merged_relative_CN", "INFO_misc", "INFO_overlaps_repeats", "INFO_real_AF", "INFO_real_AF_max", "INFO_real_AF_mean", "INFO_real_AF_min", "INFO_variantID", "INFO_worse_FILTER", "sampleID"]
    
    fun.get_integrated_SV_CNV_df_severalSamples(paths_df, integrated_calls_outdir, gff, reference_genome, threads=threads, fields_varCall=fields_varCall_SV_CNV, fields_varAnnot="all", integrated_CNperWindow_file=integrated_CNperWindow_file, add_overlapping_samples_eachSV=False)

    # get integrated small vars
    fields_varCall = ["#Uploaded_variation", "#CHROM", "POS", "REF", "ALT", "QUAL", "CALLEDALGS", "HC_GT", "INREPEATS", "ISSNP", "NCALLED", "NPASS", "PASSALGS", "bt_GT", "common_GT", "fb_GT", "mean_AD", "mean_DP", "mean_fractionReadsCov_PASS_algs", "mean_fractionReadsCov_called_algs", "relative_CN"]

    # redefine run_ploidy2_ifHaploid
    if ploidy==2: run_ploidy2_ifHaploid = False
    else: run_ploidy2_ifHaploid = True

    fun.get_integrated_small_vars_df_severalSamples(paths_df, integrated_calls_outdir, ploidy, gff, run_ploidy2_ifHaploid=run_ploidy2_ifHaploid, threads=threads, fields_varCall=fields_varCall, fields_varAnnot="all")


print("SUCCESS!! You integrated all the variant calls")