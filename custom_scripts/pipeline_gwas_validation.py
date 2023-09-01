#!/usr/bin/env python

# This script is to run GWAS as a validation. It should be run in Candida_mine_env. This script gets all the data from the SRA for C. aurs and C. glabrata that was not analyzed in our GWAS analysis and has some susceptibility data available. This requires running all the data related to the first submission of the paper. Note that the perSVade_env is the one of v1.02.7, which is the same as the one used in the initial paper but modular.

##### DEFINE ENVIRONMENT #######

# module imports
import os
import sys
from Bio import SeqIO
import numpy as np
import multiprocessing as multiproc
import pandas as pd
import copy as cp

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    run_in_cluster = False    
    threads = 4
else:
    run_in_cluster = True    
    ParentDir = "/gpfs/projects/bsc40/mschikora"

# import functions
print("import functions")
CurDir = "%s/CandidaMine_data_generation/v1"%ParentDir
sys.path.insert(0, CurDir)
import Cmine_functions as fun

# define threads
if run_in_cluster is True:

    cluster_name = fun.get_current_clusterName_mareNostrum()
    if cluster_name=="MN4": threads = 48
    elif cluster_name=="Nord3": threads=4
    else: raise ValueError("invalid clustername")

# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir # I used version 1.02.8
#sys.path.insert(0, perSVade_dir)
perSVade_wrapper = "%s/perSVade"%perSVade_dir # I used version 1.02.7
get_trimmed_reads_for_srr_py = "%s/get_trimmed_reads_for_srr.py"%perSVade_dir

# define dirs
DataDir = "%s/data"%CurDir
ProcessedDataDir = "%s/processed_data"%CurDir # to access files generated in the paper
DataDir_gwas_validation = "%s/data_gwas_validation"%CurDir; fun.make_folder(DataDir_gwas_validation)

################################

# load the interesting GWAS hits
df_GWAS_hits_all = fun.get_tab_as_df_or_empty_df("%s/Tables_RecentEvolPaper/supplementary_tables/TableS3-High-confidence_GWAS_hits.csv"%CurDir)

# define the dataframe with the GWAS validation data (only drug resistance data)
df_metadata_new_samples = fun.get_metadata_df_GWAS_validation("%s/fetching_metadata"%DataDir_gwas_validation, CurDir)
fun.save_object(df_metadata_new_samples, "%s/df_metadata_new_samples_GWAS_validation.py"%ProcessedDataDir)

# define the interesting taxIDs
interesting_taxIDs = [5478, 498019] # glabrata and auris
#interesting_taxIDs = [498019] # auris

# define list of cmds, which will be run in the cluster as a greasy array (see section)
list_cmds = [] # this is a tuple of (type_job, cmd, env, threads)

# define the filters of the coverage
min_coverage = 40 # only Runs with >40x average coverage will be considered
min_pct_covered = 90 # only Runs where 90% of the genome is covered will be considered
min_coverage_pos = 12 # only positions with >12x coverage will be considered for variant calling / tree generation

##### OBTAIN DATA IN THE CLUSTER ######

# go through each taxID
for taxID in interesting_taxIDs:

    # define species specific things
    sciName = fun.taxID_to_sciName[taxID]; print(taxID, sciName)
    taxID_dir = "%s/%s_%i"%(DataDir, sciName, taxID)
    mitochondrial_chromosome = fun.taxID_to_mtChromosome[taxID]
    mitochondrial_code = fun.taxID_to_mDNA_code[taxID]
    gDNA_code = fun.taxID_to_gDNA_code[taxID]
    ploidy = fun.taxID_to_ploidy[taxID]
    genome = "%s/genome.fasta"%taxID_dir
    gff = "%s/annotations.gff"%taxID_dir
    taxIDs = fun.taxID_to_taxIDs[taxID]
    repeats_file = "%s.repeats.tab"%genome # run in the main analysis, obtained with perSVade

    # check that the gff is correct and it has mtDNA annotations
    if fun.get_gff_has_mtDNA_annotations(gff, taxID) is False: raise ValueError("%s lacks mtDNA annotations"%gff)

    # everything is only for p==1
    if ploidy!=1: raise ValueError("ploidy should be 1")

    # define the dir where to store all outputs of this pipeline
    species_dir = "%s/%s"%(DataDir_gwas_validation, sciName); fun.make_folder(species_dir)

    # load the metadata df with information of the initial dataset
    df_metadata = fun.load_object("%s/processed_data/metadata_df_with_resistance_and_others.py"%CurDir)
    df_metadata = df_metadata[df_metadata.species_name==sciName]

    # load the GWAS hits, which define interesting genes and drugs, filtering out some drugs where we have not enough power in the original GWAS
    print("defining interesting runs")
    df_GWAS_hits = df_GWAS_hits_all[(df_GWAS_hits_all.species==sciName) & (df_GWAS_hits_all.type_collapsing.isin({"genes", "domains", "none"}))]
    if sciName=="Candida_glabrata": df_GWAS_hits = df_GWAS_hits[df_GWAS_hits.drug!="MIF"]
    elif sciName=="Candida_auris": df_GWAS_hits = df_GWAS_hits[df_GWAS_hits.drug!="ANI"]



    # define the interesting drugs, those that had some GWAS hits
    df_metadata_new_samples_spp = df_metadata_new_samples[df_metadata_new_samples.species==sciName]
    interesting_drugs = sorted(set(df_GWAS_hits.drug))

    # for each drug, get interesting Runs. Only S/R for runs with >=5 S and R strains
    drug_to_interesting_srr_runs = {}
    for drug in interesting_drugs:

        # define SRRs that were used in the previous GWAS's (to discard them)
        sampleIDs_previous_gwas = set(fun.get_tab_as_df_or_empty_df("%s/ancestral_GWAS_drugResistance/GWAS_%s_resistance/resistance_df.tab"%(taxID_dir, drug)).sampleID.apply(int))
        runs_previous_gwas = set(df_metadata[df_metadata.sampleID.apply(int).isin(sampleIDs_previous_gwas)].Run)
        if len(sampleIDs_previous_gwas)!=len(runs_previous_gwas): raise ValueError("lens should be the same")

        # get df
        resistance_f = "%s_resistance"%drug
        df_d = df_metadata_new_samples_spp[(df_metadata_new_samples_spp[resistance_f].isin({"R", "S"})) & ~(df_metadata_new_samples_spp.Run.isin(runs_previous_gwas))]

        # keep samples
        if sum(df_d[resistance_f]=="R")>=5 and sum(df_d[resistance_f]=="S")>=5: drug_to_interesting_srr_runs[drug] = set(df_d.Run)
        else: print("WARINNG: Not >=5 S&R samples for %s-%s"%(sciName, drug))

    interesting_srr_runs = sorted(set.union(*drug_to_interesting_srr_runs.values()))
    print("There are %i/%i SRR runs to run GWAS on"%(len(interesting_srr_runs), len(df_metadata_new_samples_spp)))

    #sprint(sciName, "\n", pd.Series(drug_to_interesting_srr_runs).apply(len)); continue

    # check that the srr runs were not previously analyzed
    previous_srrs = set(os.listdir("%s/varCall_output"%taxID_dir))
    if len(previous_srrs.intersection(interesting_srr_runs))>0: raise ValueError("some srrs were already analyzed")

    # download SRR files (precompiled files with reads)
    print("Downloading SRR files") # THis was failing for some ununderstood reason. I ran one example prefetch command (without parallelizing) and then this code worked
    varcall_dir = "%s/varcall"%species_dir; fun.make_folder(varcall_dir)
    srr_to_SRRfile = {}
    srr_to_outdirReads = {}
    srr_to_varcallDir = {}
    for srr in interesting_srr_runs:
        folder_srr = "%s/%s"%(varcall_dir, srr); fun.make_folder(folder_srr)
        folder_reads = "%s/reads"%folder_srr; fun.make_folder(folder_reads)
        srr_to_outdirReads[srr] = folder_reads
        srr_to_varcallDir[srr] = folder_srr
        srr_to_SRRfile[srr] = "%s/%s.srr"%(folder_reads, srr)

    inputs_downloads = [(srr, srr_to_SRRfile[srr], False) for srr in interesting_srr_runs if fun.file_is_empty(srr_to_SRRfile[srr])]#[0:5]
    with multiproc.Pool(threads) as pool:
        list_srr_files = pool.starmap(fun.download_srr_with_prefetch, inputs_downloads, chunksize=1)
        pool.close()

    # this section keeps appending things to list_cmds. Only append types of commands where the prior cmd has not been completed

    # get reads and map them
    print("getting reads and map commands...")
    for srr in interesting_srr_runs:

        # define read files
        reads1 = "%s/%s_trimmed_reads_1.fastq.gz"%(srr_to_outdirReads[srr], srr)
        reads2 = "%s/%s_trimmed_reads_2.fastq.gz"%(srr_to_outdirReads[srr], srr)

        # get reads
        if fun.file_is_empty(reads2): 
            list_cmds.append(("getReads", "%s --srr %s --outdir %s --threads %i"%(get_trimmed_reads_for_srr_py, srr, srr_to_outdirReads[srr], 8), "perSVade_env", 8))
            continue # do not map reads if reads are not obtained

        # map reads
        outdir_align_reads = "%s/align_reads"%(srr_to_varcallDir[srr])
        if fun.file_is_empty("%s/perSVade_finished_file.txt"%outdir_align_reads):
            list_cmds.append(("alignReads", "%s align_reads -o %s --fraction_available_mem 0.5 --threads 24 --fractionRAM_to_dedicate 0.8 --verbose -r %s --min_chromosome_len 5000 -f1 %s -f2 %s"%(perSVade_wrapper, outdir_align_reads, genome, reads1, reads2), "perSVade_env", 24))


    # if some align reads need to be performed continue
    srr_to_alignReadsDir = {srr : "%s/align_reads"%(srr_to_varcallDir[srr]) for srr in interesting_srr_runs}
    if any([fun.file_is_empty("%s/perSVade_finished_file.txt"%d) for d in srr_to_alignReadsDir.values()]): 
        print("There are still reads to map, skipping further steps")
        continue

    # filter by coverage
    coverage_df = fun.get_df_coverage_per_run_gwas_validation(srr_to_alignReadsDir, sciName, "%s/coverage_per_run.tab"%species_dir, threads)
    coverage_df_filt = coverage_df[(coverage_df.mean_coverage>=min_coverage) & (coverage_df.pct_covered>=min_pct_covered)]
    print("There are %i/%i srrs that have a correct coverage"%(len(coverage_df_filt), len(coverage_df)))
    variant_calling_srrs = sorted(set(coverage_df_filt.srr))

    # run module get_stats_optimization to define the samples in which to run parameter optimization
    file_sorted_bams = "%s/variant_calling_bamfiles.tab"%species_dir
    fun.save_df_as_tab(pd.DataFrame({srr : {"sorted_bam":"%s/aligned_reads.bam.sorted"%(srr_to_alignReadsDir[srr]), "sampleID":srr}  for srr in variant_calling_srrs}).transpose().reset_index(drop=True)[["sampleID", "sorted_bam"]], file_sorted_bams)

    outdir_select_parms = "%s/get_stats_optimization_outdir"%species_dir
    final_file_select_parms = "%s/perSVade_finished_file.txt"%outdir_select_parms
    if fun.file_is_empty(final_file_select_parms):
        print("running get_stats_optimization ")
        fun.run_cmd_anyEnv("%s get_stats_optimization -o %s  --verbose --fraction_available_mem 1. --fractionRAM_to_dedicate 0.8 -thr %i -r %s -mchr %s --min_chromosome_len 100000 --samples_file %s --overlap_coverage 15 --overlap_insert_size 15 --overlap_insert_size_sd 15 --overlap_read_len 15 --rerun_sample_clustering --skip_cleaning_tmpdir"%(perSVade_wrapper, outdir_select_parms, threads, genome, mitochondrial_chromosome, file_sorted_bams), "perSVade_env") # args 

    fun.delete_folder("%s/tmp_bams"%outdir_select_parms)

    # run parameter optimization for representative samples
    df_rep_samples_optimization = fun.get_tab_as_df_or_empty_df("%s/samples_parameter_optimization.tab"%outdir_select_parms)
    representative_srrs = sorted(set(df_rep_samples_optimization.representative_sampleID))
    print("There are %i/%i SRRs with unique parameters. Running parameter optimization on them"%(len(representative_srrs), len(df_rep_samples_optimization)))

    outdir_optmize_parameters_all = "%s/optimize_parameters_rep_samples"%species_dir; fun.make_folder(outdir_optmize_parameters_all)
    all_optimize_parms_finished = True
    for srr in representative_srrs:

        # log
        #print("\n", srr, "\n", df_rep_samples_optimization[df_rep_samples_optimization.representative_sampleID==srr][["sampleID", "gDNA_median_coverage", "mtDNA_median_coverage", "read_length", "median_insert_size", "median_insert_size_sd"]])

        # define outdir
        outdir_optmize_parms = "%s/%s"%(outdir_optmize_parameters_all, srr)
        final_file_optimize_parms = "%s/perSVade_finished_file.txt"%outdir_optmize_parms
        if fun.file_is_empty(final_file_optimize_parms): 
            list_cmds.append(("optimizeParms", "%s optimize_parameters -o %s --fraction_available_mem 1.0 --threads 48 --fractionRAM_to_dedicate 0.85 --verbose -r %s --min_chromosome_len 100000 -sbam %s/aligned_reads.bam.sorted -mchr %s --repeats_file %s --regions_SVsimulations random --simulation_ploidies %s --nvars 50 --nsimulations 2 --range_filtering_benchmark theoretically_meaningful_NoFilterRepeats"%(perSVade_wrapper, outdir_optmize_parms, genome, srr_to_alignReadsDir[srr], mitochondrial_chromosome, repeats_file, {1:"haploid", 2:"diploid_hetero"}[ploidy]), "perSVade_env", 48))
            
            all_optimize_parms_finished = False
            continue

        # print accuracy
        df_acc = fun.get_tab_as_df_or_empty_df("%s/optimized_parameters_accuracy.tab"%outdir_optmize_parms)
        mean_f_value = np.mean(df_acc[df_acc.SV_type=="integrated"].Fvalue)
        if mean_f_value<0.85: print("WARINNG: mean F score in %s-%s is %.3f"%(sciName, srr, mean_f_value))

    # do not contin
    if all_optimize_parms_finished is False: continue

    # run var call of most modules
    print("getting cmds variant calling...")
    small_vars_vcfs = []
    all_vars_generated = True

    for Is,srr in enumerate(variant_calling_srrs):
        print("%s %i/%i"%(srr, Is+1, len(variant_calling_srrs)))

        # run most variant calling steps. Small var calling with p=1, skip annotation as this will be done after p=2
        if ploidy!=1: raise ValueError("this is only for p=1")
        outdir_varcall_general = "%s/general_varcall"%(srr_to_varcallDir[srr])
        parms_file = "%s/%s/optimized_parameters.json"%(outdir_optmize_parameters_all, df_rep_samples_optimization[df_rep_samples_optimization.sampleID==srr].iloc[0].representative_sampleID) # parms of the representative
        sbam_file = "%s/aligned_reads.bam.sorted"%(srr_to_alignReadsDir[srr])

        if fun.file_is_empty("%s/perSVade_finished_file.txt"%outdir_varcall_general): 
            list_cmds.append(("generalVarcall", "%s run_several_modules get_cov_genes,call_small_variants,call_CNVs,call_SVs,integrate_SV_CNV_calls,annotate_SVs %s -r %s -sbam %s --repeats_file %s -p 1 --callers HaplotypeCaller,bcftools,freebayes --min_AF 0.9 --min_coverage %i --min_chromosome_len 100000 --verbose --fraction_available_mem 0.5 --threads 24 --fractionRAM_to_dedicate 0.85 -gff %s --mitochondrial_chromosome %s --cnv_calling_algs HMMcopy,AneuFinder --window_size_CNVcalling 300 --SVcalling_parameters %s --mitochondrial_code %i --gDNA_code %i"%(perSVade_wrapper, outdir_varcall_general, genome, sbam_file, repeats_file, min_coverage_pos, gff, mitochondrial_chromosome, parms_file, mitochondrial_code, gDNA_code), "perSVade_env", 24))
            all_vars_generated = False
            continue

        # run variant calling for p=2  (to also account for some het positions)
        outdir_varcall_p2 = "%s/small_var_call_p2"%(srr_to_varcallDir[srr])
        if fun.file_is_empty("%s/perSVade_finished_file.txt"%outdir_varcall_p2): 
            list_cmds.append(("p2VarCall", "%s call_small_variants -o %s -r %s -sbam %s --repeats_file %s -p 2 --callers HaplotypeCaller,bcftools,freebayes --min_AF 0.25 --min_coverage %i --min_chromosome_len 100000 --verbose --fraction_available_mem 0.25 --threads 12 --fractionRAM_to_dedicate 0.85 --outdir_callCNVs %s/call_CNVs"%(perSVade_wrapper, outdir_varcall_p2, genome, sbam_file, repeats_file, min_coverage_pos, outdir_varcall_general), "perSVade_env", 12))
            all_vars_generated = False
            continue

        # keep vcfs
        small_vars_vcfs += ["%s/call_small_variants/variant_calling_ploidy1.tab"%outdir_varcall_general, "%s/variant_calling_ploidy2.tab"%outdir_varcall_p2]

    # debug
    if all_vars_generated is False: continue

    # integrate variants using perSVade's module integrate_several_samples. I will run with table_comparisons_file for testing, but this has no further usage
    print("integrating variants...")
    variant_integration_dir = "%s/variant_integration"%species_dir; fun.make_folder(variant_integration_dir)
    variant_paths = "%s/paths_variant_calls.tab"%variant_integration_dir
    variant_paths_df = pd.DataFrame({srr : {"sampleID":srr, "call_small_variants_p2_outdir":"%s/small_var_call_p2"%(srr_to_varcallDir[srr])} for srr in variant_calling_srrs}).transpose().reset_index(drop=True)
    variant_paths_df["call_small_variants_p1_outdir"] = variant_paths_df.sampleID.apply(lambda s:  "%s/general_varcall/%s"%(srr_to_varcallDir[s], "call_small_variants"))
    variant_paths_df["sorted_bam"] = variant_paths_df.sampleID.apply(lambda s:  "%s/aligned_reads.bam.sorted"%(srr_to_alignReadsDir[s]))
    for m in ["integrate_SV_CNV_calls"]: variant_paths_df["%s_outdir"%m] = variant_paths_df.sampleID.apply(lambda s:  "%s/general_varcall/%s"%(srr_to_varcallDir[s], m))
    fun.save_df_as_tab(variant_paths_df, variant_paths)

    table_comparisons_file = "%s/sample_comparisons.tab"%variant_integration_dir
    fun.save_df_as_tab(pd.DataFrame({srr : {"sampleID":srr, "background_sampleIDs":",".join(variant_calling_srrs[-2:])} for srr in variant_calling_srrs[0:5]}).transpose(), table_comparisons_file)

    outdir_raw_integration = "%s/raw_integration"%variant_integration_dir
    fields_small_variants = ["#Uploaded_variation", "#CHROM", "POS", "REF", "ALT", "QUAL", "CALLEDALGS", "HC_GT", "INREPEATS", "ISSNP", "NCALLED", "NPASS", "PASSALGS", "bt_GT", "common_GT", "fb_GT", "mean_AD", "mean_DP", "mean_fractionReadsCov_PASS_algs", "mean_fractionReadsCov_called_algs", "relative_CN"]
    fields_small_variants_file = "%s/fields_small_variants.txt"%variant_integration_dir
    open(fields_small_variants_file, "w").write("\n".join(fields_small_variants))

    fields_SV_CNVs = ["#CHROM", "POS", "ID", "REF", "ALT", "INFO_BPS_TYPE", "INFO_BREAKENDIDs", "INFO_BREAKEND_FILTER", "INFO_BREAKEND_QUAL", "INFO_BREAKEND_coordinates", "INFO_BREAKEND_has_poly16GC", "INFO_BREAKEND_len_inserted_sequence", "INFO_BREAKEND_length_inexactHomology", "INFO_BREAKEND_length_microHomology", "INFO_BREAKEND_overlaps_repeats", "INFO_BREAKEND_real_AF", "INFO_BREAKPOINTID", "INFO_BREAKPOINTIDs", "INFO_END", "INFO_FILTER", "INFO_QUAL", "INFO_QUAL_max", "INFO_QUAL_mean", "INFO_QUAL_min", "INFO_RELCOVERAGE", "INFO_RELCOVERAGE_TO_3", "INFO_RELCOVERAGE_TO_5", "INFO_SVTYPE", "INFO_all_FILTERs", "INFO_any_has_poly16GC", "INFO_any_overlaps_repeats", "INFO_best_FILTER", "INFO_bpIDs", "INFO_has_poly16GC", "INFO_len_inserted_sequence_max", "INFO_len_inserted_sequence_mean", "INFO_len_inserted_sequence_min", "INFO_length_event_max", "INFO_length_event_mean", "INFO_length_event_min", "INFO_length_inexactHomology", "INFO_length_inexactHomology_max", "INFO_length_inexactHomology_mean", "INFO_length_inexactHomology_min", "INFO_length_microHomology", "INFO_length_microHomology_max", "INFO_length_microHomology_mean", "INFO_length_microHomology_min", "INFO_median_coverage_corrected", "INFO_median_relative_CN_AneuFinder", "INFO_median_relative_CN_HMMcopy", "INFO_merged_relative_CN", "INFO_misc", "INFO_overlaps_repeats", "INFO_real_AF", "INFO_real_AF_max", "INFO_real_AF_mean", "INFO_real_AF_min", "INFO_variantID", "INFO_worse_FILTER"]
    fields_SV_CNVs_file = "%s/fields_SV_CNVs.txt"%variant_integration_dir
    open(fields_SV_CNVs_file, "w").write("\n".join(fields_SV_CNVs))

    if fun.file_is_empty("%s/perSVade_finished_file.txt"%outdir_raw_integration): fun.run_cmd("%s integrate_several_samples --fraction_available_mem 1.0 --fractionRAM_to_dedicate 0.5 --min_chromosome_len 100000 --threads %i --verbose --paths_table %s  -o %s --fields_small_variants %s --fields_SV_CNVs %s -r %s -mchr %s --repeats_file %s -p %i"%(perSVade_wrapper, threads, variant_paths, outdir_raw_integration, fields_small_variants_file, fields_SV_CNVs_file, genome, mitochondrial_chromosome, repeats_file, ploidy), env="perSVade_env") # missing --table_comparisons %s table_comparisons_file

    # run small variant annotation (all in one)
    outdir_small_var_annot = "%s/small_var_annotation"%(species_dir)
    if fun.file_is_empty("%s/perSVade_finished_file.txt"%outdir_small_var_annot): fun.run_cmd("%s annotate_small_vars -o %s -r %s -mchr %s --merged_vcf %s -gff %s --mitochondrial_code %i --gDNA_code %i --min_chromosome_len 100000 --verbose --fraction_available_mem 0.1 --threads %i --fractionRAM_to_dedicate 0.85"%(perSVade_wrapper, outdir_small_var_annot, genome, mitochondrial_chromosome, "%s/integrated_small_variants_annotation_fields.vcf"%outdir_raw_integration, gff, mitochondrial_code, gDNA_code, threads), env="perSVade_env")

    # run SV_CNV annotation (all in one)
    outdir_SV_CNV_annot = "%s/SV_CNV_annotation"%(species_dir)
    if fun.file_is_empty("%s/perSVade_finished_file.txt"%outdir_SV_CNV_annot): fun.run_cmd("%s annotate_SVs -o %s -r %s -mchr %s --SV_CNV_vcf	 %s -gff %s --mitochondrial_code %i --gDNA_code %i --min_chromosome_len 100000 --verbose --fraction_available_mem 0.1 --threads %i --fractionRAM_to_dedicate 0.85"%(perSVade_wrapper, outdir_SV_CNV_annot, genome, mitochondrial_chromosome, "%s/integrated_SVs_CNVs_annotation_fields.vcf"%outdir_raw_integration, gff, mitochondrial_code, gDNA_code, threads), env="perSVade_env")

    # get filtered vars
    print("filtering vars...")
    small_vars_filt_file = "%s/smallVars_filt.py"%variant_integration_dir
    SV_CNV_filt_file = "%s/SV_CNV_filt.py"%variant_integration_dir

    if fun.file_is_empty(small_vars_filt_file):
        small_vars = fun.get_tab_as_df_or_empty_df("%s/integrated_small_variants.tab"%outdir_raw_integration)
        small_vars_annot = fun.get_tab_as_df_or_empty_df("%s/annotated_variants_corrrectedGene.tab"%outdir_small_var_annot)
        fun.get_filtered_small_vars_df(small_vars, small_vars_annot, small_vars_filt_file, replace=False)
        del small_vars; del small_vars_annot

    if fun.file_is_empty(SV_CNV_filt_file):
        SV_CNV = fun.get_tab_as_df_or_empty_df("%s/integrated_SVs_CNVs.tab"%outdir_raw_integration)
        SV_CNV_annot = fun.get_tab_as_df_or_empty_df("%s/annotated_variants_corrrectedGene.tab"%outdir_SV_CNV_annot)
        if ploidy!=1: raise ValueError("p should be 1")

        min_minAF = 0.2
        min_maxAF = 0.8
        max_relative_CN_deletion = 0.0
        min_relative_CN_duplication = 2.0
        max_relative_coverage_deletion = 0.1
        min_relative_coverage_duplication = 1.7

        fun.get_filtered_SV_CNV_df(SV_CNV, SV_CNV_annot, SV_CNV_filt_file, replace=False, min_minAF=min_minAF, min_maxAF=min_maxAF, min_lenCNV=600, max_relative_CN_deletion=max_relative_CN_deletion, min_relative_CN_duplication=min_relative_CN_duplication, max_relative_coverage_deletion=max_relative_coverage_deletion, min_relative_coverage_duplication=min_relative_coverage_duplication)
        del SV_CNV; del SV_CNV_annot


    # gwas cmds
    for drug in sorted(drug_to_interesting_srr_runs):
        #if sciName!="Candida_glabrata" or drug!="VRC": continue

        # define dirs
        outdir_gwas = "%s/GWAS_results/%s"%(species_dir, drug);
        os.makedirs(outdir_gwas, exist_ok=True)

        # load the df interpro annotations
        df_interpro = fun.load_InterProAnnotation("%s/InterproScan_annotation/interproscan_annotation.out"%taxID_dir)

        # define the gene features df for all species
        gene_features_df = fun.get_gene_features_df_all_species(DataDir, fun.get_species_to_gff(CurDir), ProcessedDataDir, replace=False)
        gene_features_df = gene_features_df[gene_features_df.species==sciName]

        # get cmds
        list_cmds += fun.get_GWAS_cmds_one_drug_validation_genes(outdir_gwas, drug, df_metadata_new_samples_spp, genome, df_GWAS_hits[df_GWAS_hits.drug==drug], ploidy, srr_to_varcallDir, min_coverage_pos, variant_calling_srrs, sciName, gff, threads, species_dir, df_interpro, gene_features_df)

# run in cluster
if len(list_cmds)>0:
    print("Running jobs in cluster")

    # get as a df
    df_cmds = pd.DataFrame(list_cmds, columns=["name_job", "cmd", "env", "threads"])

    # define jobs dir
    jobsdir_all = "%s/running_jobs_gwas_validation"%CurDir; fun.make_folder(jobsdir_all)

    # run one command for each name job
    for name_job in sorted(set(df_cmds.name_job)):
        print(name_job)

        # get df
        df_c = df_cmds[df_cmds.name_job==name_job]
        if len(set(df_c.env))!=1: raise ValueError("env should be the same")
        if len(set(df_c.threads))!=1: raise ValueError("threads should be the same")
        env = df_c.env.iloc[0]
        threads_per_job = df_c.threads.iloc[0]

        # define dirs
        jobsdir = "%s/%s"%(jobsdir_all, name_job); fun.make_folder(jobsdir)
        stddir = "%s/STDfiles"%jobsdir; fun.delete_folder(stddir); fun.make_folder(stddir)

        # define the jobs filename
        jobs_filename = "%s/jobs.%s"%(jobsdir, name_job)
        for f in [x for x in os.listdir(jobsdir) if x.startswith(fun.get_file(jobs_filename))]: fun.remove_file("%s/%s"%(jobsdir, f))
        open(jobs_filename, "w").write("\n".join(["source /gpfs/projects/bsc40/mschikora/anaconda3/etc/profile.d/conda.sh && conda activate %s > /dev/null && %s > %s/job.%i.std 2>&1"%(env, cmd, stddir, (I + 1)) for I, cmd in enumerate(df_c.cmd)]))
       
        # submit jobs
        os.chdir(jobsdir)
        fun.run_jobarray_file_MN4_greasy(jobs_filename, name_job, time="02:00:00", queue="debug", threads_per_job=threads_per_job, nodes=16) # max 16 nodes
        #fun.run_jobarray_file_MN4_greasy(jobs_filename, name_job, time="48:00:00", queue="bsc_ls", threads_per_job=threads_per_job, nodes=2)

    print("exiting execution since some jobs were submitted...")
    sys.exit(0)

#######################################







print("SUCCESS: Validation GWAS pipeline performed")