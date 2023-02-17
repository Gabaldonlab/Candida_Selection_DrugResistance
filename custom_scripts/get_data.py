#!/usr/bin/env python

# This script is useful to download all the data from candidamine in an optimised way. The idea is to download all the 'srr' files with prefetch in local and then with fastq-dump parallel in the cluster. There is info in the functions file about the procedence of all genomes. It runs jobs of perSVade v0.5 for the small variant calling and v0.6 (the exact same small variant calling) for the SV calling. It should be run with the Candida_mine_env

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
    
# define whether to submit jobs
submit_jobs = True

# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir
sys.path.insert(0, perSVade_dir)

# define dirs
CurDir = "%s/CandidaMine_data_generation/v1"%ParentDir
sys.path.insert(0, CurDir)
print("import functions")
import Cmine_functions as fun

# define threads
if run_in_cluster is True:

    cluster_name = fun.get_current_clusterName_mareNostrum()
    if cluster_name=="MN4": threads = 48
    elif cluster_name=="Nord3": threads=4


DataDir = "%s/data"%CurDir
PlotsDir = "%s/plots"%CurDir
JobsDir = "%s/running_jobs"%CurDir; fun.make_folder(JobsDir)
ProcessedDataDir = "%s/processed_data"%CurDir; fun.make_folder(ProcessedDataDir)
manually_curated_data = "%s/manually_curated_data"%CurDir

# define paths
perSVade_py = "%s/perSVade.py"%perSVade_dir # I used version 0.5, which is equivalent to 1.02.7 but more buggy
get_trimmed_reads_for_srr_py = "%s/get_trimmed_reads_for_srr.py"%perSVade_dir

####################################

# define the interesting taxIDs
interesting_taxIDs = [5476, 5480, 273371, 5482, 5478, 498019] # all (excluding metapsilosis (273372))
#interesting_taxIDs = [498019] # only auris

# get the taxIDs under these taxIDs

# initialize the cmds
all_cmds_getReads = []
all_cmds_bwamem = []
all_cmds_VarCall = []
all_cmds_treeRunning = []
all_cmds_treeRunningRandomHetSNPs= []
all_cmds_admixture = []
all_cmds_pairwiseSNPdistance = []
all_cmds_fastGEAR = []
all_cmds_fastGEAR_predefinedClades = []
all_cmds_fastGEAR_predefinedClades_RandomHetSNPs = []
all_cmds_pairwiseSNPdistance_fastGEAR_regions = []
all_cmds_hogwash_GWAS_drugResistance = []
all_cmds_ancestral_GWAS = []
all_cmds_ASRmutations = []

# define the resources
threads_getReads = 16
threads_bwamem = 16
threads_treeRunning = 16
threads_treeRunningRandomHetSNPs = 16 # on 16 it works well
threads_VarCall = 8; fraction_available_mem = 1.0 # nord3
#threads_VarCall = 48; fraction_available_mem = None; # MN4
threads_admixture = 16
threads_pairwiseSNPdistance = 8
threads_fastGEAR = 16
threads_fastGEAR_predefinedClades = 4 # with 8 would be enough, we use 16 to use all the node
threads_fastGEAR_predefinedClades_RandomHetSNPs = 3 # 8 may be necessary. 4 is necessary in MN and nord3
threads_pairwiseSNPdistance_fastGEAR_regions = 8
threads_hogwash_GWAS_drugResistance = 2 # it was enough with 600 Mb on average
threads_ASRmutations = 48
threads_ancestral_GWAS = 12 # run_ancestral_GWAS_gwas


# define interesting files
allSRRs_coverage = "%s/allSRRs_coverage.tab"%DataDir

# define the type of running (edit manually)
# "get_bams", "plot_coverage", "delete_reads", "run_varcall", "save_small_varcall_into_folder", "run_varcall", "get_SNPs_trees", "get_filtered_varcalls", run_pop_structure, "run_validate_bam_file", "get_pairwise_SNPdistance_per_window", "run_interproscan", "run_fastGEAR", "get_SNPs_trees_resamplingHeteroSNPs", "run_fastGEAR_predefinedClades", "run_fastGEAR_predefinedClades_RandomHetSNPs", "get_pairwise_SNPdistance_per_window_fastGEAR_regions", "run_hogwash_GWAS_drugResistance", "clean", 'run_ASRmutations', run_ancestral_GWAS_muts, run_ancestral_GWAS_gwas, "run_ancestral_GWAS_gwas_resampled"
steps_to_run = {"run_ancestral_GWAS_gwas_resampled"}

# define the filters of the coverage
min_coverage = 40
min_pct_covered = 90
min_coverage_pos = 12

# go through each taxID
for taxID in interesting_taxIDs:

    #### GENERAL FILES ####

    # get a folder to store everything
    sciName = fun.taxID_to_sciName[taxID]
    taxID_dir = "%s/%s_%i"%(DataDir, sciName, taxID); fun.make_folder(taxID_dir)
    print(sciName)

    #if sciName!="Candida_orthopsilosis": continue

    # define parms
    mitochondrial_chromosome = fun.taxID_to_mtChromosome[taxID]
    mitochondrial_code = fun.taxID_to_mDNA_code[taxID]
    gDNA_code = fun.taxID_to_gDNA_code[taxID]
    ploidy = fun.taxID_to_ploidy[taxID]
    
    # define the files
    genome = "%s/genome.fasta"%taxID_dir
    gff = "%s/annotations.gff"%taxID_dir
    readsDir = "%s/reads"%taxID_dir; fun.make_folder(readsDir)
    varCall_dir = "%s/varCall_output"%taxID_dir; fun.make_folder(varCall_dir)

    # get all taxIDs under this one
    taxIDs = fun.taxID_to_taxIDs[taxID]

    # get the data
    sra_data_dir = "%s/obtainingSRAdata"%taxID_dir; fun.make_folder(sra_data_dir)
    fileprefix_sra = "%s/sra"%sra_data_dir
    all_SRA_runInfo_df = fun.get_allWGS_runInfo_fromSRA_forTaxIDs(fileprefix_sra, taxIDs, genome, replace=False, min_coverage=10).set_index("Run", drop=False)
    all_SRA_runInfo_df = all_SRA_runInfo_df[~all_SRA_runInfo_df.Run.isin(fun.srrs_to_remove)]

    # get the previous repeats (needs to be retested if generated twice)
    previous_repeats_table = fun.get_repeat_maskerDF_from_perSVade(genome, perSVade_dir, threads=threads, replace=False)

    # remove SRRs
    for srr in fun.srrs_to_remove: fun.delete_folder("%s/%s"%(readsDir, srr))

    # define mappings
    srr_to_outdirReads = {srr : "%s/%s"%(readsDir, srr) for srr in all_SRA_runInfo_df.Run}
    srr_to_SRRfile = {srr : "%s/%s.srr"%(srr_to_outdirReads[srr], srr) for srr in all_SRA_runInfo_df.Run}
    srr_to_trimmedReads = {srr : tuple(["%s/%s_trimmed_reads_%i.fastq.gz"%(srr_to_outdirReads[srr], srr, N) for N in [1, 2]])  for srr in all_SRA_runInfo_df.Run}
    srr_to_sorted_bam = {srr : "%s/%s/aligned_reads.bam.sorted"%(varCall_dir, srr) for srr in all_SRA_runInfo_df.Run}
    srr_to_varCallOutdir =  {srr : "%s/%s"%(varCall_dir, srr) for srr in all_SRA_runInfo_df.Run}

    # define the srrs for which prefetch sould be done. These are reads where none of the SRRfile, trimmed reads or bam files are obtained
    srr_to_expectedFilesToNotRunPrefetch = {srr : {srr_to_SRRfile[srr], srr_to_trimmedReads[srr][1], srr_to_sorted_bam[srr]}  for srr in srr_to_outdirReads}
    srrs_withNeedPrefetch = [srr for srr, expected_files in srr_to_expectedFilesToNotRunPrefetch.items() if all([fun.file_is_empty(f) for f in expected_files])] 

    # define the SRRs for which to run fastqdump and trimming of reads
    srr_to_expectedFilesToNotRunFastqdump = {srr : {srr_to_trimmedReads[srr][1], srr_to_sorted_bam[srr]}  for srr in srr_to_outdirReads}
    srrs_withNeedFastqdump = [srr for srr, expected_files in srr_to_expectedFilesToNotRunFastqdump.items() if all([fun.file_is_empty(f) for f in expected_files])] 

    # define the SRRs for which you need to run the bam file obtention
    srrs_withNeedBWAmem = [srr for srr, sorted_bam in srr_to_sorted_bam.items() if fun.file_is_empty(sorted_bam)] 


    # there are two species, C. parapsilosis and tropicalis, for which the 'annotations.gff' has no entries for the mtDNA. I will run augustus on these to get the annotations.
    if sciName in {"Candida_parapsilosis", "Candida_tropicalis"}: gff = fun.get_gff_with_mtDNA_annotations_throughAugustus(gff, taxID, genome, replace=False, threads=threads)

    # check that the gff has annotations
    if fun.get_gff_has_mtDNA_annotations(gff, taxID) is False: raise ValueError("%s lacks mtDNA annotations"%gff)



    ###############################

    ########## READ OBTENTION AND ALIGNMENT ##########
    if "get_bams" in steps_to_run:

        # prefetch (local)
        if len(srrs_withNeedPrefetch)>0:
            print("There are %i/%i srrs to run prefetch on"%(len(srrs_withNeedPrefetch), len(srr_to_sorted_bam)))            
            if run_in_cluster is True: raise ValueError("You should not be running in the cluster")

            threads_available = multiproc.cpu_count()

            # define the inputs of the downloads
            inputs_downloads = [(srr, srr_to_SRRfile[srr], False) for srr in srrs_withNeedPrefetch]
            with multiproc.Pool(threads_available) as pool:
                list_srr_files = pool.starmap(fun.download_srr_with_prefetch, inputs_downloads)
                pool.close()

            continue

        # fastqdump (cluster)
        if len(srrs_withNeedFastqdump)>0:
            print("There are %i/%i srrs to run fastqdump and trimmomatic on"%(len(srrs_withNeedFastqdump), len(srr_to_sorted_bam)))

            for srr in srrs_withNeedFastqdump: all_cmds_getReads.append("%s --srr %s --outdir %s --threads %i"%(get_trimmed_reads_for_srr_py, srr, srr_to_outdirReads[srr], threads_getReads))

            continue

        # bwa mem
        if len(srrs_withNeedBWAmem)>0:
            print("There are %i/%i srrs to run bwa mem on"%(len(srrs_withNeedBWAmem), len(srr_to_sorted_bam)))

            for srr in srrs_withNeedBWAmem: all_cmds_bwamem.append("%s --ref %s --threads %i -o %s -f1 %s -f2 %s --mitochondrial_chromosome %s --StopAfter_bamFileObtention --skip_repeat_analysis"%(perSVade_py, genome, threads_bwamem, srr_to_varCallOutdir[srr],  srr_to_trimmedReads[srr][0],  srr_to_trimmedReads[srr][1], mitochondrial_chromosome))

            continue

    ##################################################

    # all the below steps require the bam to be calculated, so check
    if len(srrs_withNeedBWAmem)!=0: raise ValueError("You should have all the bams calculated")

    # get coverage df
    filename = "%s/coverage_all_srrs.tab"%taxID_dir
    coverage_df = fun.get_coverage_df_several_srrs(srr_to_sorted_bam, filename, taxID, sciName, replace=False)

    # remove bad srrs
    coverage_df = coverage_df[~(coverage_df.srr.isin(fun.srrs_to_remove))]

    ###### PLOT COVERAGE ########
    if "plot_coverage" in steps_to_run: fun.plot_coverage_allSRRs(coverage_df, "%s/%s_coverage.pdf"%(PlotsDir, sciName),  min_coverage=min_coverage, min_pct_covered=min_pct_covered)
    #############################

    #############################

    ####### DELETE READS #######
    if "delete_reads" in steps_to_run: fun.delete_folder(readsDir)
    ############################

    # keep only the SRRs hthat have a good coverage
    coverage_df_filt = coverage_df[(coverage_df.mean_coverage>=min_coverage) & (coverage_df.pct_covered>=min_pct_covered)]
    #print("There are %i/%i srrs that have a correct coverage"%(len(coverage_df_filt), len(coverage_df)))

    ####### RUN VARIANT CALLING ##########

    # this includes small variant calling (also for diploid for haploid organisms), CNV and SV calling
    if "run_varcall" in steps_to_run:
        print("running variant calling")

        # go through each srr
        for Isrr, srr in enumerate(coverage_df_filt.srr):
            perSVade_outdir = srr_to_varCallOutdir[srr]

            # define the final file
            final_file = "%s/perSVade_finished_file.txt"%perSVade_outdir
            #fun.remove_file(final_file) # debug (this is slow)

            # remove files to replace some species
            #if sciName in {"Candida_tropicalis", "Candida_parapsilosis"}: fun.remove_file(final_file)

            #  keep to run
            if fun.file_is_empty(final_file): 

                varcall_cmd = "%s --ref %s --threads %i -o %s -sbam %s  --mitochondrial_chromosome %s -gff %s --run_smallVarsCNV --caller all --coverage %i --mitochondrial_code %i --gDNA_code %i --ploidy %i --remove_smallVarsCNV_nonEssentialFiles --consider_repeats_smallVarCall --previous_repeats_table %s --run_ploidy2_ifHaploid --min_chromosome_len 100000 --verbose --nvars 50 --nsimulations 2 --simulation_ploidies auto --range_filtering_benchmark theoretically_meaningful_NoFilterRepeats --min_CNVsize_coverageBased 600 --window_size_CNVcalling 300 --cnv_calling_algs HMMcopy,AneuFinder --fractionRAM_to_dedicate 0.75 --replace_var_annotation"%(perSVade_py, genome, threads_VarCall, perSVade_outdir, srr_to_sorted_bam[srr], mitochondrial_chromosome, gff, min_coverage_pos, mitochondrial_code, gDNA_code, ploidy, previous_repeats_table)

                if fraction_available_mem is not None: varcall_cmd += " --fraction_available_mem %.2f"%fraction_available_mem

                all_cmds_VarCall.append(varcall_cmd)

                # --fraction_available_mem 0.5 may help
                # in Nord3, fraction_available_mem 1.0 makes total sense
                # --StopAfter_smallVarCall, --skip_SVcalling (only run small varcall)
                # --StopAfter_bamFileObtention to only run bamFileObtention
                # --replace_SV_CNVcalling 
                # --replace_var_annotation

    ######################################

    ######## RUN INTERPROSCAN ###########
    if "run_interproscan" in steps_to_run: fun.generate_interproscan_annotation_from_gff(gff, genome, "%s/InterproScan_annotation"%taxID_dir, taxID, threads=threads, replace=False)

    #####################################

    ######## GET THE SINGLE-FILE FILTERED VARIANT CALLING FILES ##########

    # you need to run ./integrate_varcalls.py (which integrates all the datasets) before this

    # define dirs
    integrated_varcallsDir = "%s/integrated_varcalls"%taxID_dir 
    small_vars_filt_file = "%s/smallVars_filt.py"%integrated_varcallsDir
    SV_CNV_filt_file = "%s/SV_CNV_filt.py"%integrated_varcallsDir

    # map each species to the samplesIDs
    srr_to_sampleID  = dict(fun.get_tab_as_df_or_empty_df("%s/samples_data.tab"%taxID_dir).set_index("srr")["sampleID"])

    # define a df with the sorted bams
    df_sorted_bams = pd.DataFrame({srr_to_sampleID[srr] : {"sampleID":srr_to_sampleID[srr], "sorted_bam":srr_to_sorted_bam[srr]} for srr in list(coverage_df_filt.srr)}).transpose()

    # this generates a single file with annotations and a single file with variants (all of them). It requires the run of integrate_varcalls.py (with the perSVade_env) before
    if "get_filtered_varcalls" in steps_to_run:
        print("getting filtered variant calling files")

        # get the filtered small variants
        if fun.file_is_empty(small_vars_filt_file):

            small_vars = fun.get_tab_as_df_or_empty_df("%s/smallVars.tab"%integrated_varcallsDir)
            small_vars_annot = fun.get_tab_as_df_or_empty_df("%s/smallVars_annot.tab"%integrated_varcallsDir)
            fun.get_filtered_small_vars_df(small_vars, small_vars_annot, small_vars_filt_file, replace=False)

            del small_vars
            del small_vars_annot

        # get filtered SVs
        if fun.file_is_empty(SV_CNV_filt_file):

            print("loading SV_CNV files")
            SV_CNV = fun.get_tab_as_df_or_empty_df("%s/SV_CNV.tab"%integrated_varcallsDir)
            SV_CNV_annot = fun.get_tab_as_df_or_empty_df("%s/SV_CNV_annot.tab"%integrated_varcallsDir)

            # define parameters depending on the ploidy
            if ploidy==1: # all samples are haploids according to the heterozygous SNPs analysis

                min_minAF = 0.2
                min_maxAF = 0.8
                max_relative_CN_deletion = 0.0
                min_relative_CN_duplication = 2.0
                max_relative_coverage_deletion = 0.1
                min_relative_coverage_duplication = 1.7

            elif ploidy==2:

                min_minAF = 0.1
                min_maxAF = 0.3
                max_relative_CN_deletion = 0.5
                min_relative_CN_duplication = 1.5
                max_relative_coverage_deletion = 0.6
                min_relative_coverage_duplication = 1.3

            else: raise ValueError("ploidy should be 1 or 2")

            SV_CNV_filt = fun.get_filtered_SV_CNV_df(SV_CNV, SV_CNV_annot, SV_CNV_filt_file, replace=False, min_minAF=min_minAF, min_maxAF=min_maxAF, min_lenCNV=600, max_relative_CN_deletion=max_relative_CN_deletion, min_relative_CN_duplication=min_relative_CN_duplication, max_relative_coverage_deletion=max_relative_coverage_deletion, min_relative_coverage_duplication=min_relative_coverage_duplication)

            del SV_CNV_filt
            del SV_CNV
            del SV_CNV_annot

    ######################################################################


    ########### GET A DF WITH THE SNPS THAT ARE IN UNIFORM POSITIONS ############

    # generate a df with the SNPs that are in positions with coverage >12x in all samples, have no heterozygous SNPs and have no IN/DELs
    haploid_snps_df_uniFormPositions_file = "%s/homoSNPs_positions_noINDELS_noHetSNPs_cov>%ix.py"%(integrated_varcallsDir, min_coverage_pos)

    outdir_generateHaploidSNPs_df = "%s/generating_haploid_snps_df_uniFormPositions_file"%integrated_varcallsDir
    fun.generate_haploid_snps_df_uniFormPositions_file(small_vars_filt_file, df_sorted_bams, haploid_snps_df_uniFormPositions_file, min_coverage_pos, outdir_generateHaploidSNPs_df, ploidy, threads=threads)

    fun.delete_folder(outdir_generateHaploidSNPs_df)

    # generate a fasta file with the haploid SNPs
    haploid_snps_multifasta = "%s/homoSNPs_positions_noINDELS_noHetSNPs_cov_atLeast_%ix_sequence.fasta"%(integrated_varcallsDir, min_coverage_pos)

    # define the sorted samples
    sorted_samples = sorted(set(df_sorted_bams.sampleID))

    if fun.file_is_empty(haploid_snps_multifasta):
        snps_df = fun.load_object(haploid_snps_df_uniFormPositions_file).set_index("sampleID", drop=False)
        fun.generate_multifasta_from_snps_df(snps_df, haploid_snps_multifasta, sorted_samples, threads=threads, generate_one_aln_each_chrom=True, pickRandomHetSNPs=False)

    # generate a concatenated alignment with only variable positions
    haploid_snps_multifasta_onlyVariableSites = "%s.onlyVariableSites.fasta"%haploid_snps_multifasta
    fun.get_multifasta_onlyVariableSites("%s.positions_df.py"%haploid_snps_multifasta, haploid_snps_multifasta_onlyVariableSites, sorted_samples, replace=False)

    # generate a dict mapping each chromosome to a multifasta like haploid_snps_multifasta_onlyVariableSites but for the same chromosome. This already only contains positions with variable positions
    all_possible_chroms = set(fun.get_chr_to_len(genome))
    def get_haploid_snps_multifasta_file(c): return "%s.%s.fasta"%(haploid_snps_multifasta, c)
    chrom_to_haploid_snps_multifasta = {c : get_haploid_snps_multifasta_file(c)  for c in all_possible_chroms if not fun.file_is_empty(get_haploid_snps_multifasta_file(c))}

    #############################################################################

    ######### SAVE ALL THE VARIANT CALLS INTO CandidaMine_v1_all_data/ #######

    if "save_small_varcall_into_folder" in steps_to_run:
        print("saving small variant calls into one folder in sepparate files")

        # define the dest_dir
        dest_dir = "%s/CandidaMine_v1_all_data"%DataDir
        fun.make_folder(dest_dir)

        # initialize the variant_annotation_df
        df_annotation_all = pd.DataFrame()
        variant_annotation_file = "%s/%s-variant_annotation.tab"%(dest_dir, sciName)

        # define all SRRS
        all_srrs_varCall = list(coverage_df_filt.srr)

        # go though each srr 
        for Isrr, srr in enumerate(all_srrs_varCall):

            # define the files
            if ploidy==1: origin_vcf = "%s/%s/smallVars_CNV_output/variants_atLeast2PASS_ploidy1.vcf"%(varCall_dir, srr)
            elif ploidy==2: origin_vcf = "%s/%s/smallVars_CNV_output/variants_atLeast2PASS_ploidy2.withMultiAlt.vcf"%(varCall_dir, srr)
            else: raise ValueError("%i is not a valid ploidy"%ploidy)

            dest_vcf = "%s/%s-%s-variants_ploidy%i.vcf"%(dest_dir, sciName, srr, ploidy)
            
            # check
            if fun.file_is_empty(origin_vcf): raise ValueError("%s is empty"%origin_vcf)

            # move 
            if fun.file_is_empty(dest_vcf):
                print("getting srr %i/%i"%(Isrr, len(all_srrs_varCall)))

                dest_vcf_tmp = "%s.tmp"%dest_vcf
                fun.run_cmd("rsync %s %s"%(origin_vcf, dest_vcf_tmp))
                os.rename(dest_vcf_tmp, dest_vcf)

            # load the df annotation
            if fun.file_is_empty(variant_annotation_file):
                print("getting annotations %i/%i"%(Isrr+1, len(all_srrs_varCall)))

                # load
                df_annotation = pd.read_csv("%s/%s/smallVars_CNV_output/variant_annotation_ploidy%i.tab"%(varCall_dir, srr, ploidy), sep="\t")

                # define the new vars
                if len(df_annotation_all)==0: current_vars = set()
                else: current_vars = set(df_annotation_all["#Uploaded_variation"])
                all_vars = set(df_annotation["#Uploaded_variation"])
                new_vars = all_vars.difference(current_vars)
                print("There are %i/%i new variants"%(len(new_vars), len(all_vars)))

                # get the df only with the new vars
                df_annotation = df_annotation[df_annotation["#Uploaded_variation"].isin(new_vars)]
                df_annotation_all = df_annotation_all.append(df_annotation)

        # save the variant annotation file
        if fun.file_is_empty(variant_annotation_file):
            print("saving variant annotation")

            variant_annotation_file_tmp = "%s.tmp"%variant_annotation_file
            df_annotation_all.to_csv(variant_annotation_file_tmp, sep="\t", index=False, header=True)
            os.rename(variant_annotation_file_tmp, variant_annotation_file)

        # get the genome and gff
        genome = "%s/genome.fasta"%taxID_dir
        dest_genome = "%s/%s-genome.fasta"%(dest_dir, sciName)
        if fun.file_is_empty(dest_genome):
            dest_genome_tmp = "%s.tmp"%dest_genome
            fun.run_cmd("rsync %s %s"%(genome, dest_genome_tmp))
            os.rename(dest_genome_tmp, dest_genome)

        dest_gff = "%s/%s-annotations.gff"%(dest_dir, sciName)
        if fun.file_is_empty(dest_gff):
            dest_gff_tmp = "%s.tmp"%dest_gff
            fun.run_cmd("rsync %s %s"%(gff, dest_gff_tmp))
            os.rename(dest_gff_tmp, dest_gff)

        # metadata table (only if existing. Probably just Candida glabrata)
        metadata_table = "%s/obtainingSRAdata/sra_metadata_manually_curated.tab"%(taxID_dir)
        if not fun.file_is_empty(metadata_table):
            dest_metadata_table = "%s/%s-manual_metadata.tab"%(dest_dir, sciName)
            if fun.file_is_empty(dest_metadata_table):
                dest_metadata_table_tmp = "%s.tmp"%dest_metadata_table
                fun.run_cmd("rsync %s %s"%(metadata_table, dest_metadata_table_tmp))
                os.rename(dest_metadata_table_tmp, dest_metadata_table)

    ##########################################################################


    ######## GENERATE THE PHYLOGENETIC TREES ###########

    if "get_SNPs_trees" in steps_to_run:
        print("getting tree for SNPs")

        outdir_tree = "%s/generate_tree_from_SNPs"%taxID_dir
        #fun.delete_folder(outdir_tree)
        fun.make_folder(outdir_tree)

        # get unrooted tree form the aln
        outfileprefix= "%s/iqtree_unroted"%outdir_tree
        final_file = "%s.iqtree"%outfileprefix
        if fun.file_is_empty(final_file):

            # remove previous files
            for f in os.listdir(outdir_tree):
                if f.startswith(fun.get_file(outfileprefix)): fun.remove_file("%s/%s"%(outdir_tree, f))

            # keep cmd
            cmd_tree = "iqtree -s '%s' -pre %s --mem 25G -m TEST+ASC -T AUTO -B 1000"%(haploid_snps_multifasta_onlyVariableSites, outfileprefix) # -m  GTR+ASC would do ascertainment bias correction (actually the ASC does). ASC is necessary if you input an alignment that does not span the whole genome (just some positions)
            all_cmds_treeRunning.append(cmd_tree)

    ####################################################


    ###### GENERATE PHYLOGENETIC TREE RESAMPLING ###########
    if "get_SNPs_trees_resamplingHeteroSNPs" in steps_to_run:
        print("getting jobs for get_SNPs_trees_resamplingHeteroSNPs")

        # this generates several 100x for different randomly picked heterozygous SNPs. This does not work in Nord3 because of the old GCLIB version

        # get a df with with diploid SNPs in positions that are covered and have no indels in all samples
        homo_and_hetero_snps_df_correctPositions_file = "%s/homo_and_hetero_SNPs_positions_noINDELS_cov_atLeast_%ix.py"%(integrated_varcallsDir, min_coverage_pos)

        outdir_SNPs = "%s/generating_homo_and_hetero_snps_df_correctPositions_file"%integrated_varcallsDir
        fun.generate_homo_and_hetero_snps_df_correctPositions(small_vars_filt_file, df_sorted_bams, homo_and_hetero_snps_df_correctPositions_file, min_coverage_pos, outdir_SNPs, threads=threads)
        fun.delete_folder(outdir_SNPs)

        # define a dir with all the trees
        outdir_tree_resampling = "%s/generate_tree_from_SNPs_resamplingHetSNPs"%taxID_dir; fun.make_folder(outdir_tree_resampling)
        outdir_resamplings = "%s/resamplings"%outdir_tree_resampling; fun.make_folder(outdir_resamplings)

        # make a file with all samples
        samples_file = "%s/samples.tab"%outdir_tree_resampling
        df_sorted_bams[["sampleID"]].to_csv(samples_file, sep="\t", index=False, header=True)

        # init whether all the resampled trees where genearetd
        all_resample_trees_generated = True
        all_resampled_trees = []

        # generate 100 bootstraps
        for I in range(1, 100+1):

            outdir_I = "%s/resample_%i"%(outdir_resamplings, I)
            #fun.delete_folder(outdir_I)
            final_file = "%s/tree_was_generated.txt"%outdir_I

            if fun.file_is_empty(final_file):

                cmd_treeResampling = "%s/get_tree_from_snps_df_resamplingHeteroSNPs.py --outdir %s --ref %s --threads %i --snps_df_file %s --samples_file %s"%(CurDir, outdir_I, genome, threads_treeRunningRandomHetSNPs, homo_and_hetero_snps_df_correctPositions_file, samples_file)
                all_cmds_treeRunningRandomHetSNPs.append(cmd_treeResampling)

                # record that some trees were not generated
                all_resample_trees_generated = False

            # keep the resampled tree
            all_resampled_trees.append("%s/iqtree_unroted.treefile"%outdir_I)

        # if all the resamples were generated, get the ML tree as the reference and use the bootstraps to get the support
        if all_resample_trees_generated is True: 

            print("generating consensus tree")
            #fun.generate_bestMLtree_withBootstrap_from_resampledTrees(all_resampled_trees, outdir_tree_resampling)
            fun.generate_consensus_withBootstrap_from_resampledTrees(all_resampled_trees, outdir_tree_resampling, replace=False)

    ########################################################

    ######## RUN FASTGEAR ############

    if "run_fastGEAR" in steps_to_run:
        print("running fastGEAR on all samples to infer recombination. Only for the homozygous positions")

        # get a multifasta with the positions with SNPs from haploid_snps_df_uniFormPositions_file
        outdir_fastGEAR = "%s/run_fastGEAR"%taxID_dir
        #fun.delete_folder(outdir_fastGEAR)
        fun.make_folder(outdir_fastGEAR)

        # run fastGEAR
        final_file_fastGEAR = "%s/fastGEAR_done.txt"%outdir_fastGEAR # this needs to be redefined
        if fun.file_is_empty(final_file_fastGEAR):

            # clean
            for f in ["output", "output_dir.mat", "output_snpData.mat"]: fun.delete_file_or_folder("%s/%s"%(outdir_fastGEAR, f))

            # prepare the SPECFILE for fastGEAR
            specFile = "%s/technical_specs.txt"%outdir_fastGEAR
            specFile_lines = ["15 # Number of iterations",
                              "%s # Upper bound for the number of clusters (possibly multiple values)"%(" ".join([str(x) for x in range(0, len(df_sorted_bams)+6, 5) if x>0])),
                              "1 # Run clustering for all upper bounds (0=no / 1=yes)",
                              "- # File containing a partition for strains",
                              "1 # 1=produce reduced output, 0=produce complete output"]

            open(specFile, "w").write("\n".join(specFile_lines)+"\n")

            # generate the fastGEAR command
            fastGEAR_cmd = "export LD_LIBRARY_PATH=/gpfs/projects/bsc40/mschikora/anaconda3/envs/fastGEAR_env/lib && %s/software/fastGEAR/fastGEARpackageLinux64bit/run_fastGEAR.sh '%s/software/fastGEAR/MCR_dir/v901' '%s' '%s/results.mat' '%s'"%(ParentDir, ParentDir, haploid_snps_multifasta, outdir_fastGEAR, specFile)
            all_cmds_fastGEAR.append(fastGEAR_cmd)

    ##################################

    ####### INFER POPULATION STRUCTURE #######
    if "run_pop_structure" in steps_to_run:

        # Run admixture to infer the population structure. This requires having the filtered sampleIDs

        print("running population structure")

        # define all samples (without the bad samples) and outlayers
        all_samples = set(srr_to_sampleID.values()).difference(fun.sciName_to_badSamples[sciName])
        outlayer_samples = fun.sciName_to_outlayerNumericSampleIDs[sciName]

        # check that outlayers exist
        missing_outlayers = outlayer_samples.difference(all_samples)
        if len(missing_outlayers)>0: raise ValueError("there are missing outlayers")

        for typePopStructureAnalysis, outdir_pop_structure in [("all_samples", "%s/pop_structure_inference_allSamples"%taxID_dir), ("noOutlayerSamples", "%s/pop_structure_inference"%taxID_dir)]:
            print(typePopStructureAnalysis)

            # debug the fat that there are no outlayers
            if len(outlayer_samples)==0 and typePopStructureAnalysis=="noOutlayerSamples": continue

            # define all the samples 
            if typePopStructureAnalysis=="all_samples": interesting_samples = all_samples
            elif typePopStructureAnalysis=="noOutlayerSamples": interesting_samples = all_samples.difference(outlayer_samples)
            else: raise ValueError("%s is not valid"%typePopStructureAnalysis)

            # define diploids and haploids
            if ploidy==2: 

                diploid_samples = interesting_samples
                haploid_samples = set()

            elif ploidy==1:

                diploid_samples = fun.sciName_to_diploidNumericSampleIDs[sciName]
                haploid_samples = interesting_samples.difference(diploid_samples)

            else: raise ValueError("You should define the haploid and diploid samples ")

            print("There are %i/%i diploid/haploid samples"%(len(diploid_samples), len(haploid_samples)))

            # get the jobs for admixture
            all_cmds_admixture += fun.get_admixture_jobs_severalKvals(small_vars_filt_file, diploid_samples, haploid_samples, outdir_pop_structure, df_sorted_bams, mitochondrial_chromosome, threads=threads, min_coverage_pos=min_coverage_pos, threads_admixture=threads_admixture)

    ##########################################


    ########## CALCULATE PAIRWISE SNP DISTANCES ##############

    if "get_pairwise_SNPdistance_per_window" in steps_to_run:

        # this should be ran in MN4 only
        
        # define all samples
        all_sorted_samples = sorted(df_sorted_bams.sampleID)

        # define jobs
        outdir_pairwiseSNP_distances = "%s/calculating_pairwiseSNP_distances"%taxID_dir; fun.make_folder(outdir_pairwiseSNP_distances)
        all_cmds_pairwiseSNPdistance_species = fun.get_pairwiseSNPdistances_per_window_jobs(small_vars_filt_file, outdir_pairwiseSNP_distances, mitochondrial_chromosome, df_sorted_bams, genome, ploidy, threads_pairwiseSNPdistance=threads_pairwiseSNPdistance, replace=False, min_coverage_pos=min_coverage_pos, threads=threads, window_size=10000, min_pct_covered=95)

        # add
        if len(all_cmds_pairwiseSNPdistance_species)>0: all_cmds_pairwiseSNPdistance += all_cmds_pairwiseSNPdistance_species
        
    ##########################################################

    ####### DEFINE THINGS THAT DEPEND ON THE TREES #########
    
    # define the metadata df (generated by analysis/descriptive_analysis.py, necessary for defining clades)
    if len({"run_fastGEAR_predefinedClades", "run_fastGEAR_predefinedClades_RandomHetSNPs", "run_hogwash_GWAS_drugResistance", "run_ancestral_GWAS_muts", "run_ancestral_GWAS_gwas", "run_ancestral_GWAS_gwas_resampled"}.intersection(steps_to_run))>0:
        print("getting files for metadata_df")

        # load the admixture df
        admixture_Qs_df, admixture_Kstats_df, Fst_divergence_df, admixture_Qs_se_df = fun.get_admixture_dfs(DataDir, ProcessedDataDir)

        # get the mapping between species and  numbers
        species_to_srr_to_sampleID = fun.get_species_to_srr_to_sampleID(CurDir)

        # map each species to the tree
        species_to_tree = fun.get_species_to_tree(CurDir)

        # get the distances per window
        df_pairwise_SNP_distances_perWindow_summaryStats = fun.get_df_pairwise_SNP_distances_perWindow_summaryStats(DataDir, ProcessedDataDir, species_to_tree, replace=False, threads=threads)

        # get the coverage per gene
        df_coverage_per_gene = fun.load_df_coverage_per_gene(CurDir, "%s/df_coverage_per_gene_all.py"%ProcessedDataDir , species_to_srr_to_sampleID, replace=False).set_index("species", drop=False)

        # get the metadata df (requires the previous running of ./get_metadata.py)
        metadata_df = pd.read_csv("%s/data/sra_metadata_manually_curated_all_species.tab"%CurDir, sep="\t")
        metadata_df = fun.get_df_metadata_with_resistance_and_others(metadata_df, manually_curated_data, species_to_srr_to_sampleID, DataDir, admixture_Qs_df, ProcessedDataDir, df_pairwise_SNP_distances_perWindow_summaryStats, species_to_tree, PlotsDir, df_coverage_per_gene, replace=False, make_plots=False, threads=threads)

        # keep only the metadata df for this sample
        metadata_df = metadata_df[metadata_df.species_name==sciName]
        if len(metadata_df)==0: raise ValueError("metadata_df can't be empty")

    ########################################################

    #### RUN FASTGEAR FOR PREDEFINED CLADES AND HOMOZYGOUS POSITIONS #######

    if "run_fastGEAR_predefinedClades" in steps_to_run:

        print("running fastGEAR with predefindec clades. This is using homoSNPs.")

        # define the files
        outdir_fastGEAR_predefinedClades = "%s/run_fastGEAR_predefinedClades_homoSNPs"%taxID_dir

        # get the cmd and add
        all_cmds_fastGEAR_predefinedClades += fun.get_cmds_fastGEAR_predefined_clades(outdir_fastGEAR_predefinedClades, chrom_to_haploid_snps_multifasta, metadata_df, sciName, genome, ploidy, threads, taxID_dir)

    ########################################################################


    ####### RUN FASTGEAR FOR PREDEFINED CLADES AND RANDOM SAMPLING OF HET. SAMPLES ########

    # this is similar to run_fastGEAR_predefinedClades, but on each of the resampled alignments of get_SNPs_trees_resamplingHeteroSNPs (all the alignments are the same)

    if "run_fastGEAR_predefinedClades_RandomHetSNPs" in steps_to_run:

        print("running fastGEAR on each of the randomly sampled hetero SNPs")

        # define the dir where everything will go
        outdir_fastGEAR_predefinedClades_randomHetSNPs = "%s/run_fastGEAR_predefinedClades_RandomHetSNPs"%taxID_dir
        fun.make_folder(outdir_fastGEAR_predefinedClades_randomHetSNPs)

        # define the inputs of the function that generates all the cmds
        inputs_fn = [(outdir_fastGEAR_predefinedClades_randomHetSNPs, I, taxID_dir, fun.cp.deepcopy(metadata_df), sciName, genome, ploidy, all_possible_chroms) for I in range(1, 100+1)]

        # get the cmd for the first resample (this is for the tree generation)
        x = inputs_fn [0]
        all_cmds_fastGEAR_predefinedClades_RandomHetSNPs += fun.get_cmds_run_fastGEAR_predefinedClades_RandomHetSNPs_one_resample(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7])

        # get cmds for the other resamples
        with  fun.multiproc.Pool(threads) as pool:

            all_cmds_fastGEAR_predefinedClades_RandomHetSNPs += fun.make_flat_listOflists(pool.starmap(fun.get_cmds_run_fastGEAR_predefinedClades_RandomHetSNPs_one_resample, inputs_fn[1:]))

            pool.close()
            pool.terminate()

    #######################################################################################




    ########## CALCULATE PAIRWISE SNP DISTANCES FASTGEAR REGIONS ##############

    if "get_pairwise_SNPdistance_per_window_fastGEAR_regions" in steps_to_run:

        # this should be ran in MN4 only
        
        # define all samples
        all_sorted_samples = sorted(df_sorted_bams.sampleID)

        # define outdir
        outdir_pairwiseSNP_distances_fastGEAR_regions = "%s/calculating_pairwiseSNP_distances_fastGEAR_regions"%taxID_dir; fun.make_folder(outdir_pairwiseSNP_distances_fastGEAR_regions)

        # get the windows in fastGEAR that are interesting (different for haploids and diploids)
        df_windows_fastGEAR = fun.get_df_windows_fastGEAR_any_recombination(taxID_dir, outdir_pairwiseSNP_distances_fastGEAR_regions, ploidy, min_coverage_pos, genome, replace=False, threads=threads)

        # define jobs
        all_cmds_pairwiseSNPdistance_fastGEAR_regions += fun.get_pairwiseSNPdistances_per_window_jobs(small_vars_filt_file, outdir_pairwiseSNP_distances_fastGEAR_regions, mitochondrial_chromosome, df_sorted_bams, genome, ploidy, threads_pairwiseSNPdistance=threads_pairwiseSNPdistance_fastGEAR_regions, replace=False, min_coverage_pos=min_coverage_pos, threads=threads, window_size=None, min_pct_covered=95, df_windows=df_windows_fastGEAR)

    ##########################################################################

    ######### RUN HOGWASH TO ASSOCIATE DRUG RESISTANCE AND MUTATIONS ###########

    if "run_hogwash_GWAS_drugResistance" in steps_to_run:

        # get the cmds to run hogwash for either 50,000 variants or all the variants related to a gene/pathway (We have typically been working with <500 samples and <500,000 genomic variants). This feeds from files generated in the analysis/descriptive_analysis.py

        raise ValueError("Hogwash is not working well, so that I should not use it.")

        # define dir
        outdir_hogwash = "%s/hogwash_GWAS_drugResistance"%taxID_dir; fun.make_folder(outdir_hogwash)

        # define the gene features df for all species
        gene_features_df = fun.get_gene_features_df_all_species(DataDir, fun.get_species_to_gff(CurDir), ProcessedDataDir, replace=False)
        gene_features_df = gene_features_df[gene_features_df.species==sciName]

        # define the mapping of reactome pathways to descriptions, obtained at 04/10/2021
        reactome_pathways_file =  "%s/annotation_files/ReactomePathways.txt"%(DataDir) 
        reactome_pathwayRelations_file = "%s/annotation_files/ReactomePathwaysRelation.txt"%DataDir

        # define the obo file to transfer GO terms
        obo_file = "%s/annotation_files/go-basic_30062021.obo"%(DataDir) # this was got from http://purl.obolibrary.org/obo/go/go-basic.obo

        # define the samples that are interesting for the GWAS
        species_to_drug_to_samplesForGWAS = fun.get_species_to_drug_to_samplesForGWAS_clinicalIsolates(metadata_df)
        if sciName in species_to_drug_to_samplesForGWAS: 
            
            # define the interesting samples    
            drug_to_samplesForGWAS = species_to_drug_to_samplesForGWAS[sciName]

            # get cmds
            cmds_hogwash_GWAS_drugResistance_species = fun.get_cmds_hogwash_GWAS_drugResistance(metadata_df, outdir_hogwash, taxID_dir, drug_to_samplesForGWAS, ploidy, sciName, gene_features_df, gff, obo_file, reactome_pathways_file, reactome_pathwayRelations_file, DataDir, threads=threads, threads_hogwash_GWAS_drugResistance=threads_hogwash_GWAS_drugResistance, replace=False)

            # integrate results or get single file
            if len(cmds_hogwash_GWAS_drugResistance_species)==0: pass
            else: all_cmds_hogwash_GWAS_drugResistance += cmds_hogwash_GWAS_drugResistance_species


    ############################################################################


    ####### ANCESTRAL GWAS ###########
    if len(steps_to_run.intersection({"run_ancestral_GWAS_muts", "run_ancestral_GWAS_gwas", "run_ancestral_GWAS_gwas_resampled"}))>0:

        # chose the step
        ancestral_GWAS_steps = steps_to_run.intersection({"run_ancestral_GWAS_muts", "run_ancestral_GWAS_gwas", "run_ancestral_GWAS_gwas_resampled"})
        if len(ancestral_GWAS_steps)!=1: raise ValueError("There can only be one")
        ancestral_GWAS_step = next(iter(ancestral_GWAS_steps))

        # hogwash did not work, so that I developed my own pipeline

        # define dir
        outdir_ancestralGWAS = "%s/ancestral_GWAS_drugResistance"%taxID_dir; fun.make_folder(outdir_ancestralGWAS)

        # define the gene features df for all species
        gene_features_df = fun.get_gene_features_df_all_species(DataDir, fun.get_species_to_gff(CurDir), ProcessedDataDir, replace=False)
        gene_features_df = gene_features_df[gene_features_df.species==sciName]


        # define the mapping of reactome pathways to descriptions, obtained at 04/10/2021
        reactome_pathways_file =  "%s/annotation_files/ReactomePathways.txt"%(DataDir) 
        reactome_pathwayRelations_file = "%s/annotation_files/ReactomePathwaysRelation.txt"%DataDir

        # define the obo file to transfer GO terms
        obo_file = "%s/annotation_files/go-basic_30062021.obo"%(DataDir) # this was got from http://purl.obolibrary.org/obo/go/go-basic.obo

        # define the samples that are interesting for the GWAS
        species_to_drug_to_samplesForGWAS = fun.get_species_to_drug_to_samplesForGWAS_clinicalIsolates(metadata_df)
        if sciName in species_to_drug_to_samplesForGWAS: 
            
            # define the interesting samples    
            drug_to_samplesForGWAS = species_to_drug_to_samplesForGWAS[sciName]

            # get the cmds for the normal GWAS 
            if ancestral_GWAS_step in {"ancestral_GWAS_muts", "run_ancestral_GWAS_gwas"}:

                all_cmds_ancestral_GWAS += fun.get_cmds_ancestral_GWAS_drugResistance(metadata_df, outdir_ancestralGWAS, taxID_dir, drug_to_samplesForGWAS, ploidy, sciName, gene_features_df, gff, obo_file, reactome_pathways_file, reactome_pathwayRelations_file, DataDir, genome, ancestral_GWAS_step, small_vars_filt_file, df_sorted_bams, min_coverage_pos, threads=threads, replace=False, threads_ancestral_GWAS=threads_ancestral_GWAS)

            elif ancestral_GWAS_step=="run_ancestral_GWAS_gwas_resampled": all_cmds_ancestral_GWAS += fun.get_cmds_ancestral_GWAS_drugResistance_resampledPhenotypes(sciName, outdir_ancestralGWAS, drug_to_samplesForGWAS, threads_ancestral_GWAS=threads_ancestral_GWAS, threads=threads)

            else: raise ValueError("invalid step")



    ##################################

    ######### RUN THE ASR OF EACH MUTATION ###########
    if "run_ASRmutations" in steps_to_run:
        print("run ASRmutations...")

        # define the dirs
        outdir_ASRmutations = "%s/generating_ASR_allMuts"%taxID_dir; fun.make_folder(outdir_ASRmutations)
        vars_df_file_asr = "%s/merged_vars_df.py"%outdir_ASRmutations
        integrated_file_ASRmutations = "%s/variants_df_with_ASRdata.py"%outdir_ASRmutations

        if fun.file_is_empty(integrated_file_ASRmutations):

            # generate the variants
            fun.get_vars_df_file_for_ASR_mutations(taxID_dir, vars_df_file_asr, taxID)

            # map each species to the tree
            tree_object = fun.get_species_to_tree(CurDir)[sciName]

            # get the cmd to run the ASR or integrate
            cmds_get_ASR = fun.get_all_ASR_mutations_cmds(tree_object, vars_df_file_asr, outdir_ASRmutations, threads_ASRmutations)

            if cmds_get_ASR is not None: all_cmds_ASRmutations += cmds_get_ASR
            else: fun.integrate_ASR_mutations(outdir_ASRmutations, integrated_file_ASRmutations, threads)

    ##################################################


    ########## CLEAN ##########

    # this is to clean various things
    if "clean" in steps_to_run:

        # define things to remove from taxID_dir
        things_to_remove = ["hogwash_GWAS_drugResistance", "run_fastGEAR", "annotations.gff_corrected.gff.bed_index1", "annotations.gff_corrected.gff.bed_index1.regions_lenght10000", "annotations.gff_corrected.gff_with_biotype.gff", "annotations.gff_corrected.gff_with_biotype.gff_clean.gff", "annotations.gff_corrected.gff_with_biotype.gff_clean.gz", "annotations.gff_corrected.gff_with_biotype.gff_clean.gz.tbi", "genome.fasta.windows100000bp.bed", "genome.fasta.windows100000bp.bed.generating.stderr", "genome.fasta.windows10000bp.bed", "genome.fasta.windows200000bp.bed", "genome.fasta.windows50000bp.bed", "genome.fasta.windows5000bp.bed"]

        # remove
        for f in things_to_remove: 
            print("remove %s"%f)
            fun.delete_file_or_folder("%s/%s"%(taxID_dir, f))

        # remove the reference_genome_dir
        for srr, varCallDir in srr_to_varCallOutdir.items(): 
            print("removing files from varCall_output %s"%srr)

            for f in ["reference_genome_dir", "smallVars_CNV_output/bcftools_ploidy1_out", "smallVars_CNV_output/bcftools_ploidy2_out", "smallVars_CNV_output/freebayes_ploidy1out", "smallVars_CNV_output/freebayes_ploidy2_out", "smallVars_CNV_output/HaplotypeCaller_ploidy1_out", "smallVars_CNV_output/HaplotypeCaller_ploidy2_out"]: fun.delete_file_or_folder("%s/%s"%(varCallDir, f))

    ###########################

#### RUN JOBS ####
inputs_sbatch = [("getReads", all_cmds_getReads, threads_getReads),
                 ("bwamem", all_cmds_bwamem, threads_bwamem),
                 ("varCall", all_cmds_VarCall, threads_VarCall),
                 ("treeSNPs", all_cmds_treeRunning, threads_treeRunning),
                 ("admixture", all_cmds_admixture, threads_admixture),
                 ("pairwiseSNPdistance", all_cmds_pairwiseSNPdistance, threads_pairwiseSNPdistance),
                 ("fastGEAR", all_cmds_fastGEAR, threads_fastGEAR),
                 ("treeSNPs_hetSNPs", all_cmds_treeRunningRandomHetSNPs, threads_treeRunningRandomHetSNPs),
                 ("fastGEAR_predefinedClades", all_cmds_fastGEAR_predefinedClades, threads_fastGEAR_predefinedClades),
                 ("fastGEAR_predefinedClades_hetSNPs", all_cmds_fastGEAR_predefinedClades_RandomHetSNPs, threads_fastGEAR_predefinedClades_RandomHetSNPs),
                 ("pairwiseSNPdistance_fastGEAR_regions", all_cmds_pairwiseSNPdistance_fastGEAR_regions, threads_pairwiseSNPdistance_fastGEAR_regions),
                 #("hogwash_GWAS_drugResistance", all_cmds_hogwash_GWAS_drugResistance, threads_hogwash_GWAS_drugResistance),

                 ("ancestral_GWAS_drugResistance", all_cmds_ancestral_GWAS, threads_ancestral_GWAS), # original
                 #("ancestral_GWAS_drugResistance_run_in_debug", all_cmds_ancestral_GWAS, threads_ancestral_GWAS), # debug_run

                 ("ASRmutations", all_cmds_ASRmutations, threads_ASRmutations)]



for typeJob, cmd_list_all, threads_per_job in inputs_sbatch:

    if len(cmd_list_all)>0:

        ############# GET THE JOBS FILES ##############

        njobs = len(cmd_list_all)
        print("running %s in parallel for %i jobs"%(typeJob, njobs))

        # define the jobs filename prefix
        dir_typeJob = "%s/%s"%(JobsDir, typeJob); fun.make_folder(dir_typeJob)
        jobs_filename_prefix = "%s/jobs.%s"%(dir_typeJob, typeJob)

        # remove all files with this prefix
        print("removing files with prefix")
        for f in os.listdir(dir_typeJob):
            file = "%s/%s"%(dir_typeJob, f) 
            if file.startswith(jobs_filename_prefix): fun.delete_file_or_folder(file)

        # init the list of job files
        jobs_filenames_list = []

        # define the STDdirs
        print("removing stddir")
        stddir = "%s/STDfiles"%dir_typeJob; fun.delete_folder(stddir); fun.make_folder(stddir)

        # define the number of chunks
        chunk_size = 5000
        nchunks = int(njobs/chunk_size) + 1

        # initialize the number of jobs
        currently_traversed_jobs = 0

        # go through chunks one million jobs
        for Ichunk, cmd_list in enumerate(fun.chunks(cmd_list_all, chunk_size)):
            print("working on chunk %i/%i"%(Ichunk, nchunks))

            # define the env dir activating cmd before each job
            typeJob_to_env = {"getReads":"perSVade_env", "bwamem":"perSVade_env", "varCall":"perSVade_env", "treeSNPs":"Candida_mine_env", "admixture":"Candida_mine_env", "pairwiseSNPdistance":"Candida_mine_env", "fastGEAR":"fastGEAR_env", "treeSNPs_hetSNPs":"Candida_mine_env", "fastGEAR_predefinedClades":"fastGEAR_env", "fastGEAR_predefinedClades_hetSNPs":"fastGEAR_env", "pairwiseSNPdistance_fastGEAR_regions":"Candida_mine_env", "hogwash_GWAS_drugResistance":"hogwash_env", "ancestral_GWAS_drugResistance":"ancestral_GWAS_env", "ASRmutations":"ancestral_GWAS_env", "ancestral_GWAS_drugResistance_only_CaurisCglabrata":"ancestral_GWAS_env", "ancestral_GWAS_drugResistance_run_in_debug":"ancestral_GWAS_env"}

            cmd_prefix = "source /gpfs/projects/bsc40/mschikora/anaconda3/etc/profile.d/conda.sh && conda activate %s"%(typeJob_to_env[typeJob])

            # define the prefix of the std
            std_perJob_prefix = "%s/%s"%(stddir, typeJob)

            # redefine the cmds to be redirected to sttdir
            def get_redefined_cmd(x):
                Icmd, c = x
                redefined_cmd = "%s && %s > %s.%i.out 2>&1"%(cmd_prefix, c, std_perJob_prefix, currently_traversed_jobs+Icmd+1)
                return redefined_cmd.replace(ParentDir, "/gpfs/projects/bsc40/mschikora")

            cmd_list = list(map(get_redefined_cmd, enumerate(cmd_list)))

            # update the currently_traversed_jobs
            currently_traversed_jobs += len(cmd_list)

            # define the jobs_filename for this chunk
            jobs_filename = "%s.chunk%i"%(jobs_filename_prefix, Ichunk+1)

            # write the jobs and keep
            jobs_filenames_list.append(jobs_filename)
            open(jobs_filename, "w").write("\n".join(cmd_list))

        #################################################

        # if you are in local, error
        if run_in_cluster is False: running_locally_error


        ##### generarate the run file ######

        print("generate the run file")
        if fun.cluster_name=="MN4": fun.run_jobarray_file_MN4_greasy_ListJobFiles(jobs_filenames_list, typeJob, time="02:00:00", queue="debug", threads_per_job=threads_per_job, nodes=16, submit=submit_jobs) # nodes=25 is the max. Each node will have 2GB RAM

        """
        Modify the other functions
        elif fun.cluster_name=="CTE_KNL": fun.run_jobarray_file_CTE_KNL_greasy(jobs_filename, typeJob, time="48:00:00", queue="bsc_ls", threads_per_job=threads_per_job, nodes=3, submit=submit_jobs) # 4 is the max (64 cores per node)


        elif fun.cluster_name=="Nord3": fun.run_jobarray_file_Nord3_greasy(jobs_filename, typeJob, time="48:00:00", queue="bsc_ls", threads_per_job=threads_per_job, RAM_per_thread=4000, nodes=32, submit=submit_jobs) # 32 nodes is the max. 32 gave problems with 4000 Mb/core (I put 16)


        else: raise ValueError("cluster_name should be MN4 or Nord3") 
        """

        ####################################

        print("You need to wait until the %s has been performed"%(typeJob))
        sys.exit(0)

############

# cdms to run in trantors

"""

# get the fastGEAR_predefinedClades_hetSNPs
./samba/scripts/run_MNjobArray_in_trantors.py --mn_jobsfile /data/mschikora/samba/CandidaMine_data_generation/v1/running_jobs/fastGEAR_predefinedClades_hetSNPs/jobs.fastGEAR_predefinedClades_hetSNPs --jobname fastGEAR_predefinedClades_hetSNPs --cpus_per_task 4 --walltime 100:00:00 --files_with_paths txt --prespecified_profile Cmine_get_data_fastGEAR_predefinedClades_hetSNPs --main_dir CandidaMine_data_generation --threads 64 

# (the cpus_per_task can be adapted to 3 for example)

# move the fastGEAR_predefinedClades_hetSNPs from the trantors to the cluster
./samba/scripts/transfer_data_from_trantor_to_MN.py  --prespecified_profile Cmine_get_data_fastGEAR_predefinedClades_hetSNPs --folders_to_sync CandidaMine_data_generation --threads 64

"""
