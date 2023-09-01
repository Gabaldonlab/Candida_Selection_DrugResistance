#!/usr/bin/env python

# This is a script to test in the MN tree reconstruction. It should be run with tree_from_SNPs_env. It will run the tree generation for various configurations of the C. glabrata public samples.

# env
import os
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

# define dirs
ProjectDir = "%s/scripts/tree_from_SNPs"%ParentDir
CurDir = "%s/testing/test_cglabrata"%ProjectDir
CmineDataDir = "%s/CandidaMine_data_generation/v1/data/Candida_glabrata_5478"%ParentDir
ref_genome = "%s/genome.fasta"%CmineDataDir
OutDir_all = "%s/outdir"%CurDir

# define parameters
mode = "haploid"

# prepare input files
data_dict = {}
for sampleID in ["SRR8068021", "SRR8068015", "SRR8953802", "SRR8697491"]:
    varcall_dir = "%s/varCall_output/%s/smallVars_CNV_output"%(CmineDataDir, sampleID)
    data_dict[sampleID] = {"sampleID":sampleID, "haploid_vcf":"%s/variants_atLeast2PASS_ploidy1.vcf"%varcall_dir, "diploid_vcf":"%s/variants_atLeast2PASS_ploidy2.vcf"%varcall_dir, "sorted_bam":"%s/varCall_output/%s/aligned_reads.bam.sorted"%(CmineDataDir, sampleID)}

df_paths = pd.DataFrame(data_dict).transpose().reset_index(drop=True)

# go through different modes
for mode, name_mode, vcf_fields in [("diploid_homozygous", "diploid_homozygous", ["diploid_vcf"]), ("haploid", "haploid_all", ["haploid_vcf", "diploid_vcf"]), ("haploid", "haploid_no_diploid", ["haploid_vcf"]),  ("diploid", "diploid_all",  ["diploid_vcf"])]:
    #if name_mode!="diploid_all": continue
    print(name_mode)

    # define outdir
    outdir = "%s/%s"%(OutDir_all, name_mode); os.makedirs(outdir, exist_ok=True)
    paths_table = "%s/paths.tbl"%outdir

    # get paths table
    df_paths[["sampleID", "sorted_bam"] + vcf_fields].sort_values(by="sampleID").to_csv(paths_table, sep="\t", index=False, header=True)

    # run
    print("running cmd")
    cmd = "%s/scripts/get_tree.py --paths_table %s -o %s/outdir_tree --mode %s --reference_genome %s --threads %i --min_coverage_pos 12"%(ProjectDir, paths_table, outdir, mode, ref_genome, threads)
    if mode=="diploid": cmd += " --n_resampled_trees 5"
    out_stat = os.system(cmd)
    if out_stat!=0: raise ValueError("error in run")





#testing/test_cglabrata/outdir