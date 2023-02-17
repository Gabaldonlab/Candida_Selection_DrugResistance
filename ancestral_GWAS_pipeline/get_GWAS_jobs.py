#!/usr/bin/env python

##### DEFINE ENVIRONMENT #######

# general module imports
import argparse, os, time, re, sys
from argparse import RawTextHelpFormatter
from ete3 import Tree
import pandas as pd
import multiprocessing as multiproc
import numpy as np

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import ancestral_GWAS_functions as fun

################################

###### ARGS ######

description = """
Inputs the output directory of get_ASR_mutations.py (--outdir_ASR_mutations), a table with the phenotypes (--phenotypes) and a (optional) table grouping genotypes (--grouping_table). It reconstructs the ancestral phenotypes with as specified by --phenotypes_pastml_prediction_method.

This script generates a list of jobs, each of which will output the association statistics of one mutation / group with the phenotype in a table (with the run_GWAS.py script). run_GWAS.py runs GWAS for one group of mutations for a binary phenotype. It runs the phyC and syncronous algorithms similarly to hogwash, but reimplemented. The definition of ancestral 1, 0 or NaN states is based on some methods (--ASR_methods_mutations, --ASR_methods_phenotypes) and relations between them (--classification_method_mutations, --classification_method_phenotypes). Once these states are defined for genotypes and phenotypes, we will run the two methods and write in a table the statistics of association for this group. 

These are the valid methods to pass to classification_method_<>:

- 'all': Require that all methods classify as 1 or 0 to classify a node. If this does not happen it is set to NaN. This is the most secure and restrictive.
- 'all_ignoreNaN': The same as 'all', but ignoring the methods that can't classify. This has the problem that you may have most methods saying NaN but one.

"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

# mandatory args
parser.add_argument("--outdir_ASR_mutations", dest="outdir_ASR_mutations", required=True, help="The path to the output directory given to get_ASR_mutations.py")
parser.add_argument("--phenotypes", dest="phenotypes", required=True, help="The table with the phenotypes. It should be a .tab file with a header and 'sampleID' and 'phenotype' (binary 1/0).")
parser.add_argument("--outdir", dest="outdir", required=True, help="The output directory, where all files will be written.")
parser.add_argument("--phenotypes_pastml_prediction_method", dest="phenotypes_pastml_prediction_method", required=True, help="The prediction method passed to pastml to reconstruct the ancestral phenotypes.")

parser.add_argument("--ASR_methods_mutations", dest="ASR_methods_mutations", required=True, help="The methods to be considered for the ASR of the mutations.")
parser.add_argument("--ASR_methods_phenotypes", dest="ASR_methods_phenotypes", required=True, help="The methods to be considered for the ASR of the phenotypes in the GWAS step.")
parser.add_argument("--classification_method_mutations", dest="classification_method_mutations", required=True, help="How to classify internal nodes into 1 and 0 as a function of --ASR_methods_mutations.")
parser.add_argument("--classification_method_phenotypes", dest="classification_method_phenotypes", required=True, help="How to classify internal nodes into 1 and 0 as a function of --ASR_methods_phenotypes.")
parser.add_argument("--min_support", dest="min_support", required=True, type=int, help="The branches with a support below this value won't be considered. This should be a number between 0 and 100.") 
parser.add_argument("--resampled_phenotypes_pastml_out", dest="resampled_phenotypes_pastml_out", required=True, type=str, help="The link to the resampled phenotypes") 


# optional args
parser.add_argument("--replace", dest="replace", action="store_true", default=False, help="Re-run all the steps by deleting the output directory.")
parser.add_argument("--grouping_table", dest="grouping_table", type=str, default=None, help="The path with the grouping table, with a field like --varID_field and a field called 'group'.")
parser.add_argument("--varID_field", dest="varID_field", type=str, default="variantID_across_samples", help="The field in variants_table that defines the unique mutations' IDs. By default it is 'variantID_across_samples'.")
parser.add_argument("--threads", dest="threads", type=int, default=4, help="The number of threads.")
parser.add_argument("--max_pval_keepASR", dest="max_pval_keepASR", default=0.05, type=float, help="The maximum raw p-value to consider and save the ASR results.")
parser.add_argument("--skip_tree_rendering", dest="skip_tree_rendering", action="store_true", default=False, help="Skip the rendering of the tree.")
parser.add_argument("--run_without_parallelization", dest="run_without_parallelization", action="store_true", default=False, help="Skip the parallelization")
parser.add_argument("--stop_after_resampled_phenotypes_state_df", dest="stop_after_resampled_phenotypes_state_df", action="store_true", default=False, help="Stop after the generation of the states of resampled phenotypes.")
parser.add_argument("--all_pval_methods", dest="all_pval_methods", default=None, type=str, help="A string of comma-sepparated pval methods to consider. It can be any of 'RelToBranchLen,phenotypes'")
parser.add_argument("--interesting_gwas_methods", dest="interesting_gwas_methods", default='phyC,synchronous', type=str, help="A string of comma-sepparated methods to consider")

# parse all the args
opt = parser.parse_args()

##################

print("running get_GWAS_jobs.py")
start_time = time.time()

# remove outdir if replace, and set replace to False
if opt.replace is True: 
    print("deleting output folder %s"%opt.outdir)
    fun.delete_folder(opt.outdir)

# check that the min_support is correct
if opt.min_support<0 or opt.min_support>100: raise ValueError("min_support should be between 0 and 100")

# define full paths
opt.outdir = fun.get_fullpath(opt.outdir)
opt.outdir_ASR_mutations = fun.get_fullpath(opt.outdir_ASR_mutations)
opt.phenotypes = fun.get_fullpath(opt.phenotypes)
if opt.grouping_table is not None: fun.get_fullpath(opt.grouping_table)

# modify the all_pval_methods if specified
if opt.all_pval_methods is not None:

    all_pval_methods = opt.all_pval_methods.split(",")
    for m in all_pval_methods:
        if m not in ["RelToBranchLen", "phenotypes"]: raise ValueError("invalid %s in --all_pval_methods"%m)

    fun.all_pval_methods = all_pval_methods

else: fun.all_pval_methods = ["RelToBranchLen", "phenotypes"]

# define the fields of the GWAS
fields_GWAS_stats_df = ['method', 'phenotype_f', 'genotype_f', 'nodes_withGeno', 'nodes_withoutGeno', 'nodes_unkownGeno', 'nodes_withPheno', 'nodes_withoutPheno', 'nodes_unkownPheno', 'nodes_GenoAndPheno', 'nodes_noGenoAndPheno', 'nodes_GenoAndNoPheno', 'nodes_noGenoAndNoPheno', 'epsilon', 'OR', 'chi_square', 'group_name', 'ASR_methods_phenotypes', 'ASR_methods_mutations', 'classification_method_phenotypes', 'classification_method_mutations']
pval_fields = []
for pval_method in all_pval_methods: pval_fields += ["%s_%s"%(x, pval_method) for x in ["pval_chi_square", "pval_OR", "pval_GenoAndPheno", "pval_epsilon"]]
fields_GWAS_stats_df += pval_fields

# pass to fun
fun.fields_GWAS_stats_df = fields_GWAS_stats_df
fun.pval_fields = pval_fields

# check the interesting_gwas_methods
opt.interesting_gwas_methods = opt.interesting_gwas_methods.split(",")
for m in opt.interesting_gwas_methods:
    if m not in ["phyC", "synchronous"]: raise ValueError("invalid %s in --interesting_gwas_methods"%m)

# make the outdir
fun.make_folder(opt.outdir)
fun.print_with_runtime("running GWAS on outdir '%s'"%opt.outdir)

# define final files of the integration, and break if they exist
integrated_ASR_file_raw = "%s/integrated_ASR_raw.tab"%opt.outdir
integrated_ASR_file_raw_tmp = "%s.tmp"%integrated_ASR_file_raw
integrated_ASR_file = "%s/integrated_ASR.tab"%opt.outdir

integrated_GWAS_stats_file_raw = "%s/integrated_GWAS_stats_raw.tab"%opt.outdir
integrated_GWAS_stats_file_raw_tmp = "%s.tmp"%integrated_GWAS_stats_file_raw
integrated_GWAS_stats_file = "%s/integrated_GWAS_stats.tab"%opt.outdir

if not fun.file_is_empty(integrated_GWAS_stats_file):
    print("The final file %s of this script is already existing, skip."%(integrated_GWAS_stats_file))
    sys.exit(0)

# load the df with the phenotypes
phenotypes_df = fun.get_tab_as_df_or_empty_df(opt.phenotypes)[["sampleID", "phenotype"]]
phenotypes_df["sampleID"] = phenotypes_df.sampleID.apply(str)
strange_phenotypes = set(phenotypes_df.phenotype).difference({0, 1})
if len(strange_phenotypes)>0: raise ValueError("There are strange phenotypes: %s"%strange_phenotypes)

# check that the pastml prediction methods are correct
allowed_predictionMethods = {'MPPA', 'MAP', 'JOINT', 'DOWNPASS', 'ACCTRAN', 'DELTRAN', 'ALL', 'ML', 'MP'}
if opt.phenotypes_pastml_prediction_method not in allowed_predictionMethods: raise ValueError("The --phenotypes_pastml_prediction_method '%s' is not correct. It should be among %s"%(opt.phenotypes_pastml_prediction_method, allowed_predictionMethods))

# redefine the inputs of run_GWAS.py
all_methods = ['MPPA', 'MAP', 'JOINT', 'DOWNPASS', 'ACCTRAN', 'DELTRAN']
if opt.ASR_methods_mutations=="ALL": opt.ASR_methods_mutations = ",".join(all_methods)
if opt.ASR_methods_phenotypes=="ALL": opt.ASR_methods_phenotypes = ",".join(all_methods)

# define the tree (as used in the mutations) and check that there are no politomies
correct_treefile = "%s/correct_tree.nw"%opt.outdir_ASR_mutations
fun.check_that_tree_has_no_politomies(Tree(correct_treefile))

# map each node name to the support, and leafs have 100 support
tree = fun.get_tree_with_internalNodeNames(Tree(correct_treefile)); tree.name = "root"
node_to_support = {n.name : n.support for n in tree.traverse()}
for l in tree.get_leaves(): node_to_support[l.name] = 100.0
node_to_support["root"] = 100.0
fun.check_that_tree_has_no_politomies(Tree(correct_treefile))

# build the ancestral state of the phenotypes
phenotypes_for_pastml = "%s/phenotypes.tab"%opt.outdir
outfile_pastml_phenotypes = "%s/phenotypes_ASR.out"%opt.outdir
outfile_pastml_phenotypes_html = "%s/phenotypes_ASR.html"%opt.outdir
workdir_pastml_phenotypes = "%s/phenotypes_ASR_workdir"%opt.outdir
outfile_pastml_phenotypes_std = "%s/phenotypes_ASR.std"%opt.outdir

if fun.file_is_empty(outfile_pastml_phenotypes):
    print("ASR of phenotypes")
    fun.save_df_as_tab(phenotypes_df, phenotypes_for_pastml)
    fun.run_cmd("%s %s %s %s %s %s %s > %s 2>&1"%(fun.RUN_PASTML_PY, correct_treefile, phenotypes_for_pastml, workdir_pastml_phenotypes, outfile_pastml_phenotypes, outfile_pastml_phenotypes_html, opt.phenotypes_pastml_prediction_method, outfile_pastml_phenotypes_std))

fun.remove_file(outfile_pastml_phenotypes_std)

# get the phenotype_states_df, which is a df used in the GWAS, the same for all mutations
phenotype_states_df_file = "%s/phenotype_states_df.py"%opt.outdir
if fun.file_is_empty(phenotype_states_df_file):
    print("Getting phenotype states...")

    phenotype_states_df = fun.get_ASR_states_integrate_methods(["phenotype"], [outfile_pastml_phenotypes], opt.ASR_methods_phenotypes, opt.classification_method_phenotypes, opt.min_support, node_to_support)

    fun.save_object(phenotype_states_df, phenotype_states_df_file)

# get a df with the resampled phenotype states
nresamples_to_filePhenotypesResampling = {}
for nresamples in [1000, 10000]:

    resampled_phenotype_states_df_file = "%s/phenotype_states_df_%iresamples.py"%(opt.outdir, nresamples)
    fun.generate_resampled_phenotype_states_df(opt.resampled_phenotypes_pastml_out, resampled_phenotype_states_df_file, nresamples, opt.ASR_methods_phenotypes, opt.classification_method_phenotypes, opt.threads, correct_treefile, opt.min_support, node_to_support)

    nresamples_to_filePhenotypesResampling[nresamples] = resampled_phenotype_states_df_file


# stop after this
if opt.stop_after_resampled_phenotypes_state_df is True:
    print("Exiting after the generate_resampled_phenotype_states_df function...")
    sys.exit(0)

#### CREATE GWAS JOBS #######

# define files for each 
ASR_mutations_dir = "%s/ASR_mutations"%opt.outdir_ASR_mutations

# define the available mutations
sorted_mutations = sorted(set(map(lambda x: x.split("_")[0], os.listdir(ASR_mutations_dir))))

# map each group to the mutations
print("mapping each group to the mutations. There are %i mutations"%len(sorted_mutations))
group_to_mutations_file = "%s/group_to_mutations_dict.py"%opt.outdir
if fun.file_is_empty(group_to_mutations_file):

    group_to_mutations = fun.get_GWAS_group_to_mutations(sorted_mutations, opt.grouping_table, opt.varID_field, opt.outdir_ASR_mutations)
    fun.save_object(group_to_mutations, group_to_mutations_file)

group_to_mutations = fun.load_object(group_to_mutations_file)

# check that group_to_mutations is correct
if set(group_to_mutations.keys())=={"group", "numeric_mutation"} and all([len(m)==0 for g,m in group_to_mutations.items()]): raise ValueError("The group_to_mutations is incorrect, suggesting that it is empty: %s"%group_to_mutations) # this means that it was empty

# create a df where each line is one job
print("generating the inputs of the parallel job running...")
outdir_GWAS = "%s/GWAS_results"%opt.outdir; fun.make_folder(outdir_GWAS)
df_cmds = pd.DataFrame({"mutations_list":group_to_mutations}).reset_index().rename(columns={"index":"groupName"}).set_index("groupName", drop=False)

# sort the mutations and add the representative group name
df_cmds["mutations_list"]  = df_cmds.mutations_list.apply(sorted)
df_cmds, df_rep_groups = fun.get_df_cmds_gwas_with_rep_group_name(df_cmds)

# define all the rep groups
all_representative_groups = set(df_rep_groups.rep_group)

# define groups that you still did not ran
print("checking which jobs are to be repeated...")
CurDir = fun.get_fullpath(os.getcwd())
os.chdir(outdir_GWAS)
df_cmds["outfile_group"] = df_cmds.groupName + "_GWAS_stats.tab"
df_cmds["outfile_is_empty"] = df_cmds.outfile_group.apply(fun.file_is_empty)
os.chdir(CurDir)
groups_still_to_run = set(df_cmds[df_cmds.outfile_is_empty].groupName)

# keep only cmds of the rep groups and the outfile that is empty
df_cmds = df_cmds[(df_cmds.groupName==df_cmds.rep_groupName) & (df_cmds.groupName.isin(groups_still_to_run))]

# keep only cmds where the outfile is empty
df_cmds = df_cmds[df_cmds.outfile_is_empty]

# if there are jobs, run them
if len(df_cmds)>0:
    print("preparing inputs of run_GWAS_as_function")

    # checks
    fun.check_that_tree_has_no_politomies(Tree(correct_treefile))

    # add fields
    ngroups = len(df_cmds)
    df_cmds["Igroup"] = list(range(1, ngroups+1))
    df_cmds["pct_progress"] = ((df_cmds.Igroup / ngroups)*100).apply(lambda x: "%.3f"%x)

    # define the inputs to parallelize
    def get_input_fn_from_r(r): return (outdir_GWAS, ASR_mutations_dir, phenotype_states_df_file, correct_treefile, r.groupName, r.mutations_list, opt.ASR_methods_mutations, opt.classification_method_mutations, opt.ASR_methods_phenotypes, opt.classification_method_phenotypes,  opt.max_pval_keepASR, r.pct_progress, nresamples_to_filePhenotypesResampling, opt.outdir_ASR_mutations, opt.min_support, node_to_support, opt.interesting_gwas_methods)
    inputs_fn = df_cmds.apply(get_input_fn_from_r, axis=1).values

    # define the number of chunks
    chunk_size = opt.threads*12
    nchunks = len(list(fun.chunks([1]*len(inputs_fn), chunk_size)))

    # go through chunks of thread's size
    for Ic, inputs_fn_c in enumerate(fun.chunks(inputs_fn, chunk_size)):
        fun.print_with_runtime("\n----------\nRunning %i NR GWAS jobs in %i threads for job chunk %i/%i (%i NR mutations)"%(ngroups, opt.threads, Ic+1, nchunks, len(inputs_fn))) 
        fun.print_available_mem()

        # run in parallel
        if opt.run_without_parallelization is False:

            with multiproc.Pool(opt.threads) as pool:
                pool.starmap(fun.run_GWAS_as_function, inputs_fn_c)
                pool.close()
                pool.terminate()
                pool.join()

        # run without parallelization
        else:
            for x in inputs_fn_c: fun.run_GWAS_as_function(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15])

# create the ASR_file and gwas_stats_file files for the redundant, non-created groups
df_rep_groups = df_rep_groups[(df_rep_groups.group.isin(groups_still_to_run)) & (df_rep_groups.group!=df_rep_groups.rep_group)]
if len(df_rep_groups)>0:
    fun.print_with_runtime("copying files for %i redundant groups. Running on %i threads"%(len(df_rep_groups), opt.threads))

    # add the progress
    ngroups = len(df_rep_groups)
    df_rep_groups["Igroup"] = list(range(1, ngroups+1))
    df_rep_groups["pct_progress"] = ((df_rep_groups.Igroup / ngroups)*100).apply(lambda x: "%.3f"%x)

    # define the inputs
    inputs_fn = list(df_rep_groups.apply(lambda r: (outdir_GWAS, r.rep_group, r.group, r.pct_progress, opt.max_pval_keepASR), axis=1).values)

    # run in parallel
    if opt.run_without_parallelization is False:

        with multiproc.Pool(opt.threads) as pool:
            pool.starmap(fun.copy_gwas_results_redundant_files, inputs_fn, chunksize=1)
            pool.close()
            pool.terminate()

    else:
        for x in inputs_fn: fun.copy_gwas_results_redundant_files(x[0], x[1], x[2], x[3], x[4])

else: print("There are no redundant groups")

# log
fun.print_with_runtime("All cmds for GWAS worked well. Integrating...")

#############################

########### INTEGRATE ###########

# check that all the GWAS results were ran
groups_ASR = fun.get_set_groups_done_gwas_outdir("_ASR.tab", outdir_GWAS)
groups_GWASstats = fun.get_set_groups_done_gwas_outdir("_GWAS_stats.tab", outdir_GWAS)
all_groups = set(group_to_mutations)
if groups_ASR!=all_groups: raise ValueError("groups_ASR has not all expected. Rerun")
if groups_GWASstats!=all_groups: raise ValueError("groups_GWASstats has not all expected. Rerun")

# integrate all ASR results
print("integrating ASR results")
if fun.file_is_empty(integrated_ASR_file_raw):
    fun.run_cmd("cat %s/*_ASR.tab > %s"%(outdir_GWAS, integrated_ASR_file_raw_tmp))
    os.rename(integrated_ASR_file_raw_tmp, integrated_ASR_file_raw)

print("loading asr_df")
asr_df = pd.read_csv(integrated_ASR_file_raw, sep="\t", header=None, names=fun.fields_ASR_df)
  
# integrate all the GWAS results
print("integrating GWAS results")
if fun.file_is_empty(integrated_GWAS_stats_file_raw):
    fun.run_cmd("cat %s/*_GWAS_stats.tab > %s"%(outdir_GWAS, integrated_GWAS_stats_file_raw_tmp))
    os.rename(integrated_GWAS_stats_file_raw_tmp, integrated_GWAS_stats_file_raw)

print("loading gwas_df")
gwas_df = pd.read_csv(integrated_GWAS_stats_file_raw, sep="\t", header=None, names=fun.fields_GWAS_stats_df)

# add the maxT calculation
if "phenotypes" in fun.all_pval_methods: gwas_df = fun.get_gwas_df_with_maxT_pvals(gwas_df, outdir_GWAS, all_representative_groups, opt.outdir)

# check that all the groups are correct
if len(opt.interesting_gwas_methods)==1:
    if len(gwas_df)!=len(set(gwas_df.group_name)): raise ValueError("group_name should be unique")

elif len(opt.interesting_gwas_methods)==2:
    if set(gwas_df.groupby("group_name").apply(len))!={2}: raise ValueError("There should be 2 rows for each group")

if set(gwas_df.group_name)!=all_groups: raise ValueError("Not all groups had a correct GWAS")
if set(asr_df.group_name)!=all_groups: raise ValueError("Not all groups had a correct ASR")
if any(gwas_df[fun.pval_fields].apply(pd.isna, axis=1).apply(any)): raise ValueError("There can't be NaNs in the pvalues")

# map each numeric variant to the actual variant
vars_IDmapping_df = fun.get_tab_as_df_or_empty_df("%s/variants_IDmapping.tab"%opt.outdir_ASR_mutations)
numericVar_to_var = dict(vars_IDmapping_df.set_index("numeric_mutation")[opt.varID_field])

# map each representative mutation to an array of all the mutations
repMuts_df = fun.get_tab_as_df_or_empty_df("%s/mutation_to_representative_mutation.tab"%opt.outdir_ASR_mutations)
repMut_to_mutations = repMuts_df.sort_values(by=["rep_mutation", "mutation"]).set_index("rep_mutation").mutation.groupby("rep_mutation").apply(np.array)
repMut_to_n_mutations = repMut_to_mutations.apply(len)

# add the 'variant_name' to the dfs if there is no grouping, and also add the non-representative mutations
if opt.grouping_table is None:
    asr_df = fun.get_asr_df_with_all_mutations(asr_df, repMut_to_n_mutations, repMut_to_mutations, numericVar_to_var, opt.threads, opt.run_without_parallelization)
    gwas_df = fun.get_gwas_df_with_all_mutations(gwas_df, repMut_to_n_mutations, repMut_to_mutations, numericVar_to_var, opt.interesting_gwas_methods)

# add the FDR-corrected pval
fun.print_with_runtime("adding FDR corrected p values")
gwas_df = fun.get_gwas_df_with_fdr_corrected_pvals(gwas_df)

# checks of the pval
if any(((gwas_df[fun.pval_fields]>=0).apply(any, axis=1))!=((gwas_df[fun.pval_fields]>=0).apply(all, axis=1))): raise ValueError("All the pvals should be -1 together")

# generate plots
gwas_df_plots = gwas_df[gwas_df[fun.pval_fields[0]]>=0]
plots_dir = "%s/plots"%opt.outdir
pval_fields_plot = pval_fields + [k for k in gwas_df.keys() if k.endswith("_maxT")]
fun.generate_pval_distributions(gwas_df_plots, plots_dir, pval_fields_plot)
fun.generate_qq_plot(gwas_df_plots, plots_dir, pval_fields_plot)

# keep the treefile into outdir
fun.check_that_tree_has_no_politomies(Tree(correct_treefile))
fun.rsync_file(correct_treefile,  "%s/%s"%(opt.outdir, correct_treefile.split("/")[-1]))

# clean
fun.print_with_runtime("cleaning...")
for f in [integrated_GWAS_stats_file_raw, integrated_GWAS_stats_file_raw_tmp, integrated_ASR_file_raw, integrated_ASR_file_raw_tmp]: fun.delete_file_or_folder(f)
for f in ["GWAS_results", "group_to_mutations_dict.py", "phenotypes_ASR.html", "phenotypes_ASR.out", "phenotype_states_df_1000resamples.py", "phenotype_states_df_10000resamples.py", "phenotype_states_df.py", 'gwas_df_w_maxTpvals.py', 'integrated_phenotype_resamples_df_synchronous.tab', 'integrated_phenotype_resamples_df_phyC.tab']: fun.delete_file_or_folder("%s/%s"%(opt.outdir, f))

# save final file
fun.save_df_as_tab(asr_df, integrated_ASR_file)
fun.save_df_as_tab(gwas_df, integrated_GWAS_stats_file)

#################################

# print time
elapsed_time = float(time.time() - start_time)
fun.print_with_runtime("get_GWAS_jobs.py finished correctly in %.4f seconds in outdir %s"%(elapsed_time, opt.outdir))
