#!/usr/bin/env python

##### DEFINE ENVIRONMENT #######

# general module imports
import argparse, os, time
from argparse import RawTextHelpFormatter
import sys
from ete3 import Tree
import pandas as pd
import multiprocessing as multiproc

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import ancestral_GWAS_functions as fun

################################


###### ARGS ######

description = """
Gets a list of jobs (output in a job array file), where each job performs the ancestral state reconstruction (ASR) of one mutation, given a tree and a table with mutations.
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

# mandatory args
parser.add_argument("--variants_table", dest="variants_table", required=True, help="The table with all the variants. It should end with .py or .tab and have a header. The ending will tell the program on how to load the table into a pandas dataframe. ")
parser.add_argument("--phenotypes", dest="phenotypes", required=True, help="The table with the phenotypes. It should be a .tab file with a header and 'sampleID' and 'phenotype' (binary 1/0). It can be none to get all the samples")
parser.add_argument("--tree", dest="tree", required=True, help="The path to a tree in .nw file. If unrooted, this pipeline will root it with the midpoint method.")
parser.add_argument("--outdir", dest="outdir", required=True, help="The output directory, where all files will be written. This script creates a 'correct_tree.nw' file, which includes the tree pruned to have all the samples with mutations and multifuracting for poorly resolved nodes. It also generates a folder called 'ASR_mutations', where the results of ASR for each mutation will be saved. Note that <outdir>/ASR_mutations is the input of the script get_associations.py")
parser.add_argument("--pastml_prediction_method", dest="pastml_prediction_method", required=True, help="The prediction method passed to pastml.")
parser.add_argument("--resampled_phenotypes_pastml_out", dest="resampled_phenotypes_pastml_out", required=True, help="A file, to be passed to get_resampled_phenotypesASR.py")
parser.add_argument("--phenotypes_pastml_prediction_method", dest="phenotypes_pastml_prediction_method", required=True, help="A method, to be passed to get_resampled_phenotypesASR.py")

# optional args
parser.add_argument("--replace", dest="replace", action="store_true", default=False, help="Re-run all the steps by deleting the output directory.")
parser.add_argument("--varID_field", dest="varID_field", type=str, default="variantID_across_samples", help="The field in variants_table that defines the unique mutations' IDs. By default it is 'variantID_across_samples'.")
parser.add_argument("--threads", dest="threads", type=int, default=4, help="The number of threads.")
parser.add_argument("--max_n_mutations", dest="max_n_mutations", type=int, default=None, help="Number of maximum muts to keep")
parser.add_argument("--add_pseudocount_dist_internal_nodes", dest="add_pseudocount_dist_internal_nodes", action="store_true", default=False, help="Add pseudocount dist to all internal nodes.")
parser.add_argument("--only_write_pastml_cmds", dest="only_write_pastml_cmds", action="store_true", default=False, help="Write the pastml cmds and exit")


# parse all the args
opt = parser.parse_args()

##################

# remove outdir if replace, and set replace to False
if opt.replace is True: fun.delete_folder(opt.outdir)

# define the jobs file and remove it
jobs_file = "%s/pastml_cmds.txt"%opt.outdir
fun.remove_file(jobs_file)

# define a final file
final_file = "%s/all_jobs_finished.txt"%opt.outdir
if not fun.file_is_empty(final_file):
    print("The final file of this script is already existing, skip.")
    sys.exit(0)

# define full paths
opt.outdir = fun.get_fullpath(opt.outdir)

# make the outdir
fun.make_folder(opt.outdir)

# define a tmpdir that will be removed at the end
tmpdir = "%s/tmp"%opt.outdir
fun.make_folder(tmpdir)

# load the df with the phenotypes
if opt.phenotypes!="none":

    phenotypes_df = fun.get_tab_as_df_or_empty_df(opt.phenotypes)[["sampleID", "phenotype"]]
    phenotypes_df["sampleID"] = phenotypes_df.sampleID.apply(str)
    strange_phenotypes = set(phenotypes_df.phenotype).difference({0, 1})
    if len(strange_phenotypes)>0: raise ValueError("There are strange phenotypes: %s"%strange_phenotypes)

# check that the pastml prediction methods are correct
allowed_predictionMethods = {'MPPA', 'MAP', 'JOINT', 'DOWNPASS', 'ACCTRAN', 'DELTRAN', 'ALL', 'ML', 'MP'}
if opt.pastml_prediction_method not in allowed_predictionMethods: raise ValueError("The --pastml_prediction_method '%s' is not correct. It should be among %s"%(opt.pastml_prediction_method, allowed_predictionMethods))

# redefine the threads
available_threads = multiproc.cpu_count()
opt.threads = min([opt.threads, available_threads])

######### WRITE CORRECT TREE #########

# define the treefile
correct_treefile = "%s/correct_tree.nw"%opt.outdir
if fun.file_is_empty(correct_treefile):
    print("generating corrected tree")

    # check that the branch supports are correct
    for n in Tree(opt.tree).traverse():
        if n.support<0 or n.support>100: raise ValueError("The provided tree should habe support between 1 and 100. Suport for node %s is %i. %s. %s"%(n, n.support, n.is_root(), n.dist)) 
        if n.support==0: print("WARNING: Node has 0 support:\n%s"%n)

    # load the tree so that it is rooted and has branches with high support
    nchildren_tree = len(Tree(opt.tree).get_children())
    if nchildren_tree==2:  tree = fun.get_correct_tree(opt.tree, min_support=0) # rooted
    elif nchildren_tree>2: tree = fun.get_correct_tree_midpointRooted(opt.tree, min_support=0) # rooted
    else: raise ValueError("The tree root has %i children, which can't be handled"%nchildren_tree)
    fun.check_that_tree_has_no_politomies(tree) # check that there are no politomies

    # keep only samples that are in the phenotypes table
    if opt.phenotypes!="none":interesting_samples = set(phenotypes_df.sampleID)
    else: interesting_samples = set(tree.get_leaf_names())
    tree.prune(interesting_samples, preserve_branch_length=True)
    fun.check_that_tree_has_no_politomies(tree) # check that there are no politomies

    # print  the number of multifuracting nodes
    n_nodes = len([n for n in tree.traverse() if n.is_leaf() is False])
    n_multifurcating_nodes = len([n for n in tree.traverse() if n.is_leaf() is False and len(n.get_children())>2])
    if n_multifurcating_nodes>0: raise ValueError("There are %i/%i multifurcating internal nodes"%(n_multifurcating_nodes, n_nodes))

    # add a pseudocount of the branch length. This is necessary to account for unexisting branch lengths messing with the permutation testing
    pseudocount_dist = min([l.dist for l in tree.get_leaves() if l.dist>0])*0.1
    for l in tree.get_leaves(): l.dist = l.dist + pseudocount_dist

    # add a pseudocount to the internal nodes if required
    if opt.add_pseudocount_dist_internal_nodes is True:
        for n in tree.traverse(): 
            if n.is_leaf() is False and n.is_root() is False: n.dist = n.dist + pseudocount_dist

    # check the lengths of the tree
    for n in tree.traverse():
        if n.dist<=0 and n.is_root() is False: raise ValueError("There are <=0 branch lengths lengths: %s. %s"%(n, n.dist))

    # write tree
    tree.write(outfile=correct_treefile, format=2)

else: interesting_samples = set(Tree(correct_treefile).get_leaf_names())

######################################

######### GET SQUARE VARS DF ###########

# define a file for a df that will have one col for 'sampleID' and one col for each mutation as defined in opt.varID_field
square_vars_df_file = "%s/square_variants.py"%opt.outdir

if fun.file_is_empty(square_vars_df_file): 

    # get the vars df only for the interesting samples and also discarding variants that are in all samples.
    variants_df = fun.get_vars_df_only_subset_samples_for_GWAS(opt.variants_table, interesting_samples, tmpdir, opt.varID_field)

    # add numeric df
    sorted_vars = sorted(set(variants_df[opt.varID_field]))
    numeric_vars = list(map(lambda I: "m%i"%I, range(len(sorted_vars))))
    var_to_numericVar = dict(zip(sorted_vars, numeric_vars))
    variants_df["numeric_mutation"] = variants_df[opt.varID_field].map(var_to_numericVar)
    if any(pd.isna(variants_df["numeric_mutation"])): raise ValueError("there can't be NaNs")

    # write the mapping into a tab file
    numeric_variants_IDmapping_file = "%s/variants_IDmapping.tab"%opt.outdir
    if fun.file_is_empty(numeric_variants_IDmapping_file): fun.save_df_as_tab(variants_df[["numeric_mutation", opt.varID_field]].drop_duplicates().sort_values(by="numeric_mutation"), numeric_variants_IDmapping_file)

    # keep only these fields
    variants_df = variants_df[["sampleID", "numeric_mutation"]]

    # get the square df
    print("getting square variants df")

    # get the df with the variants filtered
    variants_df["mutation_presence"] = 1

    # get the square df
    square_vars_df = variants_df.pivot(index="sampleID", columns="numeric_mutation", values="mutation_presence").fillna(0).applymap(int)
    if fun.get_uniqueVals_df(square_vars_df)!={0, 1}: raise ValueError("The genotypes df should only be 1 and 0")

    # add missing samples as all 0s
    missing_samples = sorted(interesting_samples.difference(set(square_vars_df.index)))
    if len(missing_samples)>0:
        print("WARNING: there are %i samples with no mutations"%(len(missing_samples)))
        df_all0s = pd.DataFrame(np.zeros([len(missing_samples), len(square_vars_df.columns)]), index=missing_samples, columns=square_vars_df.columns).applymap(int)
        square_vars_df = square_vars_df.append(df_all0s)

    # get the correct order and add sampleID field
    square_vars_df = square_vars_df.loc[sorted(interesting_samples), sorted(set(variants_df.numeric_mutation))]
    square_vars_df["sampleID"] = square_vars_df.index
        
    # save
    print("saving")
    fun.save_object(square_vars_df, square_vars_df_file)

# load
square_vars_df = fun.load_object(square_vars_df_file)

########################################

##### DEFINE THE PASTML COMMANDS TO MAKE THE ASR OF EACH MUTATION ######

# define the jobs filename
print("getting pastml cmds")

# define all sorted mutations (sorted by numeric number), and filter if specified
sorted_mutations = sorted([c for c in square_vars_df.columns if c!="sampleID"], key=(lambda x: int(x[1:])))
if opt.max_n_mutations is not None: sorted_mutations = sorted_mutations[0:opt.max_n_mutations]

# get only a set of representative mutations, removing those that are exactly redundant. This also writes mutation_to_representative_mutation.tab in outdir
square_vars_df = fun.get_nr_square_vars_df(square_vars_df[sorted_mutations], opt.outdir, opt.threads)

# redefine the sorted_mutations, to only keep the NR ones.
sorted_mutations = sorted([c for c in square_vars_df.columns], key=(lambda x: int(x[1:])))

# test the timing with N muts
print("testing timing...")
ntest_muts = 1000
outdir_ASR_test = "%s/ASR_mutations_test"%tmpdir; fun.delete_folder(outdir_ASR_test)
elapsed_time_test = fun.get_all_cmds_ASR_mutations(square_vars_df, sorted_mutations[0:ntest_muts], opt.threads, outdir_ASR_test, correct_treefile, opt.pastml_prediction_method)[1]
predicted_time_s = (len(sorted_mutations)/ntest_muts)*elapsed_time_test
print("There are %i mutations. Preparing the pastml cmds on %i threads. This should take ~%.3f min (starting at %s)"%(len(sorted_mutations), opt.threads, predicted_time_s/60, fun.get_date_and_time_for_print()))

# get the cmds
outdir_ASR = "%s/ASR_mutations"%opt.outdir
print("generating mutations ASR's cmds into %s"%outdir_ASR)
all_cmds, elapsed_time = fun.get_all_cmds_ASR_mutations(square_vars_df, sorted_mutations, opt.threads, outdir_ASR, correct_treefile, opt.pastml_prediction_method)
print("cmd generation took %.3f min"%(elapsed_time/60))

# clean the tmp dir
fun.delete_folder(tmpdir)

# write the jobs into a file
if opt.only_write_pastml_cmds is True and len(all_cmds)>0:
    print("Writing pastml cmds into %s and exiting..."%jobs_file)
    jobs_file_tmp = "%s.tmp"%jobs_file
    open(jobs_file_tmp, "w").write("\n".join(all_cmds))
    os.rename(jobs_file_tmp, jobs_file)
    sys.exit(0)

# define the inputs_fn
if len(all_cmds)>0:

    print("Preparing jobs to run locally")
    njobs = len(all_cmds)
    inputs_fn = list(map(lambda I_cmd: tuple(I_cmd[1].split()[1:]+["%.3f"%((I_cmd[0]/njobs)*100)]), enumerate(all_cmds)))

    print("Running pastml jobs in %i threads for %i mutations"%(opt.threads, len(all_cmds)))
    with multiproc.Pool(opt.threads) as pool:
        pool.starmap(fun.run_pastml_py_as_a_function, inputs_fn, chunksize=1)
        pool.close()
        pool.terminate()

########################################################################

############ MAKE RESAMPLES OF PHENOTYPES FOR P VALUE CALCULATIONS ############

# run pastml on 10,000 resampled phenotypes, this is necessary for the p value calculations
if opt.phenotypes!="none":

    if fun.file_is_empty(opt.resampled_phenotypes_pastml_out): 
        print("running the resampling of phenotypes")

        fun.run_cmd("%s --phenotypes %s --resampled_phenotypes_pastml_out %s --tmpdir %s.generating_tmp --pastml_prediction_method %s --threads %i --treefile %s"%(fun.GET_RESAMPLED_PHENOTYPES_ASR_PY, opt.phenotypes, opt.resampled_phenotypes_pastml_out, opt.resampled_phenotypes_pastml_out, opt.phenotypes_pastml_prediction_method, opt.threads, correct_treefile))

###############################################################################

# generate final file
print("All cmds for ASR worked well. You can now execute get_GWAS.py.")
open(final_file, "w").write("all_jobs_finished\n")

# print the message
print("get_ASR_mutations.py finished correctly.")
