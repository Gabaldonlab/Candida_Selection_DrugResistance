#!/usr/bin/env python

# This is the script with functions

########## ENV ##########

# get environment
import sys, os, re, shutil, pickle, warnings, random
import pandas as pd
import copy as cp
import multiprocessing as multiproc
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from ete3 import Tree

# ignore SettingWithCopyWarning
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1])

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])
CondaDir =  "/".join(sys.executable.split("/")[0:-4])
EnvName = EnvDir.split("/")[-1]



#########################

#### FUNCTIONS ####

def remove_file(f):

    if os.path.isfile(f): 

        try: run_cmd("rm %s > /dev/null 2>&1"%f)
        except: pass

def delete_folder(f):

    if os.path.isdir(f): shutil.rmtree(f)


def make_folder(f):

    if not os.path.isdir(f): os.mkdir(f)

def file_is_empty(path): 
    
    """ask if a file is empty or does not exist """

    if not os.path.isfile(path):
        return_val = True
    elif os.stat(path).st_size==0:
        return_val = True
    else:
        return_val = False
            
    return return_val


def run_cmd(cmd, env=EnvName):

    """This function runs a cmd with a given env"""

    # define the cmds
    SOURCE_CONDA_CMD = "source %s/etc/profile.d/conda.sh"%CondaDir
    cmd_prefix = "%s && conda activate %s &&"%(SOURCE_CONDA_CMD, env)

    # define the running
    cmd_to_run = "%s %s"%(cmd_prefix, cmd)

    # run
    out_stat = os.system(cmd_to_run) 
    if out_stat!=0: raise ValueError("\n%s\n did not finish correctly. Out status: %i"%(cmd_to_run, out_stat))

def get_tab_as_df_or_empty_df(file):

    """Gets df from file or empty df"""

    nlines = len([l for l in open(file, "r").readlines() if len(l)>1])

    if nlines==0: return pd.DataFrame()
    else: return pd.read_csv(file, sep="\t")

def get_fullpath(x):

    """Takes a path and substitutes it but the full path"""

    if x.startswith("/"): return x
    elif x.startswith("."): path = "/".join(x.split("/")[1:])
    else: path = x

    return os.getcwd() + "/" + path


def save_object(obj, filename):
    
    """ This is for saving python objects """

    filename_tmp = "%s.tmp"%filename
    remove_file(filename_tmp)
    
    with open(filename_tmp, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

    os.rename(filename_tmp, filename)

def load_object(filename):
    
    """ This is for loading python objects  in a fast way"""
    
    return pickle.load(open(filename,"rb"))




def get_df_and_header_from_vcf(vcf_file):

    """Takes a vcf file and returns the df and the header lines (as a list)"""

    vcf_strings_as_NaNs = ['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'n/a', 'nan', 'null']

    # define the header
    header_lines = [line.strip() for line in open(vcf_file, "r", encoding='utf-8', errors='ignore') if line.startswith("##")]

    # get the vcf
    df = pd.read_csv(vcf_file, skiprows=len(header_lines), sep="\t", na_values=vcf_strings_as_NaNs, keep_default_na=False)
    all_fields = list(df.keys())
    if len(all_fields)!=10: raise ValueError("There should be 10 fields in the vcf %s. Make sure that this is a vcf for a single sample."%file)

    # add fields of the SAMPLE field
    format_fields = set(df.FORMAT.apply(lambda f: tuple(f.split(":"))))
    if len(format_fields)!=1: raise ValueError(">1 FORMAT")
    format_fields = list(next(iter(format_fields)))
    for I,f in enumerate(format_fields): df["SAMPLE_%s"%f] = df[all_fields[-1]].apply(lambda x: x.split(":")[I]).apply(str)

    return df, header_lines

def get_isSNP(r):

    """define if it is SNP"""

    if "," not in r["ALT"]: return (len(r["REF"])==1 and len(r["ALT"])==1)
    else: return (len(r["REF"])==1 and all([len(x)==1 for x in r["ALT"].split(",") ]))

def get_small_vars_df_file_from_df_paths(df_paths, calling_ploidies, tmpdir):

    """Gets a df with variants and wrotes it to a file"""

    # define file
    small_vars_df_file = "%s/small_vars_df.py"%tmpdir

    if file_is_empty(small_vars_df_file):
        print("Generating variants df...")

        # get vars
        sorted_samples = sorted(set(df_paths.sampleID))
        df_paths = cp.deepcopy(df_paths).set_index("sampleID")

        # generate a df with all variants stacked
        small_vars_filt = pd.DataFrame()
        for p in calling_ploidies:
            for sID in sorted_samples:

                # load df with all variants (unfiltered)
                df = get_df_and_header_from_vcf(df_paths.loc[sID, {1:"haploid_vcf", 2:"diploid_vcf"}[p]])[0].rename(columns={"ID":"#Uploaded_variation"})

                # check no multiallelics
                if any(df.ALT.apply(lambda x: "," in x)): raise ValueError("There can't be multiallelics in ALT.")

                # add if it is a SNP
                df["ISSNP"] = df.apply(get_isSNP, axis=1)

                # keep some fields
                df = df[['#Uploaded_variation', '#CHROM', 'POS', 'REF', 'ALT', 'ISSNP', 'SAMPLE_GT']]

                # add fields
                df["sampleID"] = sID
                df["calling_ploidy"] = p

                # append
                small_vars_filt = small_vars_filt.append(df).reset_index(drop=True)

        # save
        save_object(small_vars_filt, small_vars_df_file)

    return small_vars_df_file

def get_dir(filename): return "/".join(filename.split("/")[0:-1])

def get_file(filename): return filename.split("/")[-1]


def get_mosdepth_coverage_for_positions(sorted_bam, coverage_file, positions_df, chrom_field, replace=False, threads=4):

    """Gets the coverage of the positions_df (1-based) into coverage_file from sorted_bam with mosdepth"""

    if file_is_empty(coverage_file) or replace is True:
        #print("getting coverage into %s"%coverage_file)

        # define the prefix
        fileprefix = "%s.MDout"%coverage_file

        # get the positions bed
        windows_file = "%s.windows.tab"%fileprefix
        positions_df["end_POS"] = positions_df.POS+1
        positions_df[[chrom_field, "POS", "end_POS"]].to_csv(windows_file, sep="\t", header=False, index=False)

        # run mosdepth
        mosdepth_std = "%s.std"%fileprefix
        run_cmd("mosdepth --threads %i --by %s --no-per-base --fast-mode  %s %s > %s 2>&1"%(threads, windows_file, fileprefix, sorted_bam, mosdepth_std), env=EnvName) # mosdepth does not look at internal cigar operations or correct mate overlaps (recommended for most use-cases). It is also faster

        # remove all files that don't contribute to coverage
        for sufix in ["mosdepth.global.dist.txt", "mosdepth.region.dist.txt", "mosdepth.summary.txt", "regions.bed.gz.csi", "std", "windows.tab"]: remove_file("%s.%s"%(fileprefix, sufix))

        # unzip
        run_cmd("gunzip %s.regions.bed.gz"%fileprefix, env=EnvName)

        # keep the unzipped as the coverage file
        os.rename("%s.regions.bed"%fileprefix, coverage_file)

    # load
    df_coverage = pd.read_csv(coverage_file, sep="\t", names=[chrom_field, "POS", "POS_plus1", "coverage"])

    return df_coverage

def save_df_as_tab(df, file):

    """Takes a df and saves it as tab"""

    file_tmp = "%s.tmp"%file
    df.to_csv(file_tmp, sep="\t", index=False, header=True)
    os.rename(file_tmp, file)

def get_coverage_per_positions_from_vars_df(vars_df, coverage_df_file, outdir, df_inputs, threads=4):

    """Generates a df written in coverage_df_file with the coverage of all positions in vars_df. df_inputs should have sorted bam"""

    if file_is_empty(coverage_df_file):
        print("getting coverage per position")

        # define a unique df
        positions_df = vars_df[["#CHROM", "POS"]].drop_duplicates().sort_values(by=["#CHROM", "POS"])

        # run mosdepth to get the coverage for each sample
        all_samples = sorted(set(df_inputs.sampleID))

        # define the mosdepth_outdir
        coverage_files_dir = "%s/coverage_files"%outdir; make_folder(coverage_files_dir)

        # rename the input
        df_inputs = df_inputs.set_index("sampleID")

        # run mosdepth in parallel
        inputs_fn = []
        for sampleID in all_samples:

            # define the sorted bam
            sorted_bam = df_inputs.loc[sampleID, "sorted_bam"]
            if file_is_empty(sorted_bam): raise ValueError("%s is not found"%sorted_bam)

            # define the coverage file
            coverage_file = "%s/%s.coverage.%s.tab"%(coverage_files_dir, get_file(coverage_df_file), sampleID)

            # get the inputs
            inputs_fn.append((sorted_bam, coverage_file, cp.deepcopy(positions_df), "#CHROM", False, 1))

        # run the function in parallel
        print("Calculating coverage per position in parallel...")
        with multiproc.Pool(threads) as pool:

            list_coverage_dfs = pool.starmap(get_mosdepth_coverage_for_positions, inputs_fn)
            pool.close()
            pool.terminate()

        # add to positions_df
        for Is, coverage_df_raw in enumerate(list_coverage_dfs):

            # define vars
            sampleID = all_samples[Is]
            coverage_df = coverage_df_raw.rename(columns={"coverage":sampleID})[["#CHROM", "POS", sampleID]]

            # add to positions_df
            positions_df = positions_df.merge(coverage_df, on=["#CHROM", "POS"], validate="one_to_one")

        # clean
        delete_folder(coverage_files_dir)

        # save
        save_df_as_tab(positions_df, coverage_df_file)


    # load
    positions_df = get_tab_as_df_or_empty_df(coverage_df_file)

    return positions_df

def get_small_vars_with_atLeastSomeCoverageInAllPositions(small_vars, outdir, df_sorted_bams, min_coverage_pos, threads=4):

    """Takes a df with small vars (all SNPs) and returns the same only with those that are covered in all samples of df_sorted_bams"""

    # get a df that has the covergae for each position
    coverage_df_file = "%s/coverage_positions.tab"%outdir
    positions_df = get_coverage_per_positions_from_vars_df(small_vars, coverage_df_file, outdir, df_sorted_bams, threads=threads)

    # define all samples
    all_samples = [str(x) for x in sorted(set(df_sorted_bams.sampleID))]

    # keep the positions that have at least min_coverage_pos reads covering in all samples
    positions_df = positions_df[(positions_df[all_samples]>=min_coverage_pos).apply(all, axis=1)]
    print("There are %i positions with coverage >%i in all samples"%(len(positions_df), min_coverage_pos))

    # keep the vars df that are in that positions
    positions_df = positions_df[["#CHROM", "POS"]]
    positions_df["chrom_and_pos"] = positions_df["#CHROM"] + "_" + positions_df.POS.apply(str)
    small_vars["chrom_and_pos"] =  small_vars["#CHROM"] + "_" + small_vars.POS.apply(str)
    small_vars = small_vars[small_vars.chrom_and_pos.isin(set(positions_df.chrom_and_pos))]

    return small_vars

def get_unique_positions_df_from_small_vars_df(small_vars):

    """Takes a small_vars df and returns a df with the positions (in bed style, 0-based and including all the affected position) of the affected vars"""

    # define bed files
    bed_fields = ["#CHROM", "start", "end"]

    # return empty df
    if len(small_vars)==0: return pd.DataFrame(columns=bed_fields)

    # keep non-redundant vars
    small_vars = small_vars[["#CHROM", "POS", "REF", "ALT"]].drop_duplicates()

    # print the strange characters
    acgt = {"A", "C", "G", "T"}
    strange_chars_REF = set.union(*small_vars.REF.apply(set)).difference(acgt)
    strange_chars_ALT = set.union(*small_vars.ALT.apply(set)).difference(acgt)

    if len(strange_chars_REF)>0: print("These are the strange REF characters:%s"%strange_chars_REF)
    if len(strange_chars_ALT)>0: print("These are the strange ALT characters:%s"%strange_chars_ALT)

    # the length of the var define the length of the vars
    small_vars["len_var"] = small_vars.REF.apply(len)

    # define bed-like vars
    small_vars["start"] = small_vars.POS-1
    small_vars["end"] = small_vars.start + small_vars.len_var

    # keep positions
    pos_df = small_vars[bed_fields].sort_values(by=bed_fields).drop_duplicates()

    return pos_df

def get_snps_df_outside_targetPositions(snps_df, wrong_positions_df, tmpdir):

    """This function keeps only the small vars that are not in wrong_positions_df """

    # if there are no positions, do not change
    if len(wrong_positions_df)==0: return snps_df

    # make files
    make_folder(tmpdir)

    positions_bed = "%s/wrong_positions.bed"%tmpdir
    wrong_positions_df[["#CHROM", "start", "end"]].sort_values(by=["#CHROM", "start", "end"]).to_csv(positions_bed, sep="\t", header=False, index=False)

    unique_vars_bed = "%s/unique_vars.bed"%tmpdir
    snps_df["start"] = snps_df.POS-1
    bed_fields = ["#CHROM", "start", "POS", "#Uploaded_variation"]
    snps_df[bed_fields].sort_values(by=bed_fields).drop_duplicates().to_csv(unique_vars_bed, sep="\t", header=False, index=False)

    # run bedmap to find variants in  positions_bed overlapped by unique_vars_bed
    bedmap_outfile = "%s/vars_overlapping_wrong_positions.bed"%tmpdir
    run_cmd("bedmap --range 0 --echo-map-id --delim '\t' %s %s > %s"%(positions_bed, unique_vars_bed, bedmap_outfile), env=EnvName)

    wrong_vars = set.union(*pd.read_csv(bedmap_outfile, sep="\t", header=None, names=["wrong_var"]).wrong_var.apply(lambda x: x.split(";")).apply(set))

    snps_df = snps_df[~snps_df["#Uploaded_variation"].isin(wrong_vars)]

    # clean
    delete_folder(tmpdir)

    return snps_df


def generate_haploid_snps_df_uniFormPositions_file(small_vars_file, df_sorted_bams, haploid_snps_df_uniFormPositions_file, min_coverage_pos, outdir, tree_mode, threads=4):

    """This function generates a df with snps from the filtered small vars, only keeping positions that have at least >min_coverage_pos coverage, have no indels nor heterozygous SNPs. """

    if file_is_empty(haploid_snps_df_uniFormPositions_file):
        make_folder(outdir)
        print("generating %s on %i threads"%(haploid_snps_df_uniFormPositions_file, threads))

        # define all the samples
        all_samples = set(df_sorted_bams.sampleID)

        # load the variants df (already filtered)
        print("loading variants")
        small_vars = load_object(small_vars_file)

        # add some fields
        small_vars["SAMPLE_GT"] = small_vars.SAMPLE_GT.apply(str)

        def get_upper(x): return x.upper()
        small_vars["REF"] = small_vars.REF.apply(get_upper)
        small_vars["ALT"] = small_vars.ALT.apply(get_upper)
        small_vars["sampleID"] = small_vars.sampleID.apply(str)

        # remove variants that are not in all samples
        small_vars = small_vars[small_vars.sampleID.isin(all_samples)]

        # define the target ploidy depending on the tree mode
        tree_ploidy = {"haploid":1, "diploid_homozygous":2}[tree_mode]

        # initialize the snps_df that are homozygous and of the interesting ploiyd, and ACGT
        acgt = {"A", "C", "G", "T"}
        haploid_snps = small_vars[(small_vars.SAMPLE_GT.isin({"1/1", "1"})) & (small_vars.ISSNP) & (small_vars.REF.isin(acgt)) & (small_vars.ALT.isin(acgt)) & (small_vars.calling_ploidy==tree_ploidy)]
        ntotal_snps = len(set(haploid_snps["#Uploaded_variation"]))
        print("There are %i total haploid SNPs"%ntotal_snps)

        # filter out SNPs that are in positions with some INDEL in some sample
        positions_INDELs = get_unique_positions_df_from_small_vars_df(small_vars[(small_vars.calling_ploidy==tree_ploidy) & ~(small_vars.ISSNP)])
        haploid_snps = get_snps_df_outside_targetPositions(haploid_snps, positions_INDELs, "%s/removing_positionsWithINDELS"%outdir)

        print("There are %i/%i SNPs remaining after filtering out positions with INDELs"%(len(set(haploid_snps["#Uploaded_variation"])), ntotal_snps))

        # filter out SNPs that are in positions with some heterozygous snp in some sample
        positions_df_heteroSNPs = get_unique_positions_df_from_small_vars_df(small_vars[~(small_vars.SAMPLE_GT.isin({"1/1", "1"})) & (small_vars.calling_ploidy==2) & (small_vars.ISSNP)])
        haploid_snps = get_snps_df_outside_targetPositions(haploid_snps, positions_df_heteroSNPs, "%s/removing_positionsWithHetSNPs"%outdir)

        print("There are %i/%i SNPs remaining after filtering out positions with het SNPs"%(len(set(haploid_snps["#Uploaded_variation"])), ntotal_snps))

        # filter out SNPs that are in positions not covered in all samples
        haploid_snps = get_small_vars_with_atLeastSomeCoverageInAllPositions(haploid_snps, outdir, df_sorted_bams, min_coverage_pos, threads=threads)

        print("There are %i/%i SNPs remaining after filtering out positions with low coverage in some samples"%(len(set(haploid_snps["#Uploaded_variation"])), ntotal_snps))

        # save
        save_object(haploid_snps, haploid_snps_df_uniFormPositions_file)

        del haploid_snps
        del small_vars

def get_series_ALTsequence(sampleID, positions_df, snps_df, pickRandomHetSNPs):

    """This function takes a positions_df and the samples with snps and returns a series where the index is the chrom_and_pos and the value is an ACTG alternative allele. snps_df should have [["chrom_and_pos", "ALT"]] """

    # report
    #print("working on sample %s"%sampleID)
    
    # redefine snps_df to keep one random heterozygous SNP
    if pickRandomHetSNPs is True:

        # define the number of het SNPs
        n_hetSNPs = sum(snps_df.SAMPLE_GT=="0/1")
        n_nonHet_SNPs = sum(snps_df.SAMPLE_GT!="0/1")

        # get the SNPs so that heterozygous SNPs get randomly kept
        def get_boolean_for_GT(gt):
            if gt=="1/1": return True
            elif gt=="0/1": return bool(random.getrandbits(1))
            else: raise ValueError("GT not defined: %s"%gt) 

        snps_df = snps_df[snps_df.SAMPLE_GT.apply(get_boolean_for_GT)]

        # print the fraction of heterozygous SNPs kept
        print("You kept %i/%i het SNPs and %i/%i other SNPs in sample %s"%(sum(snps_df.SAMPLE_GT=="0/1"), n_hetSNPs, sum(snps_df.SAMPLE_GT!="0/1"), n_nonHet_SNPs, sampleID))

    # debug
    if len(snps_df)==0: raise ValueError("there are no snps in %s"%sampleID)

    # keep only the first SNP of each position if it is multiallelic
    snps_df = snps_df.sort_values(by=["chrom_and_pos", "ALT"]).drop_duplicates(subset=["chrom_and_pos"], keep="first")

    # check that the positions are unique
    if len(snps_df)!=len(set(snps_df.chrom_and_pos)): raise ValueError("there are some repeated snps_df.chrom_and_pos")

    # define the series so that the positions without SNPs have NaNs
    initial_len_positions_df = len(positions_df)
    positions_df = positions_df.merge(snps_df[["chrom_and_pos", "ALT"]], on="chrom_and_pos", how="left").rename(columns={"ALT":sampleID})

    # debug merge
    if len(positions_df)!=initial_len_positions_df: raise ValueError("some missing positions in positions_df")

    # report
    print("there are %i/%i positions with snps in sample %s"%(sum(~pd.isna(positions_df[sampleID])), len(positions_df), sampleID))

    # correct, so that the NaNs get the reference allele
    def get_corrected_nt(r):
        if pd.isna(r[sampleID]): return r["REF"]
        else: return r[sampleID]
    positions_df[sampleID] = positions_df.apply(get_corrected_nt, axis=1)

    # debug
    actg = {"A", "C", "T", "G"}
    if any(~positions_df[sampleID].isin(actg)): raise ValueError("all positions should be ACTG")

    return positions_df.set_index("chrom_and_pos")[sampleID]


def generate_multifasta_from_snps_df(snps_df, multifasta_correct_SNPs, all_samples, threads=4, pickRandomHetSNPs=False, generate_one_aln_each_chrom=False):

    """Takes a snps df and a df with SNPs and generates a multifasta of them."""

    # define a df that will have all the seqs
    positions_df_file = "%s.positions_df.py"%multifasta_correct_SNPs

    if file_is_empty(multifasta_correct_SNPs):
        print("getting multifasta of the correct SNPs on %i threads"%threads)

        # check that all the samples are in all_samples
        strange_samples = set(snps_df.index).difference(set(all_samples))
        if len(strange_samples)>0: raise ValueError("there are some undefined samples in snps_df: %s"%strange_samples)

        # define a positions_df with all the positions with some SNP
        positions_df = snps_df[["#CHROM", "POS"]].drop_duplicates().sort_values(by=["#CHROM", "POS"])
        positions_df["chrom_and_pos"] = positions_df["#CHROM"] + "_" + positions_df.POS.apply(str)

        # map each chrom_pos to the reference
        snps_df["chrom_and_pos"] =  snps_df["#CHROM"] + "_" + snps_df.POS.apply(str)
        chromPos_to_ref = dict(snps_df[["chrom_and_pos", "REF"]].drop_duplicates().set_index("chrom_and_pos")["REF"])
        positions_df["REF"] = positions_df.chrom_and_pos.apply(lambda x: chromPos_to_ref[x])

        # define samples with vars
        all_samples_withVars = set(snps_df.index)

        # get list of ALT sequences in parallel
        snps_df = snps_df[["chrom_and_pos", "ALT", "SAMPLE_GT"]]
        inputs_fn = [(sampleID, positions_df, snps_df.loc[{sampleID}], pickRandomHetSNPs) for sampleID in all_samples if sampleID in all_samples_withVars]

        print("getting SNPs series")
        if threads>1:

            with multiproc.Pool(threads) as pool:

                list_ALTseries = pool.starmap(get_series_ALTsequence, inputs_fn)
                pool.close()
                pool.terminate()

        else: list_ALTseries = list(map(lambda x: get_series_ALTsequence(x[0], x[1], x[2], x[3]), inputs_fn))

        # add the alts to positions_df
        print("adding series to positions df")
        positions_df = positions_df.sort_values(by=["#CHROM", "POS"]).set_index("chrom_and_pos", drop=False)
        for series in list_ALTseries: positions_df[series.name] = series

        # add the samples with no vars
        for s in set(all_samples).difference(all_samples_withVars): positions_df[s] = positions_df.REF

        # write the positions of the fasta
        save_object(positions_df, positions_df_file)

        """

        # save the alignment for each chromosome
        if generate_one_aln_each_chrom is True:

            for chrom in set(positions_df["#CHROM"]):
                print("saving multifasta for %s"%chrom)

                positions_df_chrom = positions_df[positions_df["#CHROM"]==chrom]
                all_records = [SeqRecord(Seq("".join(positions_df_chrom[sampleID].values)), id=str(sampleID), name="", description="") for sampleID in all_samples]
                
                SeqIO.write(all_records, "%s.%s.fasta"%(multifasta_correct_SNPs, chrom), "fasta")
        """

        # generate the multifasta and write
        print("saving records of the concatenated alignment")
        all_records = [SeqRecord(Seq("".join(positions_df[sampleID].values)), id=str(sampleID), name="", description="") for sampleID in all_samples]
        multifasta_correct_SNPs_tmp = "%s.tmp"%multifasta_correct_SNPs
        SeqIO.write(all_records, multifasta_correct_SNPs_tmp, "fasta")
        os.rename(multifasta_correct_SNPs_tmp, multifasta_correct_SNPs)

    return positions_df_file



def get_multifasta_onlyVariableSites(positions_df_file, multifasta, sorted_samples, replace=False):

    """Gets a df positions, such as returned by generate_multifasta_from_snps_df, and writes multifasta """

    if file_is_empty(multifasta) or replace is True:

        # load the positions df and sort
        print("loading positions df")
        positions_df = load_object(positions_df_file)

        # sort 
        positions_df = positions_df.sort_values(by=["#CHROM", "POS"])

        # keep only variable positions
        total_npositions = len(positions_df)
        positions_df["n_different_nucleotides"] = positions_df[sorted_samples].apply(pd.unique, axis=1).apply(len)
        positions_df = positions_df[positions_df.n_different_nucleotides>1]
        print("there are %i/%i variable sites in the alignment"%(len(positions_df), total_npositions))

        # generate the multifasta and write
        all_records = [SeqRecord(Seq("".join(positions_df[sampleID].values)), id=str(sampleID), name="", description="") for sampleID in sorted_samples]
        multifasta_tmp = "%s.tmp"%multifasta
        SeqIO.write(all_records, multifasta_tmp, "fasta")
        os.rename(multifasta_tmp, multifasta)

        del positions_df

def get_correct_tree_midpointRooted(treefile, min_support=0):

    """Gets a tree with midpointRooting"""


    # load
    tree = Tree(treefile)

    # get the midpoint rooted tree
    midpoint =  tree.get_midpoint_outgroup()
    tree.set_outgroup(midpoint)

    # collapse low support bracnhes
    for n in tree.traverse():

        # for internal nodes, if the support is below the minimum, set each of the children as children of the parent
        if n.is_leaf() is False and len(n.get_ancestors())>0 and n.support<min_support: n.delete()

    return tree


def generate_homo_and_hetero_snps_df_correctPositions(small_vars_file, df_sorted_bams, homo_and_hetero_snps_df_correctPositions_file, min_coverage_pos, outdir, threads=4):

    """This function generates a df with diploid snps from the filtered small vars, only keeping positions that have at least >min_coverage_pos coverage, have no indels. """

    if file_is_empty(homo_and_hetero_snps_df_correctPositions_file):
        make_folder(outdir)
        print("generating %s on %i threads"%(homo_and_hetero_snps_df_correctPositions_file, threads))

        # define all the samples
        all_samples = set(df_sorted_bams.sampleID)

        # load the variants df (already filtered)
        print("loading variants")
        small_vars = load_object(small_vars_file)

        # add some fields
        small_vars["SAMPLE_GT"] = small_vars.SAMPLE_GT.apply(str)

        def get_upper(x): return x.upper()
        small_vars["REF"] = small_vars.REF.apply(get_upper)
        small_vars["ALT"] = small_vars.ALT.apply(get_upper)
        small_vars["sampleID"] = small_vars.sampleID.apply(str)

        # remove variants that are not in all samples
        small_vars = small_vars[small_vars.sampleID.isin(all_samples)]

        # initialize the snps_df that are diploid 2 and of the interesting ploiyd, and ACGT
        acgt = {"A", "C", "G", "T"}
        diploid_snps = small_vars[(small_vars.SAMPLE_GT.isin({"1/1", "0/1"})) & (small_vars.ISSNP) & (small_vars.REF.isin(acgt)) & (small_vars.ALT.isin(acgt)) & (small_vars.calling_ploidy==2)]

        ntotal_snps = len(set(diploid_snps["#Uploaded_variation"]))
        print("There are %i total diploid SNPs"%ntotal_snps)

        # filter out SNPs that are in positions with some INDEL in some sample
        positions_INDELs = get_unique_positions_df_from_small_vars_df(small_vars[(small_vars.calling_ploidy==2) & ~(small_vars.ISSNP)])
        diploid_snps = get_snps_df_outside_targetPositions(diploid_snps, positions_INDELs, "%s/removing_positionsWithINDELS"%outdir)

        print("There are %i/%i SNPs remaining after filtering out positions with INDELs"%(len(set(diploid_snps["#Uploaded_variation"])), ntotal_snps))

        # filter out SNPs that are in positions not covered in all samples
        diploid_snps = get_small_vars_with_atLeastSomeCoverageInAllPositions(diploid_snps, outdir, df_sorted_bams, min_coverage_pos, threads=threads)

        print("There are %i/%i SNPs remaining after filtering out positions with low coverage in some samples"%(len(set(diploid_snps["#Uploaded_variation"])), ntotal_snps))

        # save
        save_object(diploid_snps, homo_and_hetero_snps_df_correctPositions_file)

        del diploid_snps
        del small_vars


def soft_link_files(origin, target):

    """This function takes an origin file and makes it accessible through a link (target)"""

    if file_is_empty(target):

        # rename as full paths
        origin = get_fullpath(origin)
        target = get_fullpath(target)

        # check that the origin exists
        if file_is_empty(origin): raise ValueError("The origin %s should exist"%origin)

        # remove previous lisqnk
        try: run_cmd("rm %s > /dev/null 2>&1"%target)
        except: pass

        soft_linking_std = "%s.softlinking.std"%(target)
        #print("softlinking. The std is in %s"%soft_linking_std)
        run_cmd("ln -s %s %s > %s 2>&1"%(origin, target, soft_linking_std))
        remove_file(soft_linking_std)

    # check that it worked
    if file_is_empty(target): raise ValueError("The target %s should exist"%target)

def get_chr_to_len(genome, replace=False):

    chr_to_len_file = "%s.chr_to_len.py"%genome
    chr_to_len_file_tmp = "%s.tmp"%chr_to_len_file

    if file_is_empty(chr_to_len_file) or replace is True:

        remove_file(chr_to_len_file_tmp)

        # define chromosome_to_length for a genome
        chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(genome, "fasta")}

        # save
        save_object(chr_to_len, chr_to_len_file_tmp)
        os.rename(chr_to_len_file_tmp, chr_to_len_file)

    else: chr_to_len = load_object(chr_to_len_file)

    return chr_to_len


def index_genome(genome):

    """Takes a fasta and generates a <genome>.fai file"""

    # index the genome of interest if not already done
    if file_is_empty("%s.fai"%genome): 

        faidx_std = "%s.indexing.std"%genome
        #print("running faidx. The std is in %s"%faidx_std)
        run_cmd("samtools faidx %s > %s 2>&1"%(genome, faidx_std), env=EnvName) # samtools=1.9 in the orignial one. But this one installs 1.6

        remove_file(faidx_std)

def create_sequence_dict(genome, replace=False):

    """Takes a fasta and generates the reference dict"""

    rstrip = genome.split(".")[-1]
    dictionary = "%sdict"%(genome.rstrip(rstrip)); tmp_dictionary = "%s.tmp"%dictionary

    if file_is_empty(dictionary) or replace is True:

        # remove any previously created tmp_file
        remove_file(tmp_dictionary)

        # define the std
        dictionary_std = "%s.generating.std"%dictionary
        print("Creating picard dictionary. The std is in %s"%dictionary_std)

        run_cmd("picard CreateSequenceDictionary R=%s O=%s TRUNCATE_NAMES_AT_WHITESPACE=true > %s 2>&1"%(genome, tmp_dictionary, dictionary_std), env=EnvName) # picard version

        remove_file(dictionary_std)  
        os.rename(tmp_dictionary , dictionary)

def generate_fasta_with_subset_positions(input_fasta, output_fasta, positions_bed):

    """Takes a genome and keeps only the positions in positions_bed"""

    if file_is_empty(output_fasta):
        print("generating subset positions")

        # get the fasta of each position
        fasta_all_positions = "%s.all_positions.fasta"%output_fasta
        run_cmd("bedtools getfasta -fi %s -bed %s > %s"%(input_fasta, positions_bed, fasta_all_positions), env=EnvName)

        # rewrite so that there is only the chromosome name as ID
        def get_seq_with_CorrectID(seq): 
            seq.id = seq.id.split(":")[0]
            seq.name = ""
            seq.description = ""
            return seq

        all_records = list(map(get_seq_with_CorrectID, SeqIO.parse(fasta_all_positions, "fasta")))
        SeqIO.write(all_records, fasta_all_positions, "fasta")

        # collapse records of the same type
        output_fasta_tmp = "%s.tmp"%output_fasta
        cmd_get_fasta_with_positions = 'cat %s | paste - - | datamash -g 1 collapse 2 | tr -d "," | tr "\t" "\n" > %s'%(fasta_all_positions, output_fasta_tmp)
        run_cmd(cmd_get_fasta_with_positions)

        # check that all the chromosomes are there as expected
        chroms_input = set(get_chr_to_len(input_fasta))
        chroms_output_tmp = set(get_chr_to_len(output_fasta_tmp))

        if chroms_input!=chroms_output_tmp: raise ValueError("ERROR in %s\n. chroms_input is %s and chroms_output_tmp is %s. They should be the same"%(output_fasta_tmp, chroms_input, chroms_output_tmp))

        # check that the file is not empty
        if file_is_empty(output_fasta_tmp): raise ValueError("output_fasta_tmp can't be empty")

        remove_file(fasta_all_positions)
        os.rename(output_fasta_tmp, output_fasta)

def rsync_file(origin, dest):

    """syncs one file to the other"""

    dest_tmp = "%s.tmp"%dest
    run_cmd("rsync %s %s"%(origin, dest_tmp))
    os.rename(dest_tmp, dest)


def get_bgzip_and_and_tabix_vcf_file(file, replace=False):

    """Takes a vcf file and returns a tabixed and gzipped file"""

    # define files
    file_gz = "%s.gz"%file
    file_tmp_gz = "%s.tmp.gz"%file
    file_gz_tbi = "%s.gz.tbi"%file
    file_tmp_gz_tbi = "%s.tmp.gz.tbi"%file

    if file_is_empty(file_gz) or file_is_empty(file_gz_tbi) or replace is True:

        # bgzip
        bgzip_stderr = "%s.generating.stderr"%file_tmp_gz
        print("bgzipping. The stderr is in %s"%bgzip_stderr)
        run_cmd("bgzip -c %s > %s 2>%s"%(file, file_tmp_gz, bgzip_stderr), env=EnvName)

        # tabix
        tabix_std = "%s.tabixing.std"%file_tmp_gz
        print("tabix-ing. The std is in %s"%tabix_std)
        run_cmd("tabix -p vcf %s > %s 2>&1"%(file_tmp_gz, tabix_std), env=EnvName)

        # remove files
        remove_file(bgzip_stderr)
        remove_file(tabix_std)

        # rename
        os.rename(file_tmp_gz, file_gz)
        os.rename(file_tmp_gz_tbi, file_gz_tbi)

    return file_gz, file_gz_tbi


def get_alternative_genome_FastaAlternateReferenceMaker(sampleID, snps_df, reference_genome, tmpdir, pickRandomHetSNPs, positions_bed):

    """Takes a df with SNPs and generates a multifasta with the genome changed, only positions in intervals_file"""

    # define files
    vcf_file = "%s/sample_%s_variants.vcf"%(tmpdir, sampleID)
    alternative_fasta = "%s.transformed_seq.fasta"%vcf_file

    # get the alternative genome
    if file_is_empty(alternative_fasta):

        if len(snps_df)==0: rsync_file(reference_genome, alternative_fasta)

        else:

            # remove randomly half of the hetero SNPs if indicated
            if pickRandomHetSNPs is True:

                # define the number of het SNPs
                n_hetSNPs = sum(snps_df.SAMPLE_GT=="0/1")
                n_nonHet_SNPs = sum(snps_df.SAMPLE_GT!="0/1")

                # get the SNPs so that heterozygous SNPs get randomly kept
                def get_boolean_for_GT(gt):
                    if gt in {"1/1", "1"}: return True
                    elif gt in {"0/1"}: return bool(random.getrandbits(1))
                    else: raise ValueError("GT not defined: %s"%gt) 

                snps_df = snps_df[snps_df.SAMPLE_GT.apply(get_boolean_for_GT)]

                # print the fraction of heterozygous SNPs kept
                print("You kept %i/%i het SNPs and %i/%i other SNPs in sample %s"%(sum(snps_df.SAMPLE_GT=="0/1"), n_hetSNPs, sum(snps_df.SAMPLE_GT!="0/1"), n_nonHet_SNPs, sampleID))

            # add vcf fields
            vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "DATA"]
            snps_df["ID"] = snps_df["#Uploaded_variation"]
            snps_df["QUAL"] = 1000
            snps_df["FILTER"] = "PASS"
            snps_df["INFO"] = "."
            snps_df["FORMAT"] = "GT"
            snps_df["DATA"] = 1

            def get_upper(x): return x.upper()
            snps_df["REF"] = snps_df.REF.apply(get_upper)
            snps_df["ALT"] = snps_df.ALT.apply(get_upper)

            # write to vcf
            vcf_df = snps_df[vcf_fields].drop_duplicates().sort_values(by=["#CHROM", "POS"])
            vcf_lines = vcf_df.to_csv(sep="\t", header=False, index=False)
            header_lines = "\n".join(["##fileformat=VCFv4.2", "\t".join(vcf_fields)])
            open(vcf_file, "w").write(header_lines + "\n" + vcf_lines)

            # index vcf
            vcf_file_gz, vcf_file_gx_tbi = get_bgzip_and_and_tabix_vcf_file(vcf_file, replace=True)

            # get the alternative genome        
            alternative_fasta_tmp = "%s.tmp.fasta"%alternative_fasta
            make_alternative_fasta_std = "%s.generating.std"%alternative_fasta

            print("getting alternative genome. STD in %s"%make_alternative_fasta_std)
            run_cmd("gatk FastaAlternateReferenceMaker -R %s -O %s -V %s > %s 2>&1"%(reference_genome, alternative_fasta_tmp, vcf_file_gz, make_alternative_fasta_std), env=EnvName)

            # check that the fasta is correct and change the names
            print("changing names and testing that everything went well")
            chrom_to_refSeq = {seq.id : seq for seq in SeqIO.parse(reference_genome, "fasta")}
            chrom_to_len = get_chr_to_len(reference_genome)

            all_records = []
            all_new_chroms = set()
            for seq in SeqIO.parse(alternative_fasta_tmp, "fasta"):

                # change the seqID
                seq.id = seq.description.split()[1].split(":")[0]
                seq.name = ""
                seq.description = ""

                # keep
                all_new_chroms.add(seq.id)
                all_records.append(seq)

                # check that the length is the same
                if len(seq)!=chrom_to_len[seq.id]: raise ValueError("there are different lengths in chrom %s"%(seq.id))

                # check that the variant positions match those of the SNPs
                alternative_seq = pd.Series(list(seq.seq), index=range(1, len(seq)+1))
                ref_seq = pd.Series(list(chrom_to_refSeq[seq.id].seq), index=range(1, len(seq)+1))

                real_variable_positions = set(ref_seq[alternative_seq!=ref_seq].index)
                expected_variable_positions = set(snps_df[snps_df["#CHROM"]==(seq.id)].POS)

                strange_real_variable_positions = real_variable_positions.difference(expected_variable_positions)
                if len(strange_real_variable_positions)>0: 

                    # define a df with the strange_real_variable_positions, keeping the ones were the REF is ACGT
                    df_seqs = pd.DataFrame({"REF":ref_seq.loc[strange_real_variable_positions], "ALT":alternative_seq.loc[strange_real_variable_positions]}).sort_values(by=["REF", "ALT"])

                    df_seqs = df_seqs[df_seqs.REF.isin({"A", "C", "T", "G"})]

                    # if there are some of these positons, raise ERROR
                    if len(df_seqs)>0:

                        print("ERROR: There are %i unexpected variable positions in chromosome %s. These are the sequences found in these positions:\n%s"%(len(df_seqs), seq.id, df_seqs))

                        raise ValueError("there are unexpected variable positions")

                missing_expected_variable_positions = expected_variable_positions.difference(real_variable_positions)
                if len(missing_expected_variable_positions)>0: raise ValueError("There are some variable positions that are missing %s: %s"%(seq.id, missing_expected_variable_positions))

            # check that the defined chromosomes are correct
            if all_new_chroms!=set(chrom_to_len): raise ValueError("there are some unexpected chromosomes")

            # write
            SeqIO.write(all_records, alternative_fasta_tmp, "fasta")

            # clean
            remove_file(make_alternative_fasta_std)
            os.rename(alternative_fasta_tmp, alternative_fasta)
        
    # get a subset of the alternative genome in positions_bed
    alternative_fasta_subset = "%s.onlySubsetPositions.fasta"%alternative_fasta
    generate_fasta_with_subset_positions(alternative_fasta, alternative_fasta_subset, positions_bed)

    print("the alternative multifasta generation for sample %s went well"%sampleID)

    return alternative_fasta_subset

def generate_multifasta_from_snps_df_file_FastaAlternateReferenceMaker(snps_df_file, reference_genome, multifasta, sorted_samples, threads=4, expected_GTs={"1/1", "0/1"}, pickRandomHetSNPs=False):

    """Generates a multifasta with the concatenated reference genome (whole sequence) and the snps file. This fasta contains only positions that have some SNP in snps_df_file."""

    if file_is_empty(multifasta):
        print("generatig multifasta")

        ######### FILTER SNPS #########

        # define a tmpdir
        tmpdir = "%s.generating"%multifasta
        #delete_folder(tmpdir)
        make_folder(tmpdir)

        # index the genome
        index_genome(reference_genome)
        create_sequence_dict(reference_genome)

        # load snps
        print("loading SNPs")
        snps_df = load_object(snps_df_file)

        # define a bed with the positions with snps
        snp_positions_df = snps_df[["#CHROM", "POS"]].sort_values(by=["#CHROM", "POS"]).drop_duplicates()
        snp_positions_bed = "%s/SNP_positions.bed"%tmpdir
        snp_positions_df["start"] = snp_positions_df.POS-1
        snp_positions_df[["#CHROM", "start", "POS"]].to_csv(snp_positions_bed, sep="\t", header=False, index=False)

        # add fields
        snps_df["sampleID"] = snps_df.sampleID.apply(str)
        snps_df = snps_df.set_index("sampleID", drop=False)

        # check that all the samples are as expected
        strange_samples = set(snps_df.index).difference(set(sorted_samples))
        if len(strange_samples)>0: raise ValueError("there are some undefined samples in snps_df: %s"%strange_samples)

        # check that the genotypes make sense
        strange_GTs = set(snps_df.SAMPLE_GT).difference(expected_GTs)
        if len(strange_GTs)>0: raise ValueError("There are unexpected GTs: %s"%strange_GTs)

        # check the GTs
        all_GTs = set(snps_df.SAMPLE_GT)
        if len(all_GTs)==1: print("WARNING: There is only one type of GTs: %s. This means that this re-sampling procedure of homozygous/heterozygous SNPs is not valid."%all_GTs)

        ##############

        ########## GET MULIFASTA #######

        # keep important fields
        snps_df = snps_df[["#Uploaded_variation", "#CHROM", "POS", "REF", "ALT", "SAMPLE_GT"]]

        # get list of individual fasta sequences transformed to include the SNPs
        inputs_fn = [(sampleID, snps_df.loc[{sampleID}], reference_genome, tmpdir, pickRandomHetSNPs, snp_positions_bed) for sampleID in sorted_samples]

        if threads>1:

            with multiproc.Pool(threads) as pool:

                list_alternative_fastaFiles = pool.starmap(get_alternative_genome_FastaAlternateReferenceMaker, inputs_fn)
                pool.close()
                pool.terminate()

        else: list_alternative_fastaFiles = list(map(lambda x: get_alternative_genome_FastaAlternateReferenceMaker(x[0], x[1], x[2], x[3], x[4], x[5]), inputs_fn))

        sampleID_to_altGenome = dict(zip(sorted_samples, list_alternative_fastaFiles))

        ################################

        ######### CONCATENATE MULTIFASTA BY CHROMS ##########

        # define the sorted chroms
        sorted_chroms = sorted(set(get_chr_to_len(reference_genome)))

        # define a dict that maps all the records
        sample_to_chrom_to_seq = {sampleID : {seq.id : seq for seq in SeqIO.parse(altGenome, "fasta")}   for sampleID, altGenome in sampleID_to_altGenome.items()}

        # generate one fasta for each chrom with all the samples
        print("generating one sample for each chromosome")
        for chrom in sorted_chroms:
            print(chrom)

            all_records = []
            all_lengths = set()
            for sampleID in sorted_samples:

                seq = sample_to_chrom_to_seq[sampleID][chrom]
                seq.id = str(sampleID)
                seq.name = ""
                seq.description = ""
                all_records.append(seq)
                all_lengths.add(len(seq.seq))

            # check and write
            if len(all_lengths)!=1: raise ValueError("there are strange lengths")
            SeqIO.write(all_records, "%s.%s.fasta"%(multifasta, chrom), "fasta")

        # generate one single multifasta with all concatenated chromosomes
        print("generating one single multifasta with all records")
        all_records = []
        all_lengths = set()

        for sampleID in sorted_samples:
            print(sampleID)

            # init a seq record
            str_seq = "".join([str(sample_to_chrom_to_seq[sampleID][chrom].seq) for chrom in sorted_chroms])
            concatenated_seq = SeqRecord(Seq(str_seq), id=str(sampleID), name="", description="")

            all_records.append(concatenated_seq)
            all_lengths.add(len(concatenated_seq.seq))

        if all_lengths!={len(snp_positions_df)}: raise ValueError("the length of the final alignment has less positions than expected")
        multifasta_tmp = "%s.tmp"%multifasta
        SeqIO.write(all_records, multifasta_tmp, "fasta")


        #####################################################

        # clean
        delete_folder(tmpdir)
        os.rename(multifasta_tmp, multifasta)



def get_multifasta_onlyVariableSites_snpSites(multifasta):

    """This function gets positions withput SNPs from MSA"""

    # define final file
    multifasta_onlyVariableSites = "%s.onlyVariableSites.fasta"%multifasta

    if file_is_empty(multifasta_onlyVariableSites):

        # get the multifasta
        multifasta_onlyVariableSites_tmp = "%s.tmp"%multifasta_onlyVariableSites
        run_cmd("snp-sites -m -o %s %s"%(multifasta_onlyVariableSites_tmp, multifasta), env=EnvName)

        os.rename(multifasta_onlyVariableSites_tmp, multifasta_onlyVariableSites)

    return multifasta_onlyVariableSites




def generate_consensus_withBootstrap_from_resampledTrees(trees_list, outdir, replace=False):

    """Takes the consensus tree from the tree_list it uses all the trees to calculate bootstrap supports and branch lengths"""

    ######## GET ALL THE TREES ROOTED AND ADDED TO A SINGLE FILE #########

    all_trees_file = "%s/all_merged_trees_file.txt"%outdir

    if file_is_empty(all_trees_file) or replace is True:
        print("getting individual trees in a file")

        # init tmp
        all_trees_file_tmp = "%s.tmp"%all_trees_file
        remove_file(all_trees_file_tmp)
        run_cmd("touch %s"%all_trees_file_tmp)

        for Itree, treefile in enumerate(trees_list): 
            print("adding tree %i"%(Itree+1))

            # get the rooted tree
            tree = Tree(treefile)
            tree.set_outgroup(tree.get_midpoint_outgroup())

            # check that all nodes have 2 children
            for n in tree.traverse():
                if not n.is_leaf() and len(n.get_children())!=2: raise ValueError("There are some collapsed branches in the tree")

            # add to the file with all trees
            open(all_trees_file_tmp, "a").write(tree.write(format=2)+"\n")

        os.rename(all_trees_file_tmp, all_trees_file)

    ######################################################################

    ######### GET CONSENSUS TREE WITH SUPPORT ############

    # define files
    consensus_treefile_withBL = "%s/tree_consensus_withBootstraps_and_branchLengths.nw"%outdir
    consensus_treefile = "%s.contree"%all_trees_file

    if file_is_empty(consensus_treefile_withBL) or replace is True:
        print("generating consensus tree")

        # define tmps
        consensus_treefile_withBL_tmp = "%s.tmp.nw"%consensus_treefile_withBL

        # get consensus tree (by default the --sup-min is 0, meaning that it is an extended consensus)
        run_cmd("iqtree -con -t %s --sup-min 0.0"%all_trees_file, env=EnvName)

        # add the branch lengths
        get_consensus_tree_with_branchLengths_R = "%s/get_consensus_tree_with_branchLengths.R"%(CWD)
        run_cmd("%s --input_consensus_treefile %s --output_consensus_treefile %s --all_trees_file %s"%(get_consensus_tree_with_branchLengths_R, consensus_treefile, consensus_treefile_withBL_tmp, all_trees_file), env=EnvName)

        # check that the trees are equal
        consensus_tree_withBL = Tree(consensus_treefile_withBL_tmp)
        consensus_tree = Tree(consensus_treefile)
        rf_distance = consensus_tree.robinson_foulds(consensus_tree_withBL)[0]
        if rf_distance!=0: raise ValueError("The get_consensus_tree_with_branchLengths_R changed the topology of the tree")

        # keep
        remove_file(consensus_treefile)
        os.rename(consensus_treefile_withBL_tmp, consensus_treefile_withBL)

    return consensus_treefile_withBL, all_trees_file


###################
