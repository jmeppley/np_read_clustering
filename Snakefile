"""

This workflow attempts to find clusters of similar sequences in a set of
nanopore reads. It is intended to be used for samples that have been size
selected for viral-like particles.

## Overview:
### initial clustering
 * split reads into overlapping size windows
 * all-v-all minimap2 in each window to generate a distance matrix
 * run MCL on concatenated distance matrices

### refine clusters
 * run all-v-all lastal in each good cluster
 * cluster with MCL

### polish clusters
 * with medaka and racon

## Running
### Prerequisites
Create a conda environment with the provided config file (env.yaml)

    mamba env create -f conda/snake.yaml -p ./env

### Configure
The main point of configuration is the working directory: work_dir

The only stricltly necessary input is:

  * all_fasta: a fasta file of all the nanopore reads. This is assumed to be in {work_dir}/all.reads.fasta.

You may also want to specify:

 * work_dir: location to create all files (defaults to 'np_clustering')
 * name: naming prefix for the final sequences (defaults to 'SEQ')
 * pfam_hmm_path: the PFAM HMM file (for gene annotation, pfam annotations are empyt otherwise)

Specify these and any other optional parameters on the command line with --configure key=value key=value ...

or with a configuration file (YAML or JSON) (see the snakemake documentation)

### Snakemake
Run the workflow with snakemake. For example:

    snakemake -s np_reads_clustering --config name=HOT_319 work_dir=clustering -j $THREADS -p --use-conda --conda-frontend mamba

## Outputs of note

From the first pass:
 * A pdf with (pages and pages) of cluster size histograms
 * A table of cluster stats

From the second pass:
 * A pdf of subcluster histograms for each cluster
 * a table of subcluster stats for each cluster
 * A pdf with a synteny plot for each good subcluster showing aligned gene locations and
   PFAM annotations

From polishing:
 * a final fasta file for each subcluster
 * comparison tables of predicted gene lengths for each step in polishing

Final outputs:
 * a fasta file of polished (and renamed) sequences
 * a table of stats for these seqeunces
 * a summary stats file with: number of clusters, number of subclusters, etc

"""

# print help for a common mistake
try:
    from Bio import SeqIO
    import pandas
except ImportError:
    print("""
ERROR: biopython is not installed!

Make sure you are running this workflow using a conda environment with
snakemake, pandas, and biopthon installed. We recommend using the definition file
./conda/snake.yaml with conda:

    $ conda env create -p ./env -f ./conda/snake.yaml

And make sure you activate the environment with:

    $ conda activate ./env

""")
    sys.exit(10)

### Global Params
MCL_I = "{:0.1f}".format(config.get('mcl_i', 5.0))
NAME = config.get('name', 'SEQ')
MIN_POL_READS = config.get('min_pol_reads', 9)

# just do one group (to split big DAGs)
GROUP = config.get('group', -1)
if GROUP >= 0:
    logger.warning("Processing group: {}".format(GROUP))

### File Locations and templates
WORK_DIR = config.get('work_dir', 'np_clustering')
logger.debug("Working directory is: " + WORK_DIR)

MM_THREADS = config.get('mapping_threads', 9)
MCL_THREADS = config.get('mcl_threads', 19)
POLISH_THREADS = config.get('polishing_threads', 9)
HMMER_THREADS = config.get('hmmer_threads', 4)

## STEP 1: windows -> minimap2 -> mcl
ALL_FASTA = config.get('all_fasta', f"{WORK_DIR}/all.reads.fasta")
CLUSTER_OUT = f"{WORK_DIR}/mcl_all/all.I{MCL_I}.mcl"
CLUSTER_STATS=f'{WORK_DIR}/mcl_all/cluster_stats.tsv'

include: "rules/Snakefile.minimap"


## STEP 2: filter -> lastal -> mcl
include: "rules/Snakefile.refine"
# These are templates, since it's one per cluster
CLST_REFINE_DIR = f"{WORK_DIR}/refine_lastal/group.{{group}}/cluster.{{cluster}}"
SUBCLUSTER_TEMPLATE = CLST_REFINE_DIR + f"/self.m8.gt{MFRAC_CUTOFF}.I{MCL_I}.mcl"



## STEP 3: inspect and filter
SUBCLUSTER_READS_TEMPLATE = CLST_REFINE_DIR + "/subclusters/reads/subcluster.{subcluster}.fasta"
SUBCLUSTER_STATS_TEMPLATE = CLST_REFINE_DIR + "/subclusters/cluster_stats.tsv"

include: "rules/Snakefile.subclusters"


## STEP 4: polish
SUBCLUSTER_DIR_TEMPLATE = CLST_REFINE_DIR + '/subclusters/subcluster.{subcluster}'

include: "rules/Snakefile.polish"

## STEP 5: compile stats and sequences of polished seqs into a final table and fasta file
REPORT_FILE = f"{WORK_DIR}/final/report.txt"
include: "rules/Snakefile.finish"

### Entry points

## the whole enchilada
rule finish:
    input:
        report=REPORT_FILE,
        polished=lambda w: get_polished_output_file_names(),
        clusters_pdf=f'{WORK_DIR}/mcl_all/cluster_plots.pdf',

## STEP 1: windows -> minimap2 -> mcl
rule step_1:
    input:
        stats=CLUSTER_STATS,
        pdf=f'{WORK_DIR}/mcl_all/cluster_plots.pdf',

## STEP 2: filter -> lastal -> mcl
rule step_2:
    input: lambda w: get_subcluster_mcl_files(SUBCLUSTER_TEMPLATE)

## STEP 3: inspect and filter
rule step_3:
    input: lambda w: get_subcluster_mcl_files(SUBCLUSTER_STATS_TEMPLATE)

## STEP 4: polish
rule step_4:
    input:
        polished=lambda w: get_polished_output_file_names(),
        clusters_pdf=f'{WORK_DIR}/mcl_all/cluster_plots.pdf',
