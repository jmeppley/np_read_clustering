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
"""

configfile: 'config.json'

### Global Params
MCL_I = "{:0.1f}".format(config.get('mcl_i', 5.0))

GROUP_SIZE = snakemake.config.get('group_size', 1000)
# just do one cluster (for debugging?)
CLUSTER = config.get('cluster', -1)
if CLUSTER >= 0:
    GROUP = int(cluster / GROUP_SIZE)
    logger.warning("Processing just one cluster: {}".format(CLUSTER))
else:
    # just do one group (to split big DAGs)
    GROUP = config.get('group', -1)
    if GROUP >= 0:
        logger.warning("Processing group: {}".format(GROUP))

### File Locations and templates
WORK_DIR = config.get('work_dir', 'np_clustering')


## STEP 1: windows -> minimap2 -> mcl
ALL_FASTA = config['all_fasta']
ALL_SUMMARY = config['all_summary']
CLUSTER_OUT = f"{WORK_DIR}/mcl_all/all.I{MCL_I}.mcl"

include: "rules/Snakefile.minimap"


## STEP 2: filter -> lastal -> mcl
# These are templates, since it's one per cluster
CLST_REFINE_DIR = "{WORK_DIR}/refine_lastal/group.{group}/cluster.{cluster}"
SUBCLUSTER_TEMPLATE = CLST_REFINE_DIR + "/self.m8.gt{MFRAC_CUTOFF}.I{MCL_I}.mcl"

include: "rules/Snakefile.refine"


## STEP 3: inspect and filter
SUBCLUSTER_READS_TEMPLATE = CLST_REFINE_DIR + "/subclusters/reads/subcluster.{subcluster}.fasta"
SUBCLUSTER_STATS_TEMPLATE = CLST_REFINE_DIR + "/subclusters/cluster_stats.tsv"

include: "rules/Snakefile.subclusters"


## STEP 4: polish
SUBCLUSTER_DIR_TEMPLATE = CLST_REFINE_DIR + '/subclusters/subcluster.{subcluster}'

include: "rules/Snakefile.polish"

### Entry points

## the whole enchilada
rule finish:
    input: get_polished_comparison_files

## STEP 1: windows -> minimap2 -> mcl
rule step_1:
    input: CLUSTER_OUT

## STEP 2: filter -> lastal -> mcl
rule step_2:
    input: lambda w: get_subcluster_files(SUBCLUSTER_TEMPLATE)

## STEP 3: inspect and filter
rule step_3:
    input: lambda w: get_subcluster_files(SUBCLUSTER_STATS_TEMPLATE)

## STEP 4: polish
# just run finish (above)

