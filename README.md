# np_read_clustering

This workflow attempts to find clusters of similar sequences in a set of
nanopore reads. It is intended to be used for samples that have been size
selected for viral-like particles.

## Overview

### Clustering 1: windowed minimap2

The first pass uses minimap2 and mcl to find clusters of similar sequences.

Because an all-v-all comparison is system-taxing and we're looking for clusters
of complete sequences of the same length, we split the reads into buckets by
size. We use overlapping size windows so that we can combine the distance
measures and run MCL on all the reads at once.

### Clustering 2: lastal
Clusters with enough reads are processed with lastal and mcl to
find subclusters. A cutoff of 85% identity is used to build the MCL distance
matrix

### Polishing
Subclusters with enough reads and a tight size distribution are polished with
racon and medaka to produce consensus sequences.

## Installation

The workflow is executed via snakemake, which can handle the installation of
all required dependencies (with the help of conda). 

NOTE: The commands in this section (Installation) are assumed to be run from the repo directory.

### Conda
The simplest approach is to install
[conda](https://docs.conda.io/en/latest/miniconda.html) if you don't have it,
and then to use conda to install snakemake. 

#### Mamba
Mamba is an add on for conda that can isntall programs much faster than conda.
We recommend installing that, too, but you don't have to. If you do have mamba
install, replace "conda" with "mamba" in the following command.

### The snakemake env
Create a conda environment for running snakemake by using the provided
configuration file:

    $ conda env create -p ./nprc.env -f np_read_clustering/conda/snake.yaml

To use snakemake, you'll have to activate the environment in the shell (or
script) from which you want to launch the workflow:

    $ conda activate ./nprc.env

NOTE: The `-p ./nprc.env` option creates the environment in a folder named `nprc.env` in your
current directory. You can use any name and location you wish. You can also use
`-n` to name the environment and keep it in you conda installation location.
See the conda documentation for how to name and activate environments for more
detail.

### A Test Run
That's it. Now you are ready to go. Test your setup (and pre-install the rest
of the dependencies):

    $ snakemake --configfile=config.yaml -j 2 -p --use-conda --conda-frontend mamba

## Running
### Configuration
All of the configuration parameters are top level, so they can be supplied on
the command line as in the test example above.

#### Configfiles
You can also collect all the options in a config file. See the provided example and the snakeamke documentation for
more information on configurations files, but the basic ideas are: (1) they can be formatted as JSON or YAML
and (2) tell snakemake where to find the file with `--configfile=`

#### Required Parameters
The only stricltly necessary input is:

  * all_fasta: a fasta file of all the nanopore reads. This is assumed to be in {work_dir}/all.reads.fasta.

You may also want to specify:

 * work_dir: location to create all files (defaults to 'np_clustering')
 * name: naming prefix for the final sequences (defaults to 'SEQ')
 * pfam_hmm_path: the PFAM HMM file for gene annotation (pfam annotations are empty otherwise)

#### Other Parameters
See the example config.yaml for the rest of the parameters and their defaults

### parallelization and performance
#### Single node
The `-j` flag tells snakemake how many threads are availabe on your computer,
and it will run workflow steps in parallel as much as is possible (usually).

### multithreaded steps
Some steps in the workflow are mutithreaded (EG: minimap2). You can configure
how many threads these teps get in the configuration (EG: mapping_threads).

