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
Clusters with a tight size distribution are processed with lastal and mcl to
find subclusters. A cutoff of 85% identity is used to build the MCL distance
matrix

### Polishing
Subclusters with enough reads and a tight size distribution are polished with
racon and medaka to produce consensus sequences.
