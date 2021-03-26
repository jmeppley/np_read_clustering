"""
process all the mcl clusters:
 * find reads in window fasta files
 * check length distribution
 * if it passes, write out fasta file

    input:
        mcl=f'{output_dir}/mcl_all/all.I{mcl_i}.mcl',
        fasta=all_fasta
        read_lens=f'{WORK_DIR}/all.reads.lengths.tsv'
    output:
        reads=directory(f'{output_dir}/mcl_all/reads'),
        stats=f'{output_dir}/mcl_all/cluster_stats.tsv',
    params:
        sigma_cutoff=sigma_cutoff_pre,
        min_cl_size=MIN_POL_READS + 1,
"""
import pandas, numpy, os
from collections import deque
from itertools import cycle
from scipy import stats
from Bio import SeqIO

# load the read lengths from the summary file
read_lens = pandas.read_csv(snakemake.input.read_lens,
                            sep='\t', 
                            names=['read_id','sequence_length_template'], 
                            index_col='read_id', 
                            header=None).sequence_length_template.to_dict()


# process clusters to choose keepers
cluster_data = []
read_clusters = {}
sigma_cutoff = snakemake.params.sigma_cutoff
count_cutoff = snakemake.params.min_cl_size

# loop over clusters in mcl_file
with open(str(snakemake.input.mcl)) as mcl_lines:
    for i,line in enumerate(mcl_lines):
        
        # get cluster read names
        reads = set(line.strip().split())
        count = len(reads)
        
        # get cluster read length dist
        cluster_lens = numpy.array([read_lens[r] for r in reads])
        counts, bins = numpy.histogram(cluster_lens, bins=100)
        X = numpy.array([numpy.mean((bins[j], bins[j-1])) for j in range(1,len(bins))])
        mu, sigma = stats.norm.fit(cluster_lens)

        keep = (sigma <= sigma_cutoff and count >= count_cutoff)
        cluster_data.append(dict(num=i, count=count, sigma=sigma, mu=mu,
                                 keep=keep))

        if keep:
            """
            # write read list
            if not os.path.exists(str(snakemake.output.reads)):
                os.makedirs(str(snakemake.output.reads), exist_ok=True)
            with open(f"{output.reads}/cluster.{i}.reads", 'wt') as reads_out:
                reads_out.write("\n".join(reads) + "\n")
            """
            # save cluster id
            for read in reads:
                read_clusters[read] = i

cluster_table = pandas.DataFrame(cluster_data)

## assign groups
# this serves 2 purposes: 
#  1) we limit the number of files in each folder (too many files can slow
#     down snakemake)
#  2) we enable running the workflow in chunks (can perform better in some
#     cases)

keepers = cluster_table.query('keep')
num_keepers = keepers.shape[0]

# we want the number of groups, but we can get it from group_size
if 'group_size' in snakemake.config and 'num_groups' not in snakemake.config:
    group_size = snakemake.config['group_size']
    n_groups = int(numpy.ceil(num_keepers/group_size))
elif:
    n_groups = snakemake.config.get('num_groups', 100)

# assigne a group to each cluster (round-robin)
groups = cycle(range(n_groups))
cluster_groups = {c:next(groups) for c in keepers['num']}
cluster_table['group'] = [cluster_groups.get(c,None)
                          for c in cluster_table['num']

# save cluster table
cluster_table.to_csv(str(snakemake.output.stats), sep='\t',
                                      index=False)

# write fasta files
if not os.path.exists(str(snakemake.output.reads)):
    os.makedirs(str(snakemake.output.reads), exist_ok=True)

# limit number of open files with
n_open = 250
open_handle_ids = deque([])
handles = {}
def open_cluster_fasta(i):
    """ 
    checks for open handle for this scluster and returns it if found

    otherwise closes oldest handle and replaes with new handle for this cluster
    """
    # return open handle if it exists
    try:
        return handles[i]
    except KeyError:
        pass
    
    # close handle(s) if we have too many
    while len(handles) > n_open - 1:
        # drop oldest
        drop_id = open_handle_ids.popleft()
        
        # close and delete
        handles[drop_id].close()
        del handles[drop_id]
        
    group = cluster_groups[i]
    fasta_file = f"{snakemake.output.reads}/group.{group}/cluster.{i}.fasta"
    fd = os.path.dirname(fasta_file)
    if not os.path.exists(fd):
        os.path.makedirs(fd)
    handle = open(fasta_file, 'at')
    handles[i] = handle
    open_handle_ids.append(i)
    return handle
    
# loop over all reads and write out
for read in SeqIO.parse(snakemake.input.fasta, 'fasta'):
    try:
        cluster = read_clusters[read.id]
    except KeyError:
        # skip if no cluster
        continue
    
    open_cluster_fasta(cluster).write(read.format('fasta'))
                
