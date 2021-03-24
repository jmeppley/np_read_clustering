"""
process all the mcl clusters:
 * find reads in window fasta files
 * check length distribution
 * if it passes, write out fasta file

    input:
        mcl=f'{output_dir}/mcl_all/all.I{mcl_i}.mcl',
        fasta=all_fasta
    output:
        reads=directory(f'{output_dir}/mcl_all/reads'),
        stats=f'{output_dir}/mcl_all/cluster_stats.tsv',
        pdf=f'{output_dir}/mcl_all/cluster_plots.pdf'
    params:
        sigma_cutoff=sigma_cutoff_pre,
        summary=summary_table
"""
import pandas, numpy, os
from collections import deque
from scipy import stats
from Bio import SeqIO

# load the read lengths from the summary file
read_lens = pandas.read_csv(snakemake.params.summary,
                            sep='\t', 
                            usecols=['read_id','sequence_length_template'], 
                            index_col='read_id', 
                            header=0).sequence_length_template.to_dict()


cluster_data = []
read_clusters = {}
sigma_cutoff = snakemake.params.sigma_cutoff

group_size = snakemake.config.get('group_size', 1000)

def get_group(cluster):
    return int(cluster / group_size)

# loop over clusters in mcl_file
with open(str(snakemake.input.mcl)) as mcl_lines:
    for i,line in enumerate(mcl_lines):
        
        # get cluster read names
        reads = set(line.strip().split())
        
        # get cluster read length dist
        cluster_lens = numpy.array([read_lens[r] for r in reads])
        counts, bins = numpy.histogram(cluster_lens, bins=100)
        X = numpy.array([numpy.mean((bins[j], bins[j-1])) for j in range(1,len(bins))])
        mu, sigma = stats.norm.fit(cluster_lens)

        cluster_data.append(dict(num=i, count=len(reads), sigma=sigma, mu=mu))

        if sigma <= sigma_cutoff:
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
                
# write fasta files
if not os.path.exists(str(snakemake.output.reads)):
    os.makedirs(str(snakemake.output.reads), exist_ok=True)

# limit number of open files with
n_open = 250
open_handle_ids = deque([])
handles = {}
def open_cluster_fasta(i):
    """ """
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
                
#TODO make the PDF at the same time

pandas.DataFrame(cluster_data).to_csv(str(snakemake.output.stats), sep='\t')
