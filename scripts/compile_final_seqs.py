"""
Compiles the polished sequences into one file and collects some stats along the
way

From Snakemake:
    input: lambda w: get_polished_comparison_files(kind='gene')
    output:
        fasta=f"{WORK_DIR}/final/polished.seqs.fasta",
        faa=f"{WORK_DIR}/final/polished.seqs.faa",
        stats=f"{WORK_DIR}/final/polished.seqs.stats.tsv"
"""
import os, re, pandas
from Bio import SeqIO

def main(input, output, params):
    stats = {}
    files = {}
    groups = {}
    
    logger.debug("Collecting polished sequences")

    # get fasta files and gene stats
    for sc_comp_file in input:
        sc_dir = os.path.dirname(sc_comp_file)
        group, cluster, subcluster = \
            re.search(r'group.(\d+).+cluster\.(\d+).+subcluster\.(\d+)', sc_dir).groups()
        cl_sc_id = f"{cluster}_{subcluster}"
        files[cl_sc_id] = dict(
            fasta=sc_dir + "/medaka.fasta",
            faa  =sc_dir + "/medaka.faa"
        )
        gene_stats = \
            pandas.read_csv(sc_dir + "/medaka.v.drafts.gene.lengths",
                            sep='\t',
                            index_col=0) \
                  .loc['medaka']
        stats.setdefault(cluster, {})[subcluster] = \
            {f"gene_{k}": v for k,v in gene_stats.items()}
        groups[cluster] = group

    logger.debug("Found {len(files)} poilished seqs from {len(stats)} clusters")

    # get read len stats for subclusters
    for cluster in stats:
        sc_tsv = (f"{params.work_dir}/refine_lastal/group.{groups[cluster]}"
                  f"/cluster.{cluster}/subclusters/cluster_stats.tsv")

        sc_stats = pandas.read_csv(sc_tsv, sep='\t')
        for sc, row in sc_stats.iterrows():
            stats[cluster].setdefault(str(sc),{}).update(dict(
                read_len_mean=row['mu'],
                read_len_dev=row['sigma'],
                read_count=row['N']))

    # convert stats to table
    df = pandas.DataFrame({f"{cl}_{sc}": sc_stats
                        for cl, cl_stats in stats.items()
                        for sc, sc_stats in cl_stats.items()},) \
            .T
                     
    # polsihed fasta
    with open(str(output.fasta), 'wt') as fasta_out:
        for cl_sc_id in files:
            for read in SeqIO.parse(files[cl_sc_id]['fasta'], 'fasta'):
                np_read = read.id
                read.id = f"{params.name}_{cl_sc_id}"
                df.loc[cl_sc_id, 'length'] = len(read)
                N = df.loc[cl_sc_id, 'read_count']
                read.description = f"n_reads={N};rep_read={np_read}"
                fasta_out.write(read.format('fasta'))

    # write out stats table
    df.to_csv(str(output.stats), sep='\t')

    with open(str(output.faa), 'wt') as faa_out:
        for cl_sc_id in files:
            N = 0
            for gene in SeqIO.parse(files[cl_sc_id]['faa'], 'fasta'):
                N += 1
                read.id = f"{params.name}_{cl_sc_id}_{N}"
                faa_out.write(read.format('fasta'))

if __name__ == "__main__":
    main(snakemake.input, snakemake.output, snakemake.params)
