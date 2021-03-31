"""
generate a short report with:
 * number of size windows
 * number of clusters (and # > 10 and number kept)
 * number of subclusters (and # > 10 and number kept)
 * number of polished seqs (with length stats)
 
From snakemake:
    input:
        pol_stats=f"{WORK_DIR}/final/polished.seqs.stats.tsv",
        mcl_stats=f'{WORK_DIR}/mcl_all/cluster_stats.tsv',
        sc_mcls=get_sc_mcls_files,
    output:
        report=REPORT_FILE,
    params:
        min_cl_size=MIN_POL_READS + 1,
        sigma_cutoff=SIGMA_CUTOFF
"""
import pandas, numpy

def main(input, output, params):
    with open(str(output.report), 'wt') as output_handle:
        cluster_stats = pandas.read_csv(str(input.mcl_stats), sep='\t',
                                        index_col=0)
        n_clusters = cluster_stats.shape[0]
        n_gt_size_cutoff = cluster_stats.query(f'count >= {params.min_cl_size}').shape[0]
        n_kept = cluster_stats.query('keep').shape[0]

        output_handle.write(f"Cluster Search Results:\n"
                            f"  minimap2 clusters:\n"
                            f"    clusters: {n_clusters}\n"
                            f"    gt_{params.min_cl_size}: {n_gt_size_cutoff}\n"
                            f"    kept_clusters: {n_kept}\n\n")

        # count raw subclusters in mcl files
        n_scs, n_sc_gt_cutoff = 0, 0
        for sc_mcl_file in input.sc_mcls:
            with open(sc_mcl_file) as mcl_lines:
                for line in mcl_lines:
                    n_scs += 1
                    if len(line.strip().split()) > params.min_cl_size:
                        n_sc_gt_cutoff += 1

        # get stats from polished subclusters
        pol_stats = pandas.read_csv(str(input.pol_stats), sep='\t', index_col=0)
        pol_lens = pol_stats[pol_stats['length'].notna()]['length'].values
        n_sc_kept = len(pol_lens)

        output_handle.write(f"  lastal subclusters:\n"
                            f"    subclusters: {n_scs}\n"
                            f"    gt_{params.min_cl_size}: {n_sc_gt_cutoff}\n"
                            f"    kept_subclusters: {n_sc_kept}\n\n")
    
        output_handle.write(f"  polished seqs:\n"
                            f"    count: {len(pol_lens)}\n"
                            f"    mean: {pol_lens.mean()}\n"
                            f"    max: {pol_lens.max()}\n"
                            f"    min: {pol_lens.min()}\n"
                            f"    median: {numpy.median(pol_lens)}\n"
                            f"    stddev: {pol_lens.std()}\n\n")


if __name__ == "__main__":
    main(snakemake.input, snakemake.output, snakemake.params)

