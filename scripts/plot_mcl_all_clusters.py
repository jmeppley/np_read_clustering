"""
plot all the mcl clusters:
    input:
        mcl=f'{output_dir}/mcl_all/all.I{mcl_i}.mcl',
        read_lens=f'{WORK_DIR}/all.reads.lengths.tsv'
        abc=f"{WORK_DIR}/mcl_all/all.abc"
    output:
        pdf=f'{output_dir}/mcl_all/cluster_plots.pdf'
    params:
        sigma_cutoff=sigma_cutoff_pre,
        min_cl_size=MIN_POL_READS + 1,
"""
import pandas, numpy, os
from scipy import stats
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sigma_cutoff = snakemake.params.sigma_cutoff
count_cutoff = snakemake.params.min_cl_size

# load the read lengths from the summary file
read_lens = pandas.read_csv(snakemake.input.read_lens,
                            sep='\t', 
                            names=['read_id','sequence_length_template'], 
                            index_col='read_id', 
                            header=None).sequence_length_template.to_dict()

# load the all.v.all mfrac values (but just map pairs where q > h)
mfrac_dict = {tuple(sorted(i)):m
              for i,m in pandas.read_csv(snakemake.input.abc, sep='\t',
                                         names=['q', 'h', 'mfrac'],
                                         header=None, index_col=['q','h']) \
                               .mfrac.items()}

def get_mfracs(reads, mfrac_dict=mfrac_dict):
    return [mfrac_dict.get((r1, r2), 0)
            for r1 in reads
            for r2 in reads
            if r2 > r1
           ]

# load the clusters
with open(snakemake.input.mcl) as mcl_lines:
    all_clusters = [line.strip().split() for line in mcl_lines]

# plots
pdf=PdfPages(snakemake.output.pdf)

ROWS = 20
COLS = 5
N0=0

while True:
    clusters = all_clusters[N0:N0+ROWS*COLS]
    if len(clusters) == 0:
        break
    rows = int(numpy.ceil(len(clusters)/COLS))

    fig, axes = plt.subplots(rows, COLS*2, figsize=[COLS*4,rows], sharex=False,
                            squeeze=False)
    fig.subplots_adjust(hspace=.7, wspace=.6)

    cluster_iter = enumerate(clusters, start=N0)
    axes_list = axes.flatten()
    for i, cluster in enumerate(clusters):
        j = i*2
        ax1, ax2 = axes_list[j:j+2]

        # plot hist of pairwise mfracs
        h = ax1.hist(get_mfracs(cluster, mfrac_dict=mfrac_dict), bins=100,
                     range=[0,100])

        # plot hist of read lens
        cluster_lens = numpy.array([read_lens[r] for r in cluster])
        counts, bins, h_line = ax2.hist(cluster_lens, bins=100, histtype='step')
        X = numpy.array([numpy.mean((bins[j], bins[j-1])) for j in range(1,len(bins))])
        mu, sigma = stats.norm.fit(cluster_lens)

        # overlay norm dist
        best_fit_line = stats.norm.pdf(X, mu, sigma)
        best_fit_line = best_fit_line * counts.sum() / best_fit_line.sum()
        p = ax2.plot(X, best_fit_line, color='red', alpha=.5)

        keep = (sigma <= sigma_cutoff and len(cluster) >= count_cutoff)
        if keep:
            ax1.set_ylabel('keep')

        ax2.set_ylabel(f"c{i} n={len(cluster)}")
        # only put xlabels on bottom plots
        if i >= len(clusters) - COLS:
            xl = ax1.set_xlabel("score")
            xl = ax2.set_xlabel("length")

        # remove axes from top and right
        for ax in [ax1, ax2]:
            for side in ['top', 'right']:
                ax.spines[side].set_visible(False)

    # hide unused axes (on last page)
    for i in range(j+2, len(axes_list)):
        fig.delaxes(axes_list[i])

    pdf.savefig(bbox_inches='tight')
    plt.close()

    N0 += ROWS * COLS
    if len(cluster) < count_cutoff:
        break
pdf.close()

