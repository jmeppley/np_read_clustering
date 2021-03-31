#!/usr/bin/env python
"""
Given:
 * cluster reads fasta
 * cluster reads faa
 * cluster reads faa PFAM hits (non-overlapping)
 * cluster lastal results (aggregated)
 * subcluster mcl results
 * output locations (fasta dir + pdf file names)
 * params:
  * min_sub_size (10)
  * sigma_cutoff (None)
  * max_colored_genes (8)
  * gene_cmap ('gnuplot2')
  * max_synteny_reads (15)
  * read_cmap ('cool')
 
Generate 2 pdf files:
 * histograms of read lengths and pctids for subclusters
 * read synteny plots for top reads of each subcluster
 
Generate a fasta file of reads for each subcluster
"""

import numpy
import os
import pandas
import re

from collections import Counter, defaultdict
from functools import partial
from itertools import zip_longest

from scipy import stats
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt, cm, colors
from matplotlib.patches import Polygon
from matplotlib.backends.backend_pdf import PdfPages
from snakemake.rules import Namedlist
from Bio import SeqIO

from jme.jupy_tools.hit_tables import parse_blast_m8, BLAST_PLUS
from edl import blastm8

BLACK = (0, 0, 0, 1)

def grouper_trim(iterable, n):
    "Collect data into fixed-length chunks or blocks and trim last chunk (and all null values)"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return ([i for i in group if i is not None] for group in zip_longest(*args, fillvalue=None))

def main(input, output, params):
    """
    input output should be namedlists (from snakemake)
    params should be a dict (so we can fall back to defaults)
    """

    # prodigal amino acid output 
    faa_file = input.faa
    # cluster fasta reads
    fasta_file = input.fasta
    # PFAM results (non-overlapping)
    dom_tbl_U = input.domtbl
    # mcl results
    mcl_file = input.mcl
    # lastal table raw
    lastal_file = input.lastal
    # lastal table aggregated
    agg_file = input.agg
    
    ## load clusters (just a list of reads in each cluster, sorted by size)
    subclusters = load_clusters(mcl_file, params.get('min_sub_size', 10))
    
    ## load the fasta, keeping dict of lengths
    cluster_reads = {r.id:r for r in SeqIO.parse(fasta_file, 'fasta')}
    # use just the first word in the read id as a short name
    read_lens = {r.id.split('-')[0]:len(r) for r in cluster_reads.values()}

    ## plot all cluster hitsts, applying sigma cutoff
    subcluster_ids = plot_cluster_hists(subclusters, read_lens, agg_file, output.stats, output.hist_pdf, params)
    
    ## make the synteny plots for each good subcluster
    plot_subcluster_synteny(subcluster_ids, subclusters, read_lens, lastal_file, faa_file, dom_tbl_U, output.gene_pdf, params)
    
    ## write out sub cluster fasta
    os.makedirs(str(output.fasta_dir), exist_ok=True)
    for subcluster_id in subcluster_ids:
        with open(str(output.fasta_dir) + f"/subcluster.{subcluster_id}.fasta", 'wt') as fasta_out:
            for read_id in subclusters[subcluster_id]:
                fasta_out.write(cluster_reads[read_id].format('fasta'))
    
def plot_cluster_hists(subclusters,
                       read_lens,
                       agg_file,
                       stats_file,
                       pdf_file,
                       params
                      ):
    """
    For each subcluster plot:
     * histogram aod all read-read mfracs
     * histogram of all read lengths with overlaid normal dist
    """
    
    # open PDF file
    pdf = PdfPages(pdf_file)
    
    mx_len = max(read_lens.values())
    mn_len = min(read_lens.values())
    window = [mn_len, mx_len]

    sigma_cutoff = params.get('sigma_cutoff', -1)
    
    # first pass to chose subclusters to keep and plot
    cluster_stats = {}
    for i, subcluster in enumerate(subclusters):
        keep = True
        if len(subcluster) < params.get('min_sub_size', 10):
            break
        
        # calculate best normal fit to length dist
        cluster_lens = numpy.array([read_lens[r.split('-')[0]] for r in subcluster])
        counts, bins = numpy.histogram(cluster_lens, bins=100, range=window)
        #from scipy import stats
        mu, sigma = stats.norm.fit(cluster_lens)

        if sigma_cutoff > 0 and sigma > sigma_cutoff:
            keep = False
        
        # calculate the stats
        X = numpy.array([numpy.mean((bins[i], bins[i-1])) for i in range(1,len(bins))])
        tot_in, tot_out, n_in, n_out = numpy.zeros(4)
        for x, count in zip(X, counts):
            if x < mu - sigma or x > mu + sigma:
                tot_out += count
                n_out += 1
            else:
                tot_in += count
                n_in += 1
        mean_in = tot_in / n_in
        mean_out = tot_out / n_out if n_out > 0 else 0
        ratio = mean_in / mean_out
        n_ratio = n_in / (n_out + n_in)

        cluster_stats[i] = dict(zip(
            ['mu', 'sigma', 'ratio', 'n_ratio', 'N', 'keep', 'counts', 'bins', 'X'],
            [mu, sigma, ratio, n_ratio, len(subcluster), keep, counts, bins, X]
        ))

    # build cluster stats table
    write_cols = ['mu', 'sigma', 'ratio', 'n_ratio', 'N', 'keep']
    cl_st_table = pandas.DataFrame([[i,] + [d[k] for k in write_cols] 
                                    for i,d in cluster_stats.items()],
                                   columns=['index'] + write_cols)
    # write stats to file
    cl_st_table.to_csv(stats_file, sep='\t', index=None)

    # pull out list of good subclusters
    subcluster_ids = list(cl_st_table.query('keep').index)

    # load agg hits
    agg_table = pandas.read_csv(agg_file, sep='\t')
        
    # max 8 per page
    mx_rows = 8
    for page_sc_ids in grouper_trim(cluster_stats.keys(), mx_rows):
        N = len(page_sc_ids)
        fig, axes = plt.subplots(N, 4, figsize=[11 * N / mx_rows, 8.5], sharey="col", sharex="col", squeeze=False)
        fig.subplots_adjust(hspace=.7, wspace=.6)

        ax_rows = iter(axes)
        for i, subcluster_id in enumerate(page_sc_ids):
             
            axs = next(ax_rows)
            
            # remove axes from top and right
            for ax in axs:
                for side in ['top', 'right']:
                    ax.spines[side].set_visible(False)
            
            ax_sc_mf, ax_sc_id, ax_h_mf, ax_h_ln = axs
            
            # get the subset of the agg table for this subcluster
            subcluster = set(subclusters[subcluster_id])
            sub_slice = (agg_table['query'].apply(lambda q: q in subcluster)
                         & agg_table.hit.apply(lambda h: h in subcluster))
            agg_hits_cluster = agg_table[sub_slice] \
                .eval('mean_len = (hlen + qlen) / 2') \
                .eval('frac = mlen / mean_len')
            mfrac_dict = agg_hits_cluster.set_index(['query','hit']).mfrac.to_dict()

            # scatter plot mfrac and mean length
            ax_sc_mf.scatter(agg_hits_cluster.mfrac.values,
                             agg_hits_cluster.mean_len.values,
                             marker='.',
                             alpha=.5
                            )
            ax_sc_mf.set_ylabel ('mean_len')

            # scatter plot of pctid and matched fraction
            ax_sc_id.scatter(agg_hits_cluster.pctid.values,
                             agg_hits_cluster.frac.values,
                             marker='.',
                             alpha=.5
                            )
            ax_sc_id.set_ylabel ('frac aln')


            # plot hist of pairwise mfracs
            h = ax_h_mf.hist(get_mfracs(subcluster, mfrac_dict=mfrac_dict), bins=100, range=[50,100])

            # plot hist of read lens
            sc_stats = cluster_stats[subcluster_id]
            counts = sc_stats['counts']
            X = sc_stats['X']

            # recreate histogram from counts and X
            ax_h_ln.bar(X, counts, color='blue')

            # overlay norm dist
            best_fit_line = stats.norm.pdf(X, sc_stats['mu'], sc_stats['sigma'])
            best_fit_line = best_fit_line * counts.sum() / best_fit_line.sum()
            p = ax_h_ln.plot(X, best_fit_line, color='red', alpha=.5)
            
            
            ax_h_mf.set_ylabel(f"s.cl: {subcluster_id}")
            ax_h_ln.set_ylabel(f"{len(subcluster)} {int(sc_stats['sigma'])}")
            
            if i == N - 1:
                xl = ax_sc_mf.set_xlabel("score")
                xl = ax_h_ln.set_xlabel("length")
                xl = ax_sc_id.set_xlabel ('match %ID')
                xl = ax_h_mf.set_xlabel ('score')

        # close plot and go to next pdf page
        pdf.savefig(bbox_inches='tight')
        plt.close()
        
    pdf.close()

    # save stats to file, but drop extra data first
    write_cols = ['mu', 'sigma', 'ratio', 'n_ratio', 'N']
    pandas.DataFrame([[i,] + [d[k] for k in write_cols] 
                      for i,d in cluster_stats.items()],
                     columns=['index'] + write_cols).to_csv(stats_file, sep='\t', index=None)

    return subcluster_ids

def get_N_colors(N, cmap_name='Dark2'):
    """ given N and a colormap, get N evenly spaced colors"""
    color_map=plt.get_cmap(cmap_name)
    return [color_map(c) for c in numpy.linspace(0, 1, N)]

def get_scaled_color(value, minv=0, maxv=1, alpha=.75, reverse=False, cmap_name='cool'):
    colormap = plt.get_cmap(cmap_name)
    if reverse:
        maxv, minv = minv, maxv
    rangev = maxv - minv
    color = colormap((value - minv) / rangev)
    return color[:3] + (alpha,)

def get_mfracs(reads, mfrac_dict):
    return [mfrac_dict.get((r1, r2), 0)
            for r1 in reads
            for r2 in reads
            if r2 > r1
           ]

def plot_subcluster_synteny(subcluster_ids,
                            subclusters,
                            read_lens,
                            lastal_file,
                            faa_file,
                            dom_tbl_U,
                            pdf_file,
                            params
                           ):
    """
    For each subcluster:
     * identify the N genes that appear in the most reads
     * identify the M reads that have the most of the top genes
     * plot
    """
    
    ## load the gene annotations
    # first get positions from faa headers
    read_genes = {}
    for gene in SeqIO.parse(faa_file, 'fasta'):
        gene_id, start, end, strand, _ = [b.strip() for b in gene.description.split("#")]
        read, name, gene_no = re.search(r'^((\w+)-[^_]+)_(\d+)', gene_id).groups()

        read_genes.setdefault(name, []).append(dict(
            gene_id=gene_id,
            start=int(start),
            end=int(end),
            strand=int(strand),
            num=int(gene_no),
            pfam=None,
        ))

    # convert to dict of DataFrames from dict of lists of dicts
    read_genes_tables = {read:pandas.DataFrame(genes).set_index('gene_id')
                         for read, genes in read_genes.items()}
    
    # and add PFAM annotations
    for read, hits in blastm8.generate_hits(dom_tbl_U, format='hmmsearchdom'):
        read_id = read.split("-")[0]
        read_genes_table = read_genes_tables[read_id]
        
        for hit in hits:
            gene_id = hit.read

            # only assign PFAm if it's the first hit for the gene
            if pandas.isna(read_genes_table.loc[gene_id, 'pfam']):
                pfam = hit.hit
                read_genes_table.loc[gene_id, 'pfam'] = pfam

    # load all the read to read hits
    read_hits = parse_blast_m8(lastal_file, format=BLAST_PLUS)
                                
    # now open the PDF file
    pdf = PdfPages(pdf_file)

    
    # for each good subcluster
    for subcluster_id in subcluster_ids:
        subcluster = set(subclusters[subcluster_id])
        subcluster_names = {r.split('-')[0]:r for r in subcluster}

        fig = plot_subcluster_genes(subcluster_id, subcluster_names, read_genes_tables, read_hits, read_lens, params)
        
        # close plot and go to next pdf page
        pdf.savefig(bbox_inches='tight')
        plt.close()
        
    pdf.close()
    
def plot_subcluster_genes(subcluster_id, subcluster_names, read_genes_tables, read_hits, read_lens, params):
    """
    make a plot of gene positions:
     ax1 has a scatter plot of mean position by pfam
     ax2 has aligned genomes with top pfams colored
    """

    # get the positions of the named PFAMs
    pf_positions = defaultdict(list)
    for read, gene_table in read_genes_tables.items():
        if read in subcluster_names:
            # do we want to flip the read dir? (too many strand < 1)
            reverse = gene_table.eval('glen = strand * (end - start)').glen.sum() < 1
            for start, end, pfam in gene_table[['start','end','pfam']].values:
                if pandas.isna(pfam):
                    continue
                if reverse:
                    start, end = [read_lens[read] - p for p in (start, end)]
                # add mean post to list for this pfam
                pf_positions[pfam].append((end + start) / 2)

    # chose which genes to color
    N = params.get('max_colored_genes', 8)
    sorted_genes = sorted(pf_positions.keys(), key=lambda k: len(pf_positions[k]), reverse=True)
    top_N_pfams = sorted_genes[:N]
    gene_color_dict = dict(zip(top_N_pfams, get_N_colors(N, cmap_name=params.get('gene_cmap', 'Dark2'))))


    # chose which reads to draw
    M = params.get('max_synteny_reads', 20)
    def count_top_pfams_in_read(read):
        if read in read_genes_tables:
            return sum(1 for p in read_genes_tables[read].pfam.values
                       if p in top_N_pfams)
        return 0
    top_M_reads = sorted(subcluster_names, 
                         key=count_top_pfams_in_read,
                         reverse=True,
                        )[:M]
    m = len(top_M_reads)

    # calculate the sizes necessary to draw genes using the matplotlib arrow function
    align_height = (7 * (m-1) / (M-1)) #use up to 7 in
    figsize = [8.5, 4 + align_height]
    
    fig, axes = plt.subplots(2,1, figsize=figsize, gridspec_kw={'height_ratios':[4,align_height]}, sharex='col')
    fig.subplots_adjust(hspace=.1,)

    ## draw gene positions
    ax = axes[0]

    ax.set_title(f'PFAM annotations in subcluster {subcluster_id}')
    
    n = params.get('max_plotted_genes', 18)
    sorted_pf = sorted([p for p in sorted_genes[:n] if len(pf_positions[p]) > 1], 
                       key=lambda p: numpy.mean(list(pf_positions[p])))
    for i, p in enumerate(sorted_pf):
        x,y = zip(*((gp,i) for gp in pf_positions[p]))
        ax.scatter(x,y, 
                   c=len(y) * [gene_color_dict.get(p, BLACK)], 
                   ec=None, alpha=.5)
    yt = ax.set_yticks(range(len(sorted_pf)))
    ytl = ax.set_yticklabels(sorted_pf)
    for label in ytl:
        label.set_color(gene_color_dict.get(label.get_text(), BLACK))
    
    ## draw alignments
    ax = axes[-1]
    min_x = 0
    max_x = max(read_lens[r] for r in subcluster_names)
    range_x = max_x - min_x
    range_y = M
    thickness = .5
    head_length = range_x * (thickness / range_y) * (figsize[1] / figsize[0])

    cmap = params.get('read_cmap','cool')

    min_pctid = read_hits.pctid.min()
    pctid_range = 100 - min_pctid

    get_conn_color = partial(get_scaled_color, minv=min_pctid, maxv=100, alpha=.75, cmap_name=cmap)

    y = 0
    pad = .1
    prev_read = None
    for name in top_M_reads:
        read = subcluster_names[name]
        read_length = read_lens[name]
        if name in read_genes_tables:

            gene_table = read_genes_tables[name]

            # do we want to flip the read dir? (too many strand < 1)
            reverse = gene_table.eval('glen = strand * (end - start)').glen.sum() < 1

            # draw genes
            for start, end, strand, pfam in gene_table[['start','end','strand','pfam']].values:
                if reverse:
                    strand = -1 * strand
                    start = read_length - start
                    end = read_length - end

                strand = int(strand)
                hl = min(head_length, end-start)
                al = max((end - start) - hl, .0001) * strand
                ast = start if al > 0 else end
                color = gene_color_dict.get(pfam, 'k')
                plt.arrow(ast, y, al, 0, fc=color, ec=color, 
                          lw=0,
                          width=thickness, head_width=thickness, 
                          head_length=hl, 
                          head_starts_at_zero=(int(strand) > 0))
        else:
            reverse=False

        # connect matched segments for read pairs
        if prev_read is not None:
            # get hits between reads
            pair_hits = read_hits.query(f'(hit == "{read}" and query == "{prev_read}") or '
                                        f'(query == "{read}" and hit == "{prev_read}")') \
                                 .query('hit != query') \
                                 .sort_values('score', ascending=True)
            # loop over hits
            cols = ['query', 'hit', 'qstart', 'qend', 'hstart', 'hend', 'pctid']
            for query, hit, qstart, qend, hstart, hend, pctid in pair_hits[cols].values:
                # if hit was recorded the other way, flip hit/query
                if query == prev_read:
                    qstart, qend, hstart, hend = hstart, hend, qstart, qend

                # if either read is reversed, flip x coordinates 
                if reverse:
                    qstart = read_length - qstart
                    qend = read_length - qend
                if prev_rev:
                    hstart = prev_len - hstart
                    hend = prev_len - hend

                # draw connecting paralellogram
                color = get_conn_color(pctid, alpha=.9)
                xy = numpy.array([(hstart, y-1+pad),
                                  (qstart, y-pad),
                                  (qend, y-pad),
                                  (hend, y-1+pad)])
                ax.add_patch(Polygon(xy, fc=(.6,.6,.6,.2), ec=color))   

        # save read info for next one
        prev_read = read
        prev_rev = reverse
        prev_len = read_length

        # increment y value
        y += 1

    x = plt.xlim(min_x - 50, max_x + 50)
    y = plt.ylim(-.5, y - .5)

    plt.yticks(list(range(m)), top_M_reads)
    plt.xlabel('read position')

    cax = plt.axes([0.95, 0.15, 0.025, 0.4 * (align_height / 7)])
    plt.colorbar(mappable=cm.ScalarMappable(norm=colors.Normalize(min_pctid, 100), cmap=cmap), cax=cax)    
    cl = cax.set_ylabel('alignment %ID')

    return fig
    
def load_clusters(mcl_file, size_cutoff=10):
    with open(mcl_file) as mcl_lines:
        return [c for c in [line.strip().split() for line in mcl_lines] if len(c) >= size_cutoff]

# scriptify
if __name__ == "__main__":
    try:
        # assume this is called from snakemake
        input = snakemake.input
        output = snakemake.output
        params = dict(snakemake.params.items())
    except NameError:
        # TODO: fallback to argparse if we call from the command line (for testing)
        import argparse
        raise Exception("Currently only works from snakemake, sorry")

    main(input, output, params)
        
