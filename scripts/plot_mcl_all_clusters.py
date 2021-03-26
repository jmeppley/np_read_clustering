"""
plot all the mcl clusters:
    input:
        mcl=f'{output_dir}/mcl_all/all.I{mcl_i}.mcl',
        read_lens=f'{WORK_DIR}/all.reads.lengths.tsv'
    output:
        pdf=f'{output_dir}/mcl_all/cluster_plots.pdf'
    params:
        sigma_cutoff=sigma_cutoff_pre,
        min_cl_size=MIN_POL_READS + 1,
"""
import pandas, numpy, os
from matplotlib import pyplot as plt

# load the read lengths from the summary file
read_lens = pandas.read_csv(snakemake.input.read_lens,
                            sep='\t', 
                            names=['read_id','sequence_length_template'], 
                            index_col='read_id', 
                            header=None).sequence_length_template.to_dict()
