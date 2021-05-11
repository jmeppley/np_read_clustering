import os
from hit_tables import agg_hit_table

hit_table = str(snakemake.input)
fmt = str(snakemake.params.format).upper()

if os.path.getsize(hit_table) > 0:
    agg_hit_table(hit_table, format=fmt) \
        .to_csv(str(snakemake.output), sep='\t', index=None)
else:
    # input is empty, touch the output
    with open(str(snakemake.output), 'wt') as out_handle:
        pass

