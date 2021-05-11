import os
from hit_tables import agg_hit_table

hit_table = str(snakemake.input)
fmt = str(snakemake.params.format).upper()

with open(str(snakemake.output), 'wt') as out_handle:

    # identify sequences that hit longer sequences at > 95%
    subseqs = set()
    if os.path.getsize(hit_table) > 0:
        aggs = agg_hit_table(hit_table, format=fmt)
        count = 0
        values = aggs.query('hit != query')[['query', 'hit', 'matches', 'qlen', 'hlen']].values
        for query, hit, matches, qlen, hlen in values:
            if hlen < qlen:
                query, hit, qlen, hlen = hit, query, hlen, qlen
            mfracq = matches / qlen
            if mfracq > .95:
                if query.split('_')[1] != hit.split('_')[1]:
                    subseqs.add(query)
                    count += 1



