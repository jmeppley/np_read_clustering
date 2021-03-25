import pandas
from collections import defaultdict
from jme.jupy_tools.utils import HIT_TABLE_COLUMNS, PAF, BLAST_PLUS, iterable_to_stream

import re
paf_rexp = re.compile(r"""^(?P<query>\S+)\t
                          (?P<qlen>\d+)\t
                          (?P<qstart>\d+)\t
                          (?P<qend>\d+)\t
                          (?P<strand>[+-])\t
                          (?P<hit>\S+)\t
                          (?P<hlen>\d+)\t
                          (?P<hstart>\d+)\t
                          (?P<hend>\d+)\t
                          """,
                      re.VERBOSE)
blast_plus_rexp = re.compile(r"""^
                        (?P<query>\S+)\t
                        (?P<hit>\S+)\t
                        (?P<pctid>[0-9.]+)\t
                        (?P<mlen>\d+)\t
                        (?P<mismatches>\d+)\t
                        (?P<gaps>\d+)\t
                        (?P<qstart>\d+)\t
                        (?P<qend>\d+)\t
                        (?P<hstart>\d+)\t
                        (?P<hend>\d+)\t
                        (?P<evalue>\S+)\t
                        (?P<score>\S+)\t
                        (?P<qlen>\d+)\t
                        (?P<hlen>\d+)\t
                        (?P<raw_score>[0-9.])
                        """,
                             flags=re.VERBOSE)

HIT_TABLE_REXPS = {PAF: paf_rexp, BLAST_PLUS: blast_plus_rexp}

def first(S):
    return next(iter(S))

def generate_nonoverlapping_lines(hit_file, format=BLAST_PLUS, buffer=0):
    hit_rexp = HIT_TABLE_REXPS[format]
    
    with open(hit_file) as lines:
        hit_regions_by_pair = defaultdict(list)
        query_regions_by_pair = defaultdict(list)

        for line in lines:
            if line.startswith("#"):
                # skip comment lines
                continue
            
            try:
                hit_data = hit_rexp.search(line).groupdict()
                query, qstart, qend, hit, hstart, hend = \
                    [hit_data[k] for k in ['query', 'qstart', 'qend', 'hit', 'hstart', 'hend']]
            except AttributeError:
                print("ERROR on line: \n" + line)
                raise
            qstart = int(qstart)
            qend = int(qend)
            hstart = int(hstart)
            hend = int(hend)
            
            pair_key = (hit, query)
            hit_regions = hit_regions_by_pair[pair_key]
            query_regions = query_regions_by_pair[pair_key]

            for regions, start, end in [(hit_regions, hstart, hend),
                                        (query_regions, qstart, qend)]:
                if start > end:
                    start, end = end, start

                # if hit is shorter than buffer, then keep it
                if end - start <= buffer:
                    continue

                # adjust for buffer
                buf_start, buf_end = start + buffer, end - buffer

                # check overlap for all ranges
                for already_occ_range in regions:
                    if (buf_start >= already_occ_range[1] 
                        or buf_end <= already_occ_range[0]):                    
                        # does not overlap this range (try next range) 
                        continue
                    else:
                        # overlaps a previous hit...move on
                        break
                else:
                    # there was no overlap, try next pair or exit loop
                    continue
                # there was some sort of collision
                break
            else:
                # the hit was accepted by both hstart,hend and qstart,qend
                hit_regions.append((hstart,hend))
                query_regions.append((qstart, qend))
                yield line

# only read in hits that don't overlap with other hits (assume better hits are first!)
hit_table_format = str(snakemake.params.format)
hit_table_columns = HIT_TABLE_COLUMNS[hit_table_format]

use_cols = use_cols = list(range(len(hit_table_columns)))
non_ovl_hits = pandas.read_csv(iterable_to_stream(generate_nonoverlapping_lines(str(snakemake.input),
                                                                                format=hit_table_format)), 
                               usecols=use_cols, 
                               sep='\t', 
                               names=hit_table_columns)

if 'matches' not in non_ovl_hits.columns:
    non_ovl_hits = non_ovl_hits.eval('matches = pctid * mlen / 100')

# aggregate all hits by hit/query pair
agg_hits = non_ovl_hits \
    .groupby(['query','hit']) \
    .agg({'matches':sum,'mlen':sum,'qlen':first, 'hlen':first}) \
    .eval('pctid=100 * matches / mlen') \
    .eval('mfrac=200 * matches / (hlen + qlen)') \
    .reset_index()

agg_hits.to_csv(str(snakemake.output), sep='\t', index=None)
