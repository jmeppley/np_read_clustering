import logging
import io
import re
import numpy
import pandas
from contextlib import contextmanager
from collections import defaultdict

LOGGER = logging.getLogger(__name__)

BLAST_COLUMNS = [
        "query",
        "hit",
        "pctid",
        "mlen",
        "mismatches",
        "gaps",
        "qstart",
        "qend",
        "hstart",
        "hend",
        "evalue",
        "score",
]
BLAST_PLUS_COLUMNS = BLAST_COLUMNS + [
    'qlen', 'hlen', 'raw_score'
]
PAF_COLUMNS = ['query','qlen','qstart','qend','strand',
               'hit', 'hlen','hstart','hend','matches',
               'mlen','mapqv']

BLAST = 'BlastTab'.upper()
BLAST_PLUS = 'BlastTab+'.upper()
PAF = 'Paf'.upper()
PAF_ALL = 'Paf+'.upper()
LASTAL = 'lastal'.upper()
LASTAL_COLUMNS = ['score',
                  'hit', 'hstart', 'hmlen', 'hstrand', 'hlen',
                  'query', 'qstart', 'qmlen', 'qstrand', 'qlen',
                  'match_string', 'eg2', 'e']


HIT_TABLE_COLUMNS = {BLAST: BLAST_COLUMNS,
                     BLAST_PLUS: BLAST_PLUS_COLUMNS,
                     PAF: PAF_COLUMNS,
                     PAF_ALL: PAF_COLUMNS + ['tp','cm','dv','rl'],
                     LASTAL: LASTAL_COLUMNS,
                    }

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

def computeLastHitValues(blocks):
    """
    given the query length and 'blocks' string from a last hit
    return the:
        match length

    the blocks string looks something like this:
        "73,0:1,15,0:1,13,0:1,9"
        where integer elements indicate lenghts of matches and
        colon separated elements indicate lengths of gaps
    """
    matchLen = 0
    for segment in blocks.split(','):
        try:
            matchLen += int(segment)
        except ValueError:
            (hml, qml) = segment.split(":")
            mml = max(int(hml), int(qml))
            matchLen += mml

    return matchLen


def parse_blast_m8(hit_table, format=BLAST, skiprows=0, pair_overlap_buffer=-1, **cutoffs):
    """ utility for quickly loading a hit table into pandas 
    
        cutoff keys should match column names. Use negative cutoff value to accept all values.
        
        if pair_overlap_buffer is set to a nonzero value, hits that overlap previous hits 
        (in the hit_table order) between a given query/hit pair are dropped. This happens
        before any other curofs are applied.
    """

    if pair_overlap_buffer >= 0:
        # call this function again after filtering for overlaps
        if format not in HIT_TABLE_REXPS:
            raise Exception("Overlap filtering is not supported for format: " + format)
        with remove_pair_overlaps(hit_table, buffer=pair_overlap_buffer, format=format) as new_table:
            return parse_blast_m8(new_table, format=format, skiprows=skiprows, **cutoffs)
    
    column_names = HIT_TABLE_COLUMNS[format]
    use_cols = list(range(len(column_names)))
    hits = \
        pandas.read_csv(
            hit_table,
            sep="\t",
            header=None,
            comment="#",
            usecols=use_cols,
            names=column_names,
            skiprows=skiprows,
        )
    
    # format specific tweaks
    if format == LASTAL:
        hits['qsmult'] = [1 if qs == '+' else -1 
                          for qs in hits.qstrand]
        hits = \
            hits.eval('hend = hstart + hmlen - 1') \
                .eval('qend = qstart + ((qmlen - 1) * qsmult)')
        hits['mlen'] = [computeLastHitValues(blocks)
                        for blocks in hits.match_string]
        hits['evalue'] = [float(re.sub('E=','',str(e)))
                          for e in hits['e']]
    if format in [PAF, PAF_ALL]:
        # calculate pctid
        hits = hits.eval('pctid = 100 * matches / mlen')
    
    query = " and ".join(f"{column} >= {value}" 
                         if value >=0 
                         else f"{column} <= {-1*value}"
                         for column, value in cutoffs.items()
                         if numpy.abs(value) > 0
                        )
    if len(query) > 0:
        return hits.query(query)
    return hits

@contextmanager
def remove_pair_overlaps(hit_table, **params):
    """ wrapper to  generate_nonoverlapping_lines that (1) returns a file-like
    object instead of an iterator and (2) serves as a context manager """
    yield iterable_to_stream(generate_nonoverlapping_lines(hit_table, **params))


def generate_nonoverlapping_lines(hit_file, format=BLAST_PLUS, buffer=0):
    """
    Loops over lines in a hit_table keeping track of hit and query positions. Currently
    supported formats are paf and blasttab+

    Only yields lines that don't overlap with previous hits between a query/hit pair
    """
    
    if isinstance(hit_file, str):
        # open file handle and call this function again
        with open(hit_file) as lines:
            yield from generate_nonoverlapping_lines(lines, format=format, buffer=buffer)
    else:
        hit_rexp = HIT_TABLE_REXPS[format]
        hit_regions_by_pair = defaultdict(list)
        query_regions_by_pair = defaultdict(list)

        for line in hit_file:
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
            
            LOGGER.debug('processing line: \n%s', line)

            for regions, start, end in [(hit_regions, hstart, hend),
                                        (query_regions, qstart, qend)]:
                if start > end:
                    start, end = end, start
                    
                LOGGER.debug('checking for %d - %d for hit in \n%r', start, end, regions)

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
                        LOGGER.debug("overlaps %r", already_occ_range)
                        break
                else:
                    # there was no overlap, try next pair or exit loop
                    LOGGER.debug("no overlap here")
                    continue
                # there was some sort of collision
                break
            else:
                # the hit was accepted by both hstart,hend and qstart,qend
                LOGGER.debug("no overlaps")
                hit_regions.append(tuple(sorted((hstart,hend))))
                query_regions.append(tuple(sorted((qstart, qend))))
                yield line

def first(S):
    """ Simply return the first item from a collection (using next(iter(S))) """
    return next(iter(S))

def agg_hit_df(non_ovl_hits):
    # aggregate all hits by hit/query pair
    if 'matches' not in non_ovl_hits.columns:
        non_ovl_hits = non_ovl_hits.eval('matches = pctid * mlen / 100')

    agg_hits = non_ovl_hits \
        .groupby(['query','hit']) \
        .agg({'matches':sum,'mlen':sum,'qlen':first, 'hlen':first}) \
        .eval('pctid=100 * matches / mlen') \
        .eval('mfrac=200 * matches / (hlen + qlen)') \
        .reset_index()

    return agg_hits    
                
def agg_hit_table(hit_table, ovl_buffer=0, **parse_args):
    """
    for each hit/query pair return one line of data
    """
    non_ovl_hits = parse_blast_m8(hit_table, pair_overlap_buffer=ovl_buffer, **parse_args)

    return agg_hit_df(non_ovl_hits)


def iterable_to_stream(iterable, str_to_bytes=True, buffer_size=io.DEFAULT_BUFFER_SIZE):
    """
    Lets you use an iterable (e.g. a generator) that yields bytestrings as a read-only
    input stream.

    The stream implements Python 3's newer I/O API (available in Python 2's io module).
    For efficiency, the stream is buffered.
    
    src: https://stackoverflow.com/a/20260030/663466
    """
    
    if str_to_bytes:
        # encode strings as bytes
        iterable = (s.encode('utf-8') for s in iterable)

    class IterStream(io.RawIOBase):
        def __init__(self):
            self.leftover = None
        def readable(self):
            return True
        def readinto(self, b):
            try:
                l = len(b)  # We're supposed to return at most this much
                chunk = self.leftover or next(iterable)
                output, self.leftover = chunk[:l], chunk[l:]
                b[:len(output)] = output
                return len(output)
            except StopIteration:
                return 0    # indicate EOF
    return io.BufferedReader(IterStream(), buffer_size=buffer_size)
