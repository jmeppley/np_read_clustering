import pandas
from jme.jupy_tools.hit_tables import parse_blast_m8, PAF
from Bio import SeqIO

# pick a best read
hits = parse_blast_m8(str(snakemake.input.paf),format=PAF)
hit_matches = hits.groupby(['hit','query']).agg({'matches':sum})
mean_matches = {r:hit_matches.query(f'hit != query and (hit == "{r}" or query == "{r}")').matches.mean() 
                for r in set(i[0] for i in hit_matches.index).union(i[1] for i in hit_matches.index)}
best_matches = sorted(mean_matches.keys(), key=lambda r: mean_matches[r], reverse=True)
ref_read = best_matches[0]

# write out to 2 files
with open(str(snakemake.output.ref), 'wt') as ref_out:
    with open(str(snakemake.output.others), 'wt') as others_out:
        for read in SeqIO.parse(str(snakemake.input.fasta), 'fasta'):
            if read.id == ref_read:
                ref_out.write(read.format('fasta'))
            else:
                others_out.write(read.format('fasta'))
