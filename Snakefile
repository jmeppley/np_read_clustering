from Bio import SeqIO

nanopore_fasta = config.get('fasta_file', 'nanopore_reads.fasta')
output_dir = config.get('output_dir', 'spike_search')
window_size = config.get('window_size', 2_000)
overlap = config.get('window_overlap', int(window_size/2))
window_spacing = window_size - overlap
min_size = config.get('minimum_size', 0)
mcl_i = "{:0.1f}".format(config.get('mcl_i', 5.0))

def get_cluster_mcl_files(wildcards):
    return [f'{output_dir}/windows/cluster/window.{window}/self.I{mcl_i}.mcl' 
            for window in get_windows()]

rule all:
    input: get_cluster_mcl_files

rule mcl:
    input: '{file_root}.abc'
    output: '{file_root}.I{mcl_i}.mcl'
    shell:
        "mcl {input} --abc -I {wildcards.mcl_i} -o {output} > {output}.log 2>&1"

rule window_abc:
    input: '{output_dir}/windows/cluster/window.{window}/self.paf'
    output: '{output_dir}/windows/cluster/window.{window}/self.abc'
    run:
        from jme.jupy_tools.utils import parse_blast_m8, PAF
        paf_file = str(input)
        if os.path.getsize(paf_file) > 0:
            mm_hits = parse_blast_m8(paf_file, format=PAF)
            mm_hits \
                .query('mlen > .4 * (qlen + hlen)') \
                .eval('mfrac = 2 * matches / (hlen + qlen)') \
                .set_index(['query', 'hit']) \
                .mfrac \
                .to_csv(str(output), sep='\t', header=None)
        else:
            # input is empty, touch the output
            with open(str(output), 'wt') as out_handle:
                pass

rule window_minimap:
    input: '{output_dir}/windows/reads/reads.{window}.fasta'
    output: '{output_dir}/windows/cluster/window.{window}/self.paf'
    threads: 10
    shell:
        """
        minimap2 -x ava-ont -t 20 {input} {input} > {output} 2> {output}.log
        """

def get_windows():
    reads_dir = checkpoints.apply_fasta_size_window.get().output.fasta
    windows, = glob_wildcards(reads_dir + "/reads.{window}.fasta")
    return windows

checkpoint apply_fasta_size_window:
    input: nanopore_fasta
    output: 
        fasta=directory(f'{output_dir}/windows/reads'),
        counts=f'{output_dir}/windows/read_counts.txt'
    run:
        from jme.jupy_tools.utils import LogLogger
        llogger_read = LogLogger(multiplier=10)
        llogger_window = LogLogger(multiplier=2)
        handles = {}
        counts = {}
        os.makedirs(str(output.fasta), exist_ok=True)
        for i, read in enumerate(SeqIO.parse(str(input), 'fasta')):
            l = len(read)
            if l < min_size:
                continue
            
            windows = []
            
            window_num = int(l / (window_spacing))
            window_start = window_num * window_spacing
            window_end = window_start + window_size
            
            # shift down incase we'd be in the end of an earlier window
            while window_end > l + window_spacing:
                window_num -= 1
                window_start = window_num * window_spacing
                window_end = window_start + window_size

            # loop over all windows this read is in
            while window_start <= l:
                window_end = window_start + window_size
                window = f"{window_start}_{window_end}"

                counts[window] = counts.get(window, 0) + 1

                windows.append(window)

                try:
                    fasta_out = handles[window]
                except:
                    out_file = f"{output.fasta}/reads.{window}.fasta"
                    llogger_window.log("opnning file %d (%s) for window %s", (len(handles), out_file, window))
                    fasta_out = open(out_file, 'wt')
                    handles[window] = fasta_out

                fasta_out.write(read.format('fasta'))
            
                window_num += 1
                window_start = window_num * window_spacing
            
            llogger_read.log("wrote read %d (%s) of len %d to window %s", (i, read.id, l, windows))
            
        
        print(f'Closing {len(handles)} handles')
        for handle in handles.values():
            handle.close()
            
        with open(str(output.counts), 'wt') as counts_out:
            for window, count in counts.items():
                counts_out.write(f"{window}\t{count}\n")
            
            
