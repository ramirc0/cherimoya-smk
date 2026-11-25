# GC-matched negative regions for a sample's peaks.
rule negatives:
    input:
        peaks=peaks_for,
        fasta=lambda wc: fasta_of(wc.sample),
    output:
        bed=f"{OUTDIR}/{{sample}}/{{sample}}.negatives.bed",
    params:
        bin_width=lambda wc: config_of(wc.sample)["negatives"]["bin_width"],
        max_n_perc=lambda wc: config_of(wc.sample)["negatives"]["max_n_perc"],
        beta=lambda wc: config_of(wc.sample)["negatives"]["beta"],
        in_window=lambda wc: config_of(wc.sample)["negatives"]["in_window"],
        out_window=lambda wc: config_of(wc.sample)["negatives"]["out_window"],
        random_state=lambda wc: config_of(wc.sample)["random_state"],
        on_missing_contig=lambda wc: config_of(wc.sample)["peaks"]["on_missing_contig"],
    log:
        f"{LOGDIR}/negatives/{{sample}}.log",
    benchmark:
        f"{BENCHDIR}/negatives/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/negatives.py \
            -i {input.peaks:q} \
            -f {input.fasta:q} \
            -o {output.bed:q} \
            -l {params.bin_width:q} \
            -n {params.max_n_perc:q} \
            -a {params.beta:q} \
            -w {params.in_window:q} \
            -x {params.out_window:q} \
            -r {params.random_state:q} \
            --on_missing_contig {params.on_missing_contig:q}

        sync {output:q}
        """
