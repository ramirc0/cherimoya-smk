# GC-matched negative regions for a sample's peaks.
rule negatives:
    input:
        peaks=peaks_for,
        fasta=lambda wc: fasta_of(wc.sample),
    output:
        bed=f"{OUTDIR}/{{sample}}/{{sample}}.negatives.bed",
    params:
        bin_width=lambda _: config["negatives"]["bin_width"],
        max_n_perc=lambda _: config["negatives"]["max_n_perc"],
        beta=lambda _: config["negatives"]["beta"],
        in_window=lambda _: config["negatives"]["in_window"],
        out_window=lambda _: config["negatives"]["out_window"],
        on_missing_contig=lambda _: config["peaks"]["on_missing_contig"],
    log:
        f"{LOGDIR}/negatives/{{sample}}.log",
    benchmark:
        f"{BENCHDIR}/negatives/{{sample}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/negatives.py \
            -i {input.peaks:q} \
            -f {input.fasta:q} \
            -o {output.bed:q} \
            -l {params.bin_width:q} \
            -n {params.max_n_perc:q} \
            -a {params.beta:q} \
            -w {params.in_window:q} \
            -x {params.out_window:q} \
            --on_missing_contig {params.on_missing_contig:q}
        """
