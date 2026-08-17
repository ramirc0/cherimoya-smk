# GC-matched negative regions for a sample's peaks.
rule negatives:
    input:
        peaks=peaks_for,
        fasta=FASTA,
    output:
        bed=f"{OUTDIR}/{{sample}}/{{sample}}.negatives.bed",
    params:
        bin_width=config["negatives"]["bin_width"],
        max_n_perc=config["negatives"]["max_n_perc"],
        beta=config["negatives"]["beta"],
        in_window=config["negatives"]["in_window"],
        out_window=config["negatives"]["out_window"],
    log:
        f"{LOGDIR}/negatives/{{sample}}.log",
    benchmark:
        f"{BENCHDIR}/negatives/{{sample}}.tsv"
    conda:
        CONDA_ENV
    shell:
        "python workflow/scripts/negatives.py -i {input.peaks} -f {input.fasta} "
        "-o {output.bed} -l {params.bin_width} -n {params.max_n_perc} "
        "-a {params.beta} -w {params.in_window} -x {params.out_window} > {log} 2>&1"
