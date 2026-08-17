# Train the model; best checkpoint -> <prefix>.torch.
rule fit:
    input:
        fasta=FASTA,
        peaks=peaks_for,
        negatives=f"{OUTDIR}/{{sample}}/{{sample}}.negatives.bed",
        signal=lambda wc: signal_bw(wc.sample),
    output:
        model=f"{OUTDIR}/{{sample}}/{{sample}}.torch",
    params:
        name=lambda wc: prefix(wc.sample),
        flags=fit_flags(),
    threads: 8
    log:
        "logs/fit/{sample}.log",
    resources:
        mem_mb=48000,
    conda:
        CONDA_ENV
    shell:
        "python workflow/scripts/fit.py -s {input.fasta} -l {input.peaks} "
        "-neg {input.negatives} -sig {input.signal} -o {params.name} "
        "{params.flags} > {log} 2>&1"


# Score the trained model into a per-signal-group performance table.
rule evaluate:
    input:
        fasta=FASTA,
        peaks=peaks_for,
        signal=lambda wc: signal_bw(wc.sample),
        model=f"{OUTDIR}/{{sample}}/{{sample}}.torch",
    output:
        performance=f"{OUTDIR}/{{sample}}/{{sample}}.performance.tsv",
    params:
        flags=eval_flags(),
    threads: 8
    log:
        "logs/evaluate/{sample}.log",
    resources:
        mem_mb=48000,
    conda:
        CONDA_ENV
    shell:
        "python workflow/scripts/evaluate.py -s {input.fasta} -l {input.peaks} "
        "-sig {input.signal} -m {input.model} -o {output.performance} "
        "{params.flags} > {log} 2>&1"
