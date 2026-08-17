# Train the model; best checkpoint -> <prefix>.torch.
rule fit:
    input:
        fasta=FASTA,
        peaks=peaks_for,
        negatives=f"{OUTDIR}/{{sample}}/{{sample}}.negatives.bed",
        signal=lambda wc: signal_bw(wc.sample),
        blacklist=blacklist_input,
    output:
        model=f"{OUTDIR}/{{sample}}/{{sample}}.torch",
        history=f"{OUTDIR}/{{sample}}/{{sample}}.log",   # per-epoch training log
    params:
        name=lambda wc: prefix(wc.sample),
        flags=fit_flags(),
    log:
        f"{LOGDIR}/fit/{{sample}}.log",
    benchmark:
        f"{BENCHDIR}/fit/{{sample}}.tsv"
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
        history=f"{OUTDIR}/{{sample}}/{{sample}}.log",
        blacklist=blacklist_input,
    output:
        performance=f"{OUTDIR}/{{sample}}/{{sample}}.performance.tsv",
        epochs=report(
            f"{OUTDIR}/{{sample}}/{{sample}}.epochs.svg",
            category="Training curves",
            labels={"sample": "{sample}"},
        ),
    params:
        flags=eval_flags(),
    log:
        f"{LOGDIR}/evaluate/{{sample}}.log",
    benchmark:
        f"{BENCHDIR}/evaluate/{{sample}}.tsv"
    conda:
        CONDA_ENV
    shell:
        "python workflow/scripts/evaluate.py -s {input.fasta} -l {input.peaks} "
        "-sig {input.signal} -m {input.model} -o {output.performance} "
        "{params.flags} > {log} 2>&1 && "
        "python workflow/scripts/plot_epochs.py -i {input.history} "
        "-o {output.epochs} --sample {wildcards.sample} >> {log} 2>&1"
