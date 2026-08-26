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
        flags=lambda _: fit_flags(),
    log:
        f"{LOGDIR}/fit/{{sample}}.log",
    benchmark:
        f"{BENCHDIR}/fit/{{sample}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/fit.py \
            -s {input.fasta:q} \
            -l {input.peaks:q} \
            -neg {input.negatives:q} \
            -sig {input.signal:q} \
            -o {params.name:q} \
            {params.flags:q}
        """


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
        flags=lambda _: eval_flags(),
    log:
        f"{LOGDIR}/evaluate/{{sample}}.log",
    benchmark:
        f"{BENCHDIR}/evaluate/{{sample}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/evaluate.py \
            -s {input.fasta:q} \
            -l {input.peaks:q} \
            -sig {input.signal:q} \
            -m {input.model:q} \
            -o {output.performance:q} \
            {params.flags:q}

        python workflow/scripts/plot_epochs.py \
            -i {input.history:q} \
            -o {output.epochs:q} \
            --sample {wildcards.sample:q}
        """
