# Train one CV fold; best checkpoint -> <sample>/fold_<fold>/<sample>.torch.
# bigWig/peaks/negatives are fold-agnostic (shared); only the train/valid chroms differ.
rule fit:
    input:
        fasta=lambda wc: fasta_of(wc.sample),
        peaks=peaks_for,
        negatives=f"{OUTDIR}/{{sample}}/{{sample}}.negatives.bed",
        signal=lambda wc: signal_bw(wc.sample),
        control=lambda wc: control_bw(wc.sample),
        fold=lambda wc: fold_json(wc.sample, wc.fold),
        blacklist=blacklist_input,
    output:
        model=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.torch",
        history=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.log",  # per-epoch log
    params:
        name=lambda wc: fold_prefix(wc.sample, wc.fold),
        flags=lambda wc: fit_flags(wc.sample, wc.fold),
        control_flag=lambda wc: ["-c"] if control_bw(wc.sample) else [],
    log:
        f"{LOGDIR}/fit/{{sample}}.fold_{{fold}}.log",
    benchmark:
        f"{BENCHDIR}/fit/{{sample}}.fold_{{fold}}.tsv"
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
            {params.control_flag:q} {input.control:q} \
            -o {params.name:q} \
            {params.flags:q}
        """


# Score the trained model into a per-signal-group performance table.
rule evaluate:
    input:
        fasta=lambda wc: fasta_of(wc.sample),
        peaks=peaks_for,
        signal=lambda wc: signal_bw(wc.sample),
        control=lambda wc: control_bw(wc.sample),
        model=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.torch",
        history=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.log",
        fold=lambda wc: fold_json(wc.sample, wc.fold),
        blacklist=blacklist_input,
    output:
        performance=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.performance.tsv",
        counts=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.counts.tsv",
        epochs=report(
            f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.epochs.svg",
            category="Training curves",
            labels={"sample": "{sample}", "fold": "{fold}"},
        ),
        count_scatter=report(
            f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.count_scatter.svg",
            category="Count QC",
            labels={"sample": "{sample}", "fold": "{fold}"},
        ),
        epochs_png=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.epochs.png",
        count_scatter_png=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.count_scatter.png",
    params:
        flags=lambda wc: eval_flags(wc.sample, wc.fold),
        control_flag=lambda wc: ["-c"] if control_bw(wc.sample) else [],
    log:
        f"{LOGDIR}/evaluate/{{sample}}.fold_{{fold}}.log",
    benchmark:
        f"{BENCHDIR}/evaluate/{{sample}}.fold_{{fold}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/evaluate.py \
            -s {input.fasta:q} \
            -l {input.peaks:q} \
            -sig {input.signal:q} \
            {params.control_flag:q} {input.control:q} \
            -m {input.model:q} \
            -o {output.performance:q} \
            --counts_filename {output.counts:q} \
            {params.flags:q}

        python workflow/scripts/plot_epochs.py \
            -i {input.history:q} \
            -o {output.epochs:q} \
            --sample {wildcards.sample:q}

        python workflow/scripts/plot_count_scatter.py \
            -i {output.counts:q} \
            -o {output.count_scatter:q} \
            --sample {wildcards.sample:q}
        """
