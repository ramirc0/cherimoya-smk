# Hypothetical attributions of one CV fold's model over its peaks.
rule attribute:
    input:
        fasta=lambda wc: fasta_of(wc.sample),
        peaks=peaks_for,
        model=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.torch",
        fold=lambda wc: fold_json(wc.sample, wc.fold),
        blacklist=blacklist_input,
    output:
        ohe=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.attributions.ohe.npz",
        attr=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.attributions.attr.npz",
        idxs=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.attributions.idxs.npy",
    params:
        flags=lambda wc: attr_flags(wc.sample, wc.fold),
    log:
        f"{LOGDIR}/attribute/{{sample}}.fold_{{fold}}.log",
    benchmark:
        f"{BENCHDIR}/attribute/{{sample}}.fold_{{fold}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/attribute.py \
            {params.flags:q} \
            -s {input.fasta:q} \
            -l {input.peaks:q} \
            -m {input.model:q} \
            --ohe_filename {output.ohe:q} \
            --attr_filename {output.attr:q} \
            --idx_filename {output.idxs:q}
        """


# Seqlets called from one fold's attributions, in genome coordinates.
rule seqlets:
    input:
        peaks=peaks_for,
        ohe=rules.attribute.output.ohe,
        attr=rules.attribute.output.attr,
        idxs=rules.attribute.output.idxs,
    output:
        bed=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.seqlets.bed",
    params:
        flags=lambda wc: seqlet_flags(wc.sample, wc.fold),
    log:
        f"{LOGDIR}/seqlets/{{sample}}.fold_{{fold}}.log",
    benchmark:
        f"{BENCHDIR}/seqlets/{{sample}}.fold_{{fold}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/seqlets.py \
            {params.flags:q} \
            -l {input.peaks:q} \
            --ohe_filename {input.ohe:q} \
            --attr_filename {input.attr:q} \
            --idx_filename {input.idxs:q} \
            --output_filename {output.bed:q}
        """


# Nearest motifs per seqlet (TomTom) and the seqlet count per best motif.
rule annotate:
    input:
        fasta=lambda wc: fasta_of(wc.sample),
        bed=rules.seqlets.output.bed,
        motifs=config["annotate"]["motifs"],
    output:
        bed=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.seqlets_annotated.bed",
        counts=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.motif_seqlet_count.tsv",
    params:
        flags=annot_flags(),
    log:
        f"{LOGDIR}/annotate/{{sample}}.fold_{{fold}}.log",
    benchmark:
        f"{BENCHDIR}/annotate/{{sample}}.fold_{{fold}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/annotate.py \
            {params.flags:q} \
            -s {input.fasta:q} \
            -b {input.bed:q} \
            -t {input.motifs:q} \
            --n_jobs {threads} \
            --output_filename {output.bed:q} \
            --count_filename {output.counts:q}
        """


# TF-MoDISco patterns from one fold's attributions.
rule modisco_motifs:
    input:
        ohe=rules.attribute.output.ohe,
        attr=rules.attribute.output.attr,
    output:
        h5=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.modisco_results.h5",
    params:
        flags=lambda wc: modisco_flags(),
    log:
        f"{LOGDIR}/modisco_motifs/{{sample}}.fold_{{fold}}.log",
    benchmark:
        f"{BENCHDIR}/modisco_motifs/{{sample}}.fold_{{fold}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/modisco_motifs.py \
            {params.flags:q} \
            -s {input.ohe:q} \
            -a {input.attr:q} \
            -o {output.h5:q} \
            --verbose
        """


# HTML report of the patterns with their TomTom matches in annotate.motifs.
# Own env: the MEME-suite tomtom binary cannot share cherimoya.yaml's pins.
rule modisco_report:
    input:
        h5=rules.modisco_motifs.output.h5,
        motifs=config["annotate"]["motifs"],
    output:
        report=directory(f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.modisco"),
    params:
        flags=lambda wc: modisco_report_flags(),
    log:
        f"{LOGDIR}/modisco_report/{{sample}}.fold_{{fold}}.log",
    benchmark:
        f"{BENCHDIR}/modisco_report/{{sample}}.fold_{{fold}}.tsv"
    conda:
        "../envs/modisco_report.yaml"
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/modisco_report.py report \
            {params.flags:q} \
            -i {input.h5:q} \
            -o {output.report:q} \
            -s ./ \
            -m {input.motifs:q}
        """


# HTML report of each annotate.motifs consensus's effect when inserted into peaks.
rule marginalize:
    input:
        fasta=lambda wc: fasta_of(wc.sample),
        peaks=peaks_for,
        model=rules.attribute.input.model,
        fold=lambda wc: fold_json(wc.sample, wc.fold),
        motifs=config["annotate"]["motifs"],
    output:
        report=directory(f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.marginalize"),
    params:
        flags=lambda wc: marginalize_flags(wc.sample, wc.fold),
    log:
        f"{LOGDIR}/marginalize/{{sample}}.fold_{{fold}}.log",
    benchmark:
        f"{BENCHDIR}/marginalize/{{sample}}.fold_{{fold}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/marginalize.py \
            {params.flags:q} \
            -s {input.fasta:q} \
            -l {input.peaks:q} \
            -m {input.model:q} \
            -t {input.motifs:q} \
            -o {output.report:q} \
            --verbose
        """


# Mean attribution by position across the attributed peaks.
rule attribution_profile:
    input:
        ohe=rules.attribute.output.ohe,
        attr=rules.attribute.output.attr,
    output:
        plot=report(
            f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.attribution_profile.svg",
            category="Interpretation",
            labels={"sample": "{sample}", "fold": "{fold}", "plot": "attribution profile"},
        ),
        plot_png=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.attribution_profile.png",
    log:
        f"{LOGDIR}/attribution_profile/{{sample}}.fold_{{fold}}.log",
    benchmark:
        f"{BENCHDIR}/attribution_profile/{{sample}}.fold_{{fold}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_attribution_profile.py \
            --ohe_filename {input.ohe:q} \
            --attr_filename {input.attr:q} \
            -o {output.plot:q} \
            --sample {wildcards.sample:q}
        """


rule seqlet_lengths:
    input:
        bed=rules.seqlets.output.bed,
    output:
        plot=report(
            f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.seqlet_lengths.svg",
            category="Interpretation",
            labels={"sample": "{sample}", "fold": "{fold}", "plot": "seqlet lengths"},
        ),
        plot_png=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.seqlet_lengths.png",
    log:
        f"{LOGDIR}/seqlet_lengths/{{sample}}.fold_{{fold}}.log",
    benchmark:
        f"{BENCHDIR}/seqlet_lengths/{{sample}}.fold_{{fold}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_seqlet_lengths.py \
            -i {input.bed:q} \
            -o {output.plot:q} \
            --sample {wildcards.sample:q}
        """


rule motif_counts:
    input:
        counts=rules.annotate.output.counts,
    output:
        plot=report(
            f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.motif_counts.svg",
            category="Interpretation",
            labels={"sample": "{sample}", "fold": "{fold}", "plot": "top motifs"},
        ),
        plot_png=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.motif_counts.png",
    log:
        f"{LOGDIR}/motif_counts/{{sample}}.fold_{{fold}}.log",
    benchmark:
        f"{BENCHDIR}/motif_counts/{{sample}}.fold_{{fold}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_motif_counts.py \
            -i {input.counts:q} \
            -o {output.plot:q} \
            --sample {wildcards.sample:q}
        """
