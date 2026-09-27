# Hypothetical attributions of one CV fold's model over its peaks.
rule attribute:
    input:
        fasta=lambda wc: fasta_of(wc.sample),
        peaks=peaks_for,
        model=f"{OUTDIR}/{{sample}}/fold_{{fold}}/{{sample}}.torch",
        fold=lambda wc: fold_json(wc.sample, wc.fold),
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
