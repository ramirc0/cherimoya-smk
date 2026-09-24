# chrom.sizes are pre-generated per genome (config genomes.<g>.chrom_sizes,
# selected per sample by chrom_sizes_of); no chrom_sizes rule is needed.


# Decompress/copy a provided peak file to the canonical path.
rule prep_peaks:
    input:
        peaks=lambda wc: PEAKS_OF[wc.sample],
    output:
        peaks=f"{OUTDIR}/{{sample}}/{{sample}}.peaks.narrowPeak",
    log:
        f"{LOGDIR}/prep_peaks/{{sample}}.log",
    benchmark:
        f"{BENCHDIR}/prep_peaks/{{sample}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        if printf '%s' {input.peaks:q} | grep -q '\.gz$'; then
            gunzip -c {input.peaks:q} > {output.peaks:q}
        else
            cp {input.peaks:q} {output.peaks:q}
        fi
        """


# Call peaks with macs3 when a sample provides no peak file.
rule macs3:
    input:
        signal=lambda wc: SIGNAL_OF[wc.sample],
        control=macs3_control_input,
    output:
        peaks=f"{OUTDIR}/{{sample}}/{{sample}}_peaks.narrowPeak",
    params:
        fmt=macs3_format,
        gsize=lambda wc: gsize_of(wc.sample),
        q=lambda _: config["preprocess"]["callpeaks_q"],
        name=lambda wc: prefix(wc.sample),
        control=lambda wc, input: ["-c", *input.control] if input.control else [],
        max_count=lambda _: ["--max-count", "1"] if config["preprocess"]["fragments"] else [],
    log:
        f"{LOGDIR}/macs3/{{sample}}.log",
    benchmark:
        f"{BENCHDIR}/macs3/{{sample}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        macs3 callpeak \
            -f {params.fmt:q} \
            -g {params.gsize:q} \
            -n {params.name:q} \
            -q {params.q:q} \
            -t {input.signal:q} \
            {params.control:q} \
            {params.max_count:q}
        """


# bam2bw flags shared by the signal and control conversions (strand, fragments,
# read-depth). bam2bw emits `<name>.bw` with -u, else the `<name>.+.bw`/`.-.bw`
# pair, so the same flags drive both the declared outputs and the tool.
def _bam2bw_extra():
    return [
        flag
        for flag, on in [
            ("-u", config["preprocess"]["unstranded"]),
            ("-f", config["preprocess"]["fragments"]),
            ("-r", config["preprocess"]["read_depth"]),
        ]
        if on
    ]


# Convert the signal BAM/fragments to a bigWig (stranded -> (+, -) pair).
rule bam2bw:
    input:
        signal=lambda wc: SIGNAL_OF[wc.sample],
        sizes=lambda wc: chrom_sizes_of(wc.sample),
    output:
        bw=(
            [f"{OUTDIR}/{{sample}}/{{sample}}.+.bw",
             f"{OUTDIR}/{{sample}}/{{sample}}.-.bw"]
            if STRANDED else f"{OUTDIR}/{{sample}}/{{sample}}.bw"
        ),
    params:
        name=lambda wc: prefix(wc.sample),
        pos_shift=lambda _: config["preprocess"]["pos_shift"],
        neg_shift=lambda _: config["preprocess"]["neg_shift"],
        scale_factor=lambda _: config["preprocess"]["scale_factor"],
        extra=lambda _: _bam2bw_extra(),
    log:
        f"{LOGDIR}/bam2bw/{{sample}}.log",
    benchmark:
        f"{BENCHDIR}/bam2bw/{{sample}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bam2bw \
            -s {input.sizes:q} \
            -n {params.name:q} \
            -ps {params.pos_shift:q} \
            -ns {params.neg_shift:q} \
            -sf {params.scale_factor:q} \
            -p {threads} \
            {params.extra:q} \
            {input.signal:q}
        """


# Convert the control BAM/fragments to a bigWig model track (stranded -> pair).
# Only built for samples with a control; no scale factor, matching cherimoya.
rule bam2bw_control:
    input:
        control=lambda wc: CONTROL_OF[wc.sample],
        sizes=lambda wc: chrom_sizes_of(wc.sample),
    output:
        bw=(
            [f"{OUTDIR}/{{sample}}/{{sample}}.control.+.bw",
             f"{OUTDIR}/{{sample}}/{{sample}}.control.-.bw"]
            if STRANDED else f"{OUTDIR}/{{sample}}/{{sample}}.control.bw"
        ),
    params:
        name=lambda wc: f"{prefix(wc.sample)}.control",
        pos_shift=lambda _: config["preprocess"]["pos_shift"],
        neg_shift=lambda _: config["preprocess"]["neg_shift"],
        extra=lambda _: _bam2bw_extra(),
    log:
        f"{LOGDIR}/bam2bw_control/{{sample}}.log",
    benchmark:
        f"{BENCHDIR}/bam2bw_control/{{sample}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bam2bw \
            -s {input.sizes:q} \
            -n {params.name:q} \
            -ps {params.pos_shift:q} \
            -ns {params.neg_shift:q} \
            -p {threads} \
            {params.extra:q} \
            {input.control:q}
        """


# Fragment count of the signal, the depth QC covariate (full scan).
rule count_fragments:
    input:
        signal=lambda wc: SIGNAL_OF[wc.sample],
    output:
        n_fragments=f"{OUTDIR}/{{sample}}/{{sample}}.n_fragments.txt",
    params:
        fragments=lambda _: ["--fragments"] if config["preprocess"]["fragments"] else [],
        paired_end=lambda _: ["--paired_end"] if config["preprocess"]["paired_end"] else [],
    log:
        f"{LOGDIR}/count_fragments/{{sample}}.log",
    benchmark:
        f"{BENCHDIR}/count_fragments/{{sample}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/count_fragments.py \
            -i {input.signal:q} \
            -o {output.n_fragments:q} \
            {params.fragments:q} \
            {params.paired_end:q}
        """
