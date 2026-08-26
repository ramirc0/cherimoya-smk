# chrom sizes from the .fai (bam2bw -s accepts a sizes file).
rule chrom_sizes:
    input:
        fai=FAI,
    output:
        sizes=CHROM_SIZES,
    log:
        f"{LOGDIR}/chrom_sizes/refs.log",
    benchmark:
        f"{BENCHDIR}/chrom_sizes/refs.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        cut -f1,2 {input.fai:q} > {output.sizes:q}
        """


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
        gsize=lambda _: config["preprocess"]["callpeaks_gsize"],
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


# Convert the signal BAM/fragments to a bigWig.
rule bam2bw:
    input:
        signal=lambda wc: SIGNAL_OF[wc.sample],
        sizes=CHROM_SIZES,
    output:
        bw=f"{OUTDIR}/{{sample}}/{{sample}}.bw",
    params:
        name=lambda wc: prefix(wc.sample),
        pos_shift=lambda _: config["preprocess"]["pos_shift"],
        neg_shift=lambda _: config["preprocess"]["neg_shift"],
        scale_factor=lambda _: config["preprocess"]["scale_factor"],
        extra=lambda _: [
            flag
            for flag, on in [
                ("-u", config["preprocess"]["unstranded"]),
                ("-f", config["preprocess"]["fragments"]),
                ("-r", config["preprocess"]["read_depth"]),
            ]
            if on
        ],
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
