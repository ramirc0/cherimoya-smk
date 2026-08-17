# chrom sizes from the .fai (bam2bw -s accepts a sizes file).
rule chrom_sizes:
    input:
        fai=FAI,
    output:
        sizes=CHROM_SIZES,
    log:
        "logs/chrom_sizes/refs.log",
    conda:
        CONDA_ENV
    shell:
        "cut -f1,2 {input.fai} > {output.sizes} 2> {log}"


# Decompress/copy a provided peak file to the canonical path.
rule prep_peaks:
    input:
        peaks=lambda wc: PEAKS_OF[wc.sample],
    output:
        peaks=f"{OUTDIR}/{{sample}}/{{sample}}.peaks.narrowPeak",
    log:
        "logs/prep_peaks/{sample}.log",
    conda:
        CONDA_ENV
    shell:
        "if printf '%s' '{input.peaks}' | grep -q '\\.gz$'; then "
        "gunzip -c {input.peaks} > {output.peaks}; "
        "else cp {input.peaks} {output.peaks}; fi 2> {log}"


# Call peaks with macs3 when a sample provides no peak file.
rule macs3:
    input:
        signal=lambda wc: SIGNAL_OF[wc.sample],
        control=macs3_control_input,
    output:
        peaks=f"{OUTDIR}/{{sample}}/{{sample}}_peaks.narrowPeak",
    params:
        fmt=macs3_format,
        gsize=config["preprocess"]["callpeaks_gsize"],
        q=config["preprocess"]["callpeaks_q"],
        name=lambda wc: prefix(wc.sample),
        control=lambda wc, input: f"-c {input.control}" if input.control else "",
        max_count="--max-count 1" if config["preprocess"]["fragments"] else "",
    log:
        "logs/macs3/{sample}.log",
    resources:
        mem_mb=8000,
    conda:
        CONDA_ENV
    shell:
        "macs3 callpeak -f {params.fmt} -g {params.gsize} -n {params.name} "
        "-q {params.q} -t {input.signal} {params.control} {params.max_count} "
        "> {log} 2>&1"


# Convert the signal BAM/fragments to a bigWig.
rule bam2bw:
    input:
        signal=lambda wc: SIGNAL_OF[wc.sample],
        sizes=CHROM_SIZES,
    output:
        bw=f"{OUTDIR}/{{sample}}/{{sample}}.bw",
    params:
        name=lambda wc: prefix(wc.sample),
        pos_shift=config["preprocess"]["pos_shift"],
        neg_shift=config["preprocess"]["neg_shift"],
        scale_factor=config["preprocess"]["scale_factor"],
        extra=(
            ("-u " if config["preprocess"]["unstranded"] else "")
            + ("-f " if config["preprocess"]["fragments"] else "")
            + ("-r " if config["preprocess"]["read_depth"] else "")
        ).strip(),
    threads: 4
    log:
        "logs/bam2bw/{sample}.log",
    resources:
        mem_mb=8000,
    conda:
        CONDA_ENV
    shell:
        "bam2bw -s {input.sizes} -n {params.name} -ps {params.pos_shift} "
        "-ns {params.neg_shift} -sf {params.scale_factor} -p {threads} "
        "{params.extra} {input.signal} > {log} 2>&1"
