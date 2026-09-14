# Run-level QC plots. All read one tidy metrics table (gather_metrics) that joins
# every model's performance with its genome, depth, and peak count.

_ALL_PERFORMANCE = [f"{OUTDIR}/{s}/fold_{fold}/{s}.performance.tsv"
                    for s in SAMPLES for fold in FOLDS]


rule gather_metrics:
    input:
        performance=_ALL_PERFORMANCE,
        samples=config["samples"],
    output:
        metrics=f"{OUTDIR}/report/metrics.tsv",
    params:
        results_dir=OUTDIR,
    log:
        f"{LOGDIR}/gather_metrics/all.log",
    benchmark:
        f"{BENCHDIR}/gather_metrics/all.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/gather_metrics.py \
            -i {input.performance:q} \
            -s {input.samples:q} \
            -r {params.results_dir:q} \
            -o {output.metrics:q}
        """


# count_pearson vs a per-model covariate (n_peaks | n_fragments).
rule perf_vs_covariate:
    input:
        metrics=f"{OUTDIR}/report/metrics.tsv",
    output:
        plot=report(
            f"{OUTDIR}/report/count_pearson_vs_{{covariate}}.svg",
            category="Count QC",
            labels={"plot": "count_pearson vs {covariate}"},
        ),
    wildcard_constraints:
        covariate="|".join(COVARIATES),
    log:
        f"{LOGDIR}/perf_vs_covariate/{{covariate}}.log",
    benchmark:
        f"{BENCHDIR}/perf_vs_covariate/{{covariate}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_perf_vs_covariate.py \
            -i {input.metrics:q} \
            -o {output.plot:q} \
            -x {wildcards.covariate:q} \
            -y count_pearson
        """


# Under-performing models: rank plot + per-model metrics table with an outlier flag.
rule outliers:
    input:
        metrics=f"{OUTDIR}/report/metrics.tsv",
    output:
        plot=report(f"{OUTDIR}/report/outliers.svg",
                    category="Count QC", labels={"plot": "outliers"}),
        table=report(f"{OUTDIR}/report/outliers.tsv",
                     category="Count QC", labels={"table": "per-model metrics + outlier flag"}),
    log:
        f"{LOGDIR}/outliers/all.log",
    benchmark:
        f"{BENCHDIR}/outliers/all.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_outliers.py \
            -i {input.metrics:q} \
            -o {output.plot:q} \
            -t {output.table:q}
        """


# Run-level plot: distribution of model performance, faceted across metrics.


rule performance_distribution:
    input:
        performance=_ALL_PERFORMANCE,
    output:
        plot=report(
            f"{OUTDIR}/report/performance_distribution.svg",
            category="Performance",
            labels={"plot": "distribution"},
        ),
    log:
        f"{LOGDIR}/performance_distribution/all.log",
    benchmark:
        f"{BENCHDIR}/performance_distribution/all.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_performance_dist.py \
            -i {input.performance:q} \
            -o {output.plot:q}
        """
