# Run-level plot: distribution of model performance, faceted across metrics.


rule performance_distribution:
    input:
        performance=expand(
            f"{OUTDIR}/{{sample}}/{{sample}}.performance.tsv", sample=SAMPLES
        ),
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
