"""Display labels for the metrics in evaluate.py's performance TSV.

Each metric maps to (group, name, unit), in report order. "Peaks" are the
fold's test-chrom peaks; "all loci" adds its test-chrom negatives.
"""

METRICS = {
    "count_pearson": ("Counts in peaks", "Pearson", "r of log counts"),
    "count_spearman": ("Counts in peaks", "Spearman", "ρ of log counts"),
    "count_mse": (
        "Counts in peaks",
        "MSE",
        "squared error of ln counts (lower is better)",
    ),
    "all_count_pearson": ("Counts, all loci", "Pearson", "r of log counts"),
    "all_count_spearman": ("Counts, all loci", "Spearman", "ρ of log counts"),
    "all_count_mse": (
        "Counts, all loci",
        "MSE",
        "squared error of ln counts (lower is better)",
    ),
    "profile_mnll": (
        "Profile shape",
        "Multinomial NLL",
        "nats per peak (lower is better)",
    ),
    "profile_jsd": ("Profile shape", "Jensen-Shannon", "divergence (lower is better)"),
    "profile_pearson": ("Profile shape", "Pearson", "r per peak"),
    "profile_spearman": ("Profile shape", "Spearman", "ρ per peak"),
    "auroc": ("Peak vs negative", "AUROC", "scored by predicted count"),
    "auprc": ("Peak vs negative", "AUPRC", "scored by predicted count"),
}
