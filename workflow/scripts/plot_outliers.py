#!/usr/bin/env python
"""Flag under-performing models and tabulate the worst offenders.

Outliers are the lower tail of `metric` by the Tukey rule (below Q1 - k*IQR).
Writes a rank plot with the fence + flagged models highlighted, and a TSV of the
worst models (ascending) with their covariates for triage.
"""

import argparse

# Self-contained page: pandas table + click-to-sort headers + a filter box.
_TABLE_PAGE = """<!DOCTYPE html><html><head><meta charset="utf-8"><style>
body{font-family:Arial,sans-serif;font-size:12px;margin:8px}
h3{margin:4px 0}
input{margin:6px 0;padding:3px 6px}
table{border-collapse:collapse;width:100%}
th,td{padding:3px 8px;border-bottom:1px solid #ddd;text-align:left;white-space:nowrap}
th{position:sticky;top:0;background:#f4f4f4;cursor:pointer;user-select:none}
tr:hover{background:#f9f9f9}
</style></head><body>
<h3>__TITLE__</h3>
<input id="q" placeholder="filter rows..." oninput="filt()">
__TABLE__
<script>
const t=document.getElementById('wo');
t.querySelectorAll('th').forEach((h,i)=>h.onclick=()=>{
  const asc=h.dataset.asc=h.dataset.asc==='1'?'':'1';
  [...t.tBodies[0].rows].sort((a,b)=>{
    const x=a.cells[i].innerText,y=b.cells[i].innerText,nx=parseFloat(x),ny=parseFloat(y);
    const c=(!isNaN(nx)&&!isNaN(ny))?nx-ny:x.localeCompare(y);
    return asc?c:-c;
  }).forEach(r=>t.tBodies[0].appendChild(r));
});
function filt(){const q=document.getElementById('q').value.toLowerCase();
  [...t.tBodies[0].rows].forEach(r=>r.style.display=r.innerText.toLowerCase().includes(q)?'':'none');}
</script></body></html>"""


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--metrics", required=True, help="metrics.tsv.")
    parser.add_argument("-o", "--output", required=True, help="Destination SVG.")
    parser.add_argument("-t", "--table", required=True,
        help="Destination TSV of the worst offenders.")
    parser.add_argument("-T", "--table_html", default=None,
        help="Optional HTML table (sortable/filterable) for an inline report view.")
    parser.add_argument("-m", "--metric", default="count_pearson",
        help="Metric to screen on (default count_pearson).")
    parser.add_argument("-k", "--iqr_mult", type=float, default=1.5,
        help="Tukey fence multiplier (default 1.5).")
    parser.add_argument("--top", type=int, default=30,
        help="Fallback rows when no outliers are flagged (default 30).")
    return parser


def main():
    args = build_parser().parse_args()

    import numpy as np
    import pandas as pd
    from _style import apply_style
    apply_style()
    import matplotlib.pyplot as plt

    df = pd.read_csv(args.metrics, sep="\t").dropna(subset=[args.metric])
    v = df[args.metric]
    q1, q3 = v.quantile(0.25), v.quantile(0.75)
    fence = q1 - args.iqr_mult * (q3 - q1)
    df = df.sort_values(args.metric, ascending=True).reset_index(drop=True)
    is_out = df[args.metric] < fence

    # Table = every flagged outlier; fall back to the bottom --top if none flagged.
    cols = [c for c in ["sample", "dataset", "genome", "n_fragments", "n_peaks",
                        args.metric, "profile_pearson"] if c in df.columns]
    worst = df[is_out] if is_out.any() else df.head(args.top)
    worst[cols].to_csv(args.table, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(6, 3.6))
    rank = np.arange(len(df))
    ax.scatter(rank[~is_out], df[args.metric][~is_out], s=5, color="#3b6ea5",
        alpha=0.5, edgecolor="none", label="pass")
    ax.scatter(rank[is_out], df[args.metric][is_out], s=14, color="#d1495b",
        edgecolor="none", label=f"outlier (< {fence:.2f})")
    ax.axhline(fence, color="#d1495b", lw=1, ls="--")
    ax.set_xlabel(f"models ranked by {args.metric}")
    ax.set_ylabel(args.metric)
    ax.set_title(f"{args.metric}: {int(is_out.sum())} of {len(df)} below Tukey fence")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.savefig(args.output)

    # A real HTML table so the report iframes it inline (report() TSV/CSV are only
    # download links). pandas builds the table; a tiny self-contained script adds
    # click-to-sort + a filter box.
    if args.table_html:
        disp = worst[cols].copy()
        if "n_fragments" in disp:
            disp["n_fragments"] = disp["n_fragments"].map(lambda v: f"{v / 1e6:.1f}M")
        if "n_peaks" in disp:
            disp["n_peaks"] = disp["n_peaks"].astype(int)
        table = disp.to_html(index=False, table_id="wo", border=0,
            float_format=lambda v: f"{v:.3f}")
        html = (_TABLE_PAGE
                .replace("__TITLE__",
                         f"worst offenders (count_pearson &lt; {fence:.2f}, n={len(disp)})")
                .replace("__TABLE__", table))
        with open(args.table_html, "w") as fh:
            fh.write(html)

    print(f"{int(is_out.sum())} outliers (fence={fence:.3f}); wrote {args.table}")


if __name__ == "__main__":
    main()
