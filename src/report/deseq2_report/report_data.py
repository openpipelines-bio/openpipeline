"""Report data for report/deseq2_report.

Reads the output of differential_expression/deseq2 --export_normalized_counts
and computes everything the report shows (summary metrics, PCA, correlation,
heatmap, volcano/MA points, key genes, p-value histogram, results table and
methods text) as one JSON-serializable dict. report.qmd only plots it.
"""

import datetime
import html
import json
import os
import platform

import numpy as np
import pandas as pd


def esc(value):
    return html.escape(str(value), quote=False)


def sig3(value):
    """Round to 3 significant figures, keeping tiny p-values."""
    return float(f"{value:.3g}")


def r_squared(values, factor):
    """Share of the variance of `values` explained by the levels of `factor`."""
    df = pd.DataFrame({"v": values, "f": factor})
    total = ((df["v"] - df["v"].mean()) ** 2).sum()
    between = df.groupby("f")["v"].transform("mean").sub(df["v"].mean()).pow(2).sum()
    return float(between / total) if total > 0 else 0.0


def read_deseq2_output(input_dir, prefix):
    """Read the files written by differential_expression/deseq2 --export_normalized_counts."""
    base = os.path.join(input_dir, prefix)
    paths = {
        "results": f"{base}.csv",
        "samples": f"{base}_samples.csv",
        "normalized": f"{base}_normalized_counts.csv",
        "vst": f"{base}_vst.csv",
        "metadata": f"{base}_metadata.json",
    }
    missing = [p for p in paths.values() if not os.path.isfile(p)]
    if missing:
        raise FileNotFoundError(
            f"Missing DESeq2 output files: {missing}. Run differential_expression/deseq2 "
            "with --export_normalized_counts and check --input_prefix."
        )
    with open(paths["metadata"]) as f:
        metadata = json.load(f)
    samples = pd.read_csv(paths["samples"], dtype=str).set_index("sample")
    normalized = pd.read_csv(paths["normalized"])
    vst = pd.read_csv(paths["vst"])
    results = pd.read_csv(paths["results"])
    return results, samples, normalized, vst, metadata


def build_report_data(par, logger, tool="OpenPipeline report/deseq2_report"):
    """Collect everything the report shows from the DESeq2 output, as a JSON-serializable dict."""
    rng = np.random.default_rng(par["seed"])
    results, samples, normalized, vst, deseq2_meta = read_deseq2_output(
        par["input"], par["input_prefix"]
    )
    padj_threshold = (
        par["p_adj_threshold"]
        if par["p_adj_threshold"] is not None
        else deseq2_meta["p_adj_threshold"]
    )
    lfc_threshold = (
        par["log2fc_threshold"]
        if par["log2fc_threshold"] is not None
        else deseq2_meta["log2fc_threshold"]
    )
    contrast_column = deseq2_meta["contrast_column"]
    if len(deseq2_meta["contrasts"]) != 1:
        raise ValueError(
            "The report supports a single contrast, found "
            f"{[c['name'] for c in deseq2_meta['contrasts']]}"
        )
    contrast = deseq2_meta["contrasts"][0]
    test_level, ref_level = contrast["comparison_group"], contrast["control_group"]

    for col in [contrast_column, par["obs_pair"], *(par["obs_sample_label"] or [])]:
        if col and col not in samples.columns:
            raise ValueError(
                f"Column '{col}' not found in the DESeq2 sample table, which holds the design "
                "formula, contrast and cell group columns: "
                f"{[c for c in samples.columns if c not in ('size_factor', 'library_size')]}"
            )

    # ---- genes and counts ----
    gene_ids = pd.Index(normalized.pop("gene_id").astype(str))
    names = normalized.pop("gene_symbol") if "gene_symbol" in normalized else None
    vst.index = vst.pop("gene_id").astype(str)
    vst = vst.drop(columns="gene_symbol", errors="ignore")
    symbols = pd.Series(gene_ids, index=gene_ids)
    if names is not None:
        names = pd.Series(names.astype(str).to_numpy(), index=gene_ids)
        symbols = names.where(names.notna() & (names != "") & (names != "nan"), symbols)

    # ---- sample order and labels: reference group first, then test group ----
    group = samples[contrast_column].astype(str)
    label_map = dict(g.split("=", 1) for g in (par["group_labels"] or []))

    def disp(level):
        return label_map.get(level, level)

    if par["obs_sample_label"]:
        sample = samples[par["obs_sample_label"]].astype(str).agg(" ".join, axis=1)
    else:
        sample = pd.Series(samples.index, index=samples.index)
    rank = group.map({ref_level: 0, test_level: 1}).fillna(2)
    order = np.lexsort((sample.to_numpy(), rank.to_numpy()))
    samples, group, sample = samples.iloc[order], group.iloc[order], sample.iloc[order]
    pair = samples[par["obs_pair"]].astype(str).to_numpy() if par["obs_pair"] else None
    pair_label = par["pair_label"] or par["obs_pair"]
    group_disp = group.map(disp).tolist()
    sample_list = sample.tolist()

    size_factor = samples["size_factor"].astype(float).to_numpy()
    library_size = samples["library_size"].astype(float).to_numpy()
    norm = normalized[samples.index].to_numpy(dtype=float).T  # samples x genes
    # DESeq2 normalized counts are the raw counts divided by the size factor
    counts = np.round(norm * size_factor[:, None])

    # ---- DE results ----
    results["gene_id"] = results["gene_id"].astype(str)
    results["gene"] = results["gene_id"].map(symbols).fillna(results["gene_id"])
    tested = results.dropna(subset=["pvalue", "padj", "log2FoldChange"]).copy()
    # Same rule as the "significant" column of differential_expression/deseq2
    is_sig = (tested["padj"] < padj_threshold) & (
        tested["log2FoldChange"].abs() > lfc_threshold
    )
    tested["cls"] = np.where(
        ~is_sig, "ns", np.where(tested["log2FoldChange"] > 0, "up", "down")
    )
    sig = tested[tested["cls"] != "ns"]
    n_up = int((tested["cls"] == "up").sum())
    n_down = int((tested["cls"] == "down").sum())
    logger.info(
        "%s vs %s: %d genes tested, %d up, %d down",
        test_level,
        ref_level,
        len(tested),
        n_up,
        n_down,
    )

    # ---- sample-level panels on the variance-stabilized counts ----
    keep = counts.sum(axis=0) >= par["min_count"]
    kept_ids = gene_ids[keep]
    vst_kept = (
        vst.reindex(index=kept_ids, columns=samples.index).to_numpy(dtype=float).T
    )
    if not np.all(np.isfinite(vst_kept)):
        raise ValueError("The VST table does not cover all samples and genes.")

    top_var = np.argsort(vst_kept.var(axis=0, ddof=1))[::-1][: par["n_pca_genes"]]
    x = vst_kept[:, top_var] - vst_kept[:, top_var].mean(axis=0)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    pcs = u * s
    var_pct = s**2 / np.sum(s**2) * 100 if np.sum(s**2) > 0 else np.zeros_like(s)
    if pcs.shape[1] < 2:
        pcs = np.column_stack([pcs, np.zeros((pcs.shape[0], 2 - pcs.shape[1]))])
        var_pct = np.append(var_pct, np.zeros(2 - len(var_pct)))
    # Orient each PC so the test group sits on the positive side
    is_test = (group == test_level).to_numpy()
    for k in range(2):
        if pcs[is_test, k].mean() < pcs[~is_test, k].mean():
            pcs[:, k] *= -1

    # Name what each PC follows, if one factor clearly explains it
    pc_notes = []
    for k in range(2):
        fits = [("group", r_squared(pcs[:, k], group.to_numpy()))]
        if pair is not None:
            fits.append(("pair", r_squared(pcs[:, k], pair)))
        what, r2 = max(fits, key=lambda f: f[1])
        if r2 >= 0.5:
            desc = (
                f"separates {esc(disp(test_level))} from {esc(disp(ref_level))}"
                if what == "group"
                else f"follows {esc(pair_label)}"
            )
            pc_notes.append(f"PC{k + 1} {desc} (R&sup2; {r2:.2f})")
    pca_hint = "Each point is a sample." + (
        " " + "; ".join(pc_notes) + "." if pc_notes else ""
    )

    corr = pd.DataFrame(vst_kept.T).corr(method="spearman").to_numpy()

    # Heatmap: lowest-padj significant genes, up first then down, z-scored per gene
    top = sig[sig["gene_id"].isin(kept_ids)].nsmallest(par["n_heatmap_genes"], "padj")
    top = pd.concat([top[top["cls"] == "up"], top[top["cls"] == "down"]])
    col_of = pd.Series(np.arange(len(kept_ids)), index=kept_ids)
    hm = vst_kept[:, col_of[top["gene_id"]].to_numpy()].T
    sd = hm.std(axis=1, keepdims=True)
    hm_z = (hm - hm.mean(axis=1, keepdims=True)) / np.where(sd > 0, sd, 1)

    # ---- volcano / MA: all significant genes + a random subset of the rest ----
    ns_all = tested[tested["cls"] == "ns"]
    ns_pts = ns_all
    if len(ns_pts) > par["max_volcano_ns_genes"]:
        ns_pts = ns_pts.iloc[
            rng.choice(len(ns_pts), par["max_volcano_ns_genes"], replace=False)
        ]
    vol = pd.concat([ns_pts, sig])
    neglog_p = -np.log10(vol["pvalue"].clip(lower=1e-300))

    # ---- p-value histogram; pi0 (share of true nulls), Storey's estimator at lambda = 0.5 ----
    hist, _ = np.histogram(tested["pvalue"], bins=par["pvalue_bins"], range=(0, 1))
    pi0 = float(min(1.0, (tested["pvalue"] > 0.5).mean() / 0.5)) if len(tested) else 1.0
    n_nonnull = int(round((1 - pi0) * len(tested)))

    # ---- highlights ----
    best = tested.nsmallest(1, "padj").iloc[0] if len(tested) else None
    robust = sig[sig["baseMean"] >= par["min_basemean_highlight"]]
    strongest = (
        robust.loc[robust["log2FoldChange"].abs().idxmax()] if len(robust) else None
    )

    # ---- key genes: normalized counts per sample ----
    if par["highlight_genes"]:
        picks = []
        for gene in par["highlight_genes"]:
            hits = tested[(tested["gene"] == gene) | (tested["gene_id"] == gene)]
            if hits.empty:
                logger.warning("Highlight gene %s is not among the tested genes", gene)
            else:
                picks.append(hits.nsmallest(1, "padj").iloc[0])
        key = pd.DataFrame(picks)
    else:
        key = robust.reindex(
            robust["log2FoldChange"].abs().sort_values(ascending=False).index
        ).head(par["n_highlight_genes"])
    gene_pos = pd.Series(np.arange(len(gene_ids)), index=gene_ids)
    genes_panel = {
        "samples": sample_list,
        "groups": group_disp,
        "pairs": pair.tolist() if pair is not None else None,
        "pair_label": pair_label,
        "items": [
            {
                "gene": r["gene"],
                "gene_id": r["gene_id"],
                "log2fc": round(float(r["log2FoldChange"]), 2),
                "padj": sig3(r["padj"]),
                "values": np.round(norm[:, gene_pos[r["gene_id"]]] + 0.5, 1).tolist(),
            }
            for _, r in key.iterrows()
        ],
    }

    table = [
        {
            "gene": r.gene,
            "gene_id": r.gene_id,
            "baseMean": round(float(r.baseMean), 1),
            "log2fc": round(float(r.log2FoldChange), 3),
            "pvalue": sig3(r.pvalue),
            "padj": sig3(r.padj),
        }
        for r in tested.nsmallest(par["n_table_genes"], "padj").itertuples()
    ]

    # ---- KPIs and tiles ----
    group_sizes = group.value_counts()
    min_reps = int(min(group_sizes.get(ref_level, 0), group_sizes.get(test_level, 0)))
    kpis = {
        "samples": int(len(samples)),
        "genes_quantified": int((counts.sum(axis=0) > 0).sum()),
        "median_mapping": None,
        "de_genes": n_up + n_down,
        "up": n_up,
        "down": n_down,
        "padj_threshold": padj_threshold,
        "log2fc_threshold": lfc_threshold,
        "shrunk": bool(deseq2_meta["lfc_shrinkage"]),
    }
    de_label = (
        f"DE genes (padj<{padj_threshold:g}"
        + (f", |log2FC|&gt;{lfc_threshold:g}" if lfc_threshold else "")
        + ")"
    )
    tiles = [
        {"l": "Samples", "v": str(kpis["samples"])},
        {"l": "Genes detected", "v": f"{kpis['genes_quantified']:,}"},
        {"l": "Median library size", "v": f"{np.median(library_size) / 1e6:.1f} M"},
        {"l": "Genes tested", "v": f"{len(tested):,}"},
        {"l": "Min. replicates", "v": str(min_reps)},
        {
            "l": de_label,
            "v": f"{kpis['de_genes']:,}",
            "accent": True,
            "sub": f'<span class="pos">&#9650; {n_up:,}</span> &nbsp; '
            f'<span class="neg">&#9660; {n_down:,}</span>',
        },
    ]

    # ---- methods, from the DESeq2 run metadata and the report parameters ----
    sig_rule = f"padj &lt; {padj_threshold:g}" + (
        f" and |log2FC| &gt; {lfc_threshold:g}" if lfc_threshold else ""
    )
    vs = deseq2_meta["variance_stabilization"]
    cell_group = deseq2_meta.get("cell_group")
    design = esc(par["design_description"] or "")
    methods = [
        [
            "Design",
            (design + "; " if design else "")
            + f"model <code>{esc(deseq2_meta['design_formula'])}</code>"
            + (
                f"; samples with {esc(cell_group['column'])} = {esc(cell_group['value'])}"
                if cell_group
                else ""
            ),
        ],
        [
            "Counts",
            esc(par["counts_description"] or f"raw counts from {deseq2_meta['input']}")
            + f"; {kpis['genes_quantified']:,} genes detected, {int(keep.sum()):,} with "
            f"&ge; {par['min_count']} total counts used for the sample-level panels",
        ],
        [
            "Differential expression",
            f"DESeq2 {esc(deseq2_meta['versions']['DESeq2'])}, "
            f"{esc(deseq2_meta['test'])} test {esc(test_level)} vs {esc(ref_level)}, "
            f"{esc(deseq2_meta['p_adjust_method'])} adjustment; significant at {sig_rule}; "
            f"fold changes {'shrunken' if kpis['shrunk'] else 'not shrunken'}",
        ],
        [
            "Expression scale",
            f"DESeq2 variance-stabilizing transformation ({esc(vs['function_name'])}), "
            f"{'blind to' if vs['blind'] else 'using'} the design",
        ],
        [
            "Ordination",
            f"PCA on the {min(par['n_pca_genes'], int(keep.sum())):,} most variable genes",
        ],
        [
            "Correlation",
            f"Spearman between samples, all {int(keep.sum()):,} retained genes",
        ],
        [
            "Heatmap",
            f"{len(top)} significant genes with the lowest padj (up, then down), "
            "z-scored per gene",
        ],
        [
            "Volcano &amp; MA",
            f"all {len(sig):,} significant genes plus {len(ns_pts):,} of {len(ns_all):,} "
            "non-significant genes, sampled at random",
        ],
        [
            "Key genes",
            "DESeq2 normalized counts + 0.5, log scale"
            + (
                f"; lines join samples from the same {esc(pair_label)}"
                if pair is not None
                else ""
            ),
        ],
        [
            "P-values",
            f"raw p-values of the {len(tested):,} tested genes; "
            f"&pi;<sub>0</sub> = {pi0:.2f} (Storey, &lambda; = 0.5)",
        ],
    ]

    data = {
        "meta": {
            "title": par["title"],
            "project": par["project"] or "",
            "design": par["design_description"] or "",
            "pipeline": par["pipeline"] or "",
            "ref": par["reference"] or "",
            "attribution": par["attribution"] or "",
        },
        "kpis": kpis,
        "tiles": tiles,
        "summary": {
            "top": None
            if best is None
            else {
                "gene": best["gene"],
                "log2fc": round(float(best["log2FoldChange"]), 2),
                "padj": sig3(best["padj"]),
            },
            "strongest": None
            if strongest is None
            else {
                "gene": strongest["gene"],
                "log2fc": round(float(strongest["log2FoldChange"]), 2),
                "baseMean": round(float(strongest["baseMean"]), 1),
                "min_basemean": par["min_basemean_highlight"],
            },
            "pi0": round(pi0, 3),
            "n_nonnull": n_nonnull,
        },
        "pca": {
            "pc1": np.round(pcs[:, 0], 2).tolist(),
            "pc2": np.round(pcs[:, 1], 2).tolist(),
            "sample": sample_list,
            "group": group_disp,
            "var": np.round(var_pct[:2], 1).tolist(),
            "hint": pca_hint,
        },
        "volcano": {
            "x": vol["log2FoldChange"].round(3).tolist(),
            "y": neglog_p.round(3).tolist(),
            "base_mean": vol["baseMean"].round(2).tolist(),
            "gene": vol["gene"].tolist(),
            "cls": vol["cls"].tolist(),
        },
        "pvalues": {
            "counts": hist.tolist(),
            "bin_width": 1 / par["pvalue_bins"],
            "n": int(len(tested)),
            "pi0": round(pi0, 3),
            "n_nonnull": n_nonnull,
        },
        "corr": {"z": np.round(corr, 3).tolist(), "labels": sample_list},
        "heatmap": {
            "z": np.round(hm_z, 3).tolist(),
            "genes": top["gene"].tolist(),
            "samples": sample_list,
        },
        "genes": genes_panel,
        "table": table,
        "methods": methods,
        "provenance": {
            "tool": tool,
            "generated": datetime.datetime.now().isoformat(timespec="seconds"),
            "contrast": {
                "column": contrast_column,
                "test": test_level,
                "reference": ref_level,
            },
            "deseq2": deseq2_meta,
            "versions": {
                **deseq2_meta["versions"],
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "python": platform.python_version(),
            },
        },
    }

    return data
