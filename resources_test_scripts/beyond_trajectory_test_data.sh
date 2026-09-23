#!/bin/bash
# Builds the BEYOND test fixture from a real, openly licensed snRNA-seq cohort.
#
# Source
# ------
#   PsychAD RADC_Cohort, dorsolateral prefrontal cortex, 152 donors,
#   693682 nuclei x 34176 genes, 10x 3' v3.
#
#   CZ CELLxGENE Discover collection:
#     https://cellxgene.cziscience.com/collections/84ce6837-548d-4a1f-919f-0bc0d9a3952f
#   Publication: doi:10.1101/2024.10.31.24316513
#   Licence: CC-BY 4.0. Attribution is required; redistribution is not restricted.
#   The submitters certify the data as non-identifiable, so no individual-level
#   controlled-access data is involved.
#
# Why real data
# -------------
# The previous fixture was a hand-written negative-binomial simulation. Composition,
# trajectory and trait effects were all planted by the same script that the workflow was
# then asked to recover, so any structure the real method depends on but the simulation
# did not reproduce was invisible. This cohort carries the structure BEYOND is built for:
# a three-level annotation hierarchy (class 8 / subclass 27 / subtype 65) matching
# BEYOND's cell class / cell type / subpopulation, and donor-level neuropathology.
#
# The BEYOND cellular landscape is a landscape of *participants*, not of cells, so the
# proportion matrix is the unit of analysis. Subsetting cells damages it, which is why
# `proportions.csv` below is computed from all 693682 nuclei while `atlas.h5mu` carries a
# subsample. The two are consistent by construction: the same within-class normalisation
# on the same annotation.
#
# Produces: resources_test/beyond_test_data/
#   atlas.h5mu        - 32377 nuclei x ~300 genes; all 152 donors (167-226 nuclei each,
#                       stratified by subtype); obs["participant_id"], obs["cell_class"],
#                       obs["subpopulation"]. Feeds stats/calculate_label_proportions.
#   proportions.csv   - 152 donors x 65 subtypes, within-class prevalence, computed from
#                       ALL 693682 nuclei. Feeds every group-level step.
#   traits.csv        - 152 donors x 8 traits. AD_status / Parkinson_disease /
#                       Vascular_status vary; Schizophrenia and ASCVD_status are
#                       single-level in this cohort and are therefore true null traits.
#   de_EN.csv         - AD vs non-AD differential expression per cell class, computed
#   de_IN.csv           from donor pseudobulk of the nuclei above (Welch t-test on
#   de_Astro.csv        log2 CPM, BH correction).
#   gene_sets.gmt     - one set of class-specific marker genes per cell class, ranked by
#                       expression specificity in this dataset. Independent of the AD
#                       contrast above, so the enrichment step is not circular.
#
# The 6.28 GB source file is downloaded once and kept. Set BEYOND_RADC_H5AD to reuse a
# copy that is already on disk.

set -eo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

ID=beyond_test_data
OUT="resources_test/$ID"
mkdir -p "$OUT"

RADC_URL="https://datasets.cellxgene.cziscience.com/54293783-669c-410e-919d-474960f8761b.h5ad"
RADC_H5AD="${BEYOND_RADC_H5AD:-$OUT/temp_RADC_Cohort.h5ad}"

if [ ! -f "$RADC_H5AD" ]; then
  echo "Downloading PsychAD RADC_Cohort (6.28 GB) to $RADC_H5AD"
  mkdir -p "$(dirname "$RADC_H5AD")"
  curl -L --fail --retry 3 -o "$RADC_H5AD" "$RADC_URL"
else
  echo "Reusing $RADC_H5AD"
fi

python3 <<PYCODE
import os
import numpy as np
import pandas as pd
import h5py
import scipy.sparse as sp
from anndata import AnnData
from mudata import MuData

SRC = "$RADC_H5AD"
OUT = "$OUT"

SEED               = 0
CELLS_PER_DONOR    = 200   # cap; donors with fewer nuclei keep all of them
N_MARKERS_PER_CLASS = 25   # size of each gene set, and of the class part of the panel
N_TOP_EXPRESSED     = 100  # genes added to the panel regardless of class specificity
DE_CLASSES          = ["EN", "IN", "Astro"]

rng = np.random.default_rng(SEED)
f = h5py.File(SRC, "r")

# -- obs -----------------------------------------------------------------------
def _decode(values):
    return np.array(
        [v.decode() if isinstance(v, bytes) else v for v in values], dtype=object
    )


def categorical(group, key):
    """Read an AnnData column, which may be stored plainly or as a categorical."""
    g = group[key]
    if isinstance(g, h5py.Group) and "categories" in g:
        return pd.Categorical.from_codes(g["codes"][:], _decode(g["categories"][:]))
    return _decode(g[:])

OBS_KEYS = [
    "donor_id", "class", "subclass", "subtype", "sex", "genetic_ancestry",
    "development_stage", "AD_status", "Parkinson_disease", "Vascular_status",
    "DLBD_status", "Schizophrenia", "ASCVD_status",
]
obs = pd.DataFrame({k: categorical(f["obs"], k) for k in OBS_KEYS})
n_cells_total = len(obs)
print(f"Source: {n_cells_total} nuclei, {obs.donor_id.nunique()} donors, "
      f"{obs['class'].nunique()} classes, {obs.subtype.nunique()} subtypes")

# -- proportions: within-class prevalence over ALL nuclei -----------------------
#
# The reference implementation normalises within the cell-type grouping, not over all
# subpopulations of a donor (2. Cell-type analysis/3.create.proportion.matrix.R:32-37),
# so each donor's row sums to the number of classes they have nuclei in, not to 1.
counts = (
    obs.groupby(["class", "subtype", "donor_id"], observed=True)
    .size()
    .rename("n")
    .reset_index()
)
counts["prevalence"] = counts["n"] / counts.groupby(
    ["class", "donor_id"], observed=True
)["n"].transform("sum")
proportions = counts.pivot_table(
    index="donor_id", columns="subtype", values="prevalence",
    fill_value=0, aggfunc="sum", observed=True,
)
proportions.index.name = "participant_id"
proportions.to_csv(f"{OUT}/proportions.csv")
print(f"Wrote {OUT}/proportions.csv  "
      f"({proportions.shape[0]} donors x {proportions.shape[1]} subtypes, "
      f"from all {n_cells_total} nuclei)")

# -- donor traits ---------------------------------------------------------------
traits = obs.groupby("donor_id", observed=True).agg(lambda s: s.iloc[0])
traits = traits[[
    "AD_status", "Parkinson_disease", "Vascular_status", "DLBD_status",
    "Schizophrenia", "ASCVD_status", "sex", "genetic_ancestry",
]].copy()
# "68-year-old stage" -> 68, so the association tests have a continuous covariate
age = obs.groupby("donor_id", observed=True)["development_stage"].first().astype(str)
traits["age"] = age.str.extract(r"(\d+)").astype(float).to_numpy()
traits.index.name = "participant_id"
traits.to_csv(f"{OUT}/traits.csv")
print(f"Wrote {OUT}/traits.csv  ({traits.shape[0]} donors x {traits.shape[1]} traits)")

# -- stratified cell subset ------------------------------------------------------
keep = []
for donor, group in obs.groupby("donor_id", observed=True):
    if len(group) <= CELLS_PER_DONOR:
        keep.append(group.index.to_numpy())
        continue
    frac = CELLS_PER_DONOR / len(group)
    for _, sub in group.groupby("subtype", observed=True):
        k = max(1, int(round(len(sub) * frac)))
        keep.append(rng.choice(sub.index.to_numpy(), size=min(k, len(sub)), replace=False))
selected = np.sort(np.concatenate(keep))
sub_obs = obs.loc[selected]
print(f"Subset: {len(selected)} nuclei, {sub_obs.donor_id.nunique()} donors, "
      f"{sub_obs.subtype.nunique()} subtypes")

# -- read the selected rows of the CSR matrix ------------------------------------
X = f["X"]
indptr = X["indptr"][:]
data_ds, indices_ds = X["data"], X["indices"]
n_genes_total = int(X.attrs["shape"][1])

rows, cols, vals = [], [], []
for i, cell in enumerate(selected):
    start, stop = int(indptr[cell]), int(indptr[cell + 1])
    if stop <= start:
        continue
    idx = indices_ds[start:stop]
    cols.append(idx)
    vals.append(data_ds[start:stop])
    rows.append(np.full(len(idx), i, dtype=np.int32))
counts_matrix = sp.csr_matrix(
    (np.concatenate(vals),
     (np.concatenate(rows), np.concatenate(cols).astype(np.int32))),
    shape=(len(selected), n_genes_total),
    dtype=np.float32,
)
print(f"Read {counts_matrix.nnz} non-zero entries")

gene_names = np.asarray(categorical(f["var"], "feature_name"), dtype=object)

# -- log CPM, used for markers and for DE ----------------------------------------
library = np.asarray(counts_matrix.sum(axis=1)).ravel()
library[library == 0] = 1.0
cpm = counts_matrix.multiply(1e6 / library[:, None]).tocsr()
logcpm = cpm.copy()
logcpm.data = np.log2(logcpm.data + 1.0)

cell_class = sub_obs["class"].to_numpy()
classes = list(pd.unique(sub_obs["class"]))

# -- class marker genes ----------------------------------------------------------
#
# Specificity = mean log2 CPM inside the class minus mean log2 CPM outside it. Computed
# from expression only, so the sets are independent of the AD contrast used for the DE
# tables and the enrichment step is not testing the same numbers twice.
class_means = {}
for cls in classes:
    mask = cell_class == cls
    class_means[cls] = np.asarray(logcpm[mask].mean(axis=0)).ravel()
overall_mean = np.asarray(logcpm.mean(axis=0)).ravel()
n_per_class = {cls: int((cell_class == cls).sum()) for cls in classes}

markers = {}
used = set()
for cls in classes:
    rest = (overall_mean * len(selected) - class_means[cls] * n_per_class[cls]) / max(
        1, len(selected) - n_per_class[cls]
    )
    specificity = class_means[cls] - rest
    # A marker has to be expressed, not merely relatively enriched in a sparse gene
    expressed = class_means[cls] > 0.5
    order = np.argsort(-np.where(expressed, specificity, -np.inf))
    picked = []
    for j in order:
        if not np.isfinite(specificity[j]) or not expressed[j]:
            break
        name = gene_names[j]
        if name in used:
            continue
        picked.append(j)
        used.add(name)
        if len(picked) == N_MARKERS_PER_CLASS:
            break
    markers[cls] = picked
    print(f"  {cls:7s} {len(picked)} markers, e.g. {list(gene_names[picked[:4]])}")

with open(f"{OUT}/gene_sets.gmt", "w") as handle:
    for cls, idx in markers.items():
        genes = [str(gene_names[j]) for j in idx]
        handle.write("\t".join([f"{cls}_MARKERS", f"{cls} class markers"] + genes) + "\n")
print(f"Wrote {OUT}/gene_sets.gmt  ({len(markers)} sets)")

# -- gene panel for the atlas ------------------------------------------------------
# Every marker gene, plus the most-expressed genes so the panel is not only markers.
# Duplicate gene symbols would collide in var_names, so keep the first index per symbol.
panel, seen = [], set()

def add_gene(j):
    name = str(gene_names[j])
    if name in seen:
        return False
    seen.add(name)
    panel.append(int(j))
    return True

for idx in markers.values():
    for j in idx:
        add_gene(j)
n_marker_genes = len(panel)

added = 0
for j in np.argsort(-overall_mean):
    if added == N_TOP_EXPRESSED:
        break
    if add_gene(int(j)):
        added += 1

panel = np.array(sorted(panel))
print(f"Gene panel: {len(panel)} genes ({n_marker_genes} markers + {added} most expressed)")

# -- differential expression per class, over ALL nuclei ------------------------------
#
# Donor pseudobulk needs every nucleus, not the 200-per-donor subsample: with the
# subsample no gene survives BH correction, which leaves the enrichment step nothing to
# rank. X is CSR, so it is streamed in row blocks and summed into a
# (class, donor) x gene count matrix; the full matrix is never held in memory.
from scipy.stats import ttest_ind

ad_status = traits["AD_status"]
ad_donors = set(ad_status.index[ad_status == "Yes"])
ctrl_donors = set(ad_status.index[ad_status == "No"])

all_donors = list(proportions.index)
donor_index = {d: i for i, d in enumerate(all_donors)}
obs_donor_idx = obs["donor_id"].astype(str).map(donor_index).to_numpy()
obs_class = obs["class"].astype(str).to_numpy()

pseudobulk_counts = {
    cls: np.zeros((len(all_donors), n_genes_total), dtype=np.float64)
    for cls in DE_CLASSES
}
BLOCK = 20000
for a in range(0, n_cells_total, BLOCK):
    b = min(a + BLOCK, n_cells_total)
    start_ptr, stop_ptr = int(indptr[a]), int(indptr[b])
    block = sp.csr_matrix(
        (
            data_ds[start_ptr:stop_ptr],
            indices_ds[start_ptr:stop_ptr].astype(np.int32),
            indptr[a : b + 1] - start_ptr,
        ),
        shape=(b - a, n_genes_total),
    )
    for cls in DE_CLASSES:
        mask = obs_class[a:b] == cls
        if not mask.any():
            continue
        rows_here = np.flatnonzero(mask)
        donors_here = obs_donor_idx[a:b][rows_here]
        selector = sp.csr_matrix(
            (
                np.ones(len(rows_here)),
                (donors_here, rows_here),
            ),
            shape=(len(all_donors), b - a),
        )
        pseudobulk_counts[cls] += np.asarray((selector @ block).todense())
    print(f"  pseudobulk {b}/{n_cells_total} nuclei", flush=True)

f.close()


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    if ok.sum() == 0:
        return q
    ranked = np.argsort(p[ok])
    m = int(ok.sum())
    adj = p[ok][ranked] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[ranked] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


for cls in DE_CLASSES:
    pb = pseudobulk_counts[cls]
    library = pb.sum(axis=1)
    present = library > 0
    cpm = np.zeros_like(pb)
    cpm[present] = pb[present] / library[present, None] * 1e6
    logcpm_pb = np.log2(cpm + 1.0)

    is_ad = np.array([d in ad_donors and present[i] for i, d in enumerate(all_donors)])
    is_ctrl = np.array([d in ctrl_donors and present[i] for i, d in enumerate(all_donors)])
    stat, pval = ttest_ind(
        logcpm_pb[is_ad], logcpm_pb[is_ctrl], axis=0, equal_var=False
    )
    lfc = logcpm_pb[is_ad].mean(axis=0) - logcpm_pb[is_ctrl].mean(axis=0)
    base = cpm[present].mean(axis=0)
    expressed = base > 1.0

    de = pd.DataFrame(
        {
            "gene": gene_names[expressed],
            "baseMean": base[expressed],
            "log2FoldChange": lfc[expressed],
            "stat": stat[expressed],
            "pvalue": pval[expressed],
            "padj": bh(pval[expressed]),
        }
    ).set_index("gene")
    de = de[~de.index.duplicated()].sort_values("pvalue")
    de.to_csv(f"{OUT}/de_{cls}.csv")
    n_sig = int((de["padj"] < 0.05).sum())
    print(f"Wrote {OUT}/de_{cls}.csv  ({len(de)} genes, {n_sig} with padj < 0.05, "
          f"{int(is_ad.sum())} AD vs {int(is_ctrl.sum())} control donors, "
          f"all {n_cells_total} nuclei)")

# -- atlas.h5mu ---------------------------------------------------------------------
atlas_obs = pd.DataFrame(
    {
        "participant_id": sub_obs["donor_id"].astype(str).to_numpy(),
        "cell_class": sub_obs["class"].astype(str).to_numpy(),
        "subclass": sub_obs["subclass"].astype(str).to_numpy(),
        "subpopulation": sub_obs["subtype"].astype(str).to_numpy(),
    },
    index=pd.Index([f"cell_{i:06d}" for i in range(len(selected))], name="cell_id"),
)
for col in atlas_obs.columns:
    atlas_obs[col] = atlas_obs[col].astype("category")

var = pd.DataFrame(index=pd.Index(gene_names[panel].astype(str), name="gene_symbol"))
adata = AnnData(X=counts_matrix[:, panel].tocsr(), obs=atlas_obs, var=var)
adata.layers["counts"] = adata.X.copy()

mdata = MuData({"rna": adata})
mdata.write_h5mu(f"{OUT}/atlas.h5mu", compression="gzip")
size_mb = os.path.getsize(f"{OUT}/atlas.h5mu") / 1e6
print(f"Wrote {OUT}/atlas.h5mu  ({adata.n_obs} nuclei x {adata.n_vars} genes, "
      f"{atlas_obs.participant_id.nunique()} donors, {size_mb:.1f} MB)")
PYCODE

echo "Done. Contents of $OUT:"
ls -la "$OUT"
