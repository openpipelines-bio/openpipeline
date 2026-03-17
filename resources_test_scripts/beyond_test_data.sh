#!/bin/bash
# Generates synthetic multi-donor snRNA-seq atlas for testing BEYOND components.
#
# Produces: resources_test/beyond_test_data/
#   atlas.h5mu              — 144 cells × 2000 genes; 8 donors; 3 cell types;
#                             9 subpopulations (3 per type); obs/obsm ready for
#                             BEYOND components (calculate_proportions, phate, etc.)
#   proportions_output.h5mu — atlas + uns["proportions"] + obsm["proportions"]
#                             (simulates output of metadata/calculate_proportions)
#   pseudotime_output.h5mu  — proportions_output + obs["palantir_pseudotime"]
#                             (simulates output of trajectory/palantir)
#   dynamics_output.h5mu    — pseudotime_output + uns["dynamics"]
#                             (simulates output of trajectory/pseudotime_dynamics)
#   de_ExN.csv              — synthetic DESeq2 results for ExN subpopulations
#   de_InN.csv              — synthetic DESeq2 results for InN subpopulations
#   de_Ast.csv              — synthetic DESeq2 results for Ast subpopulations
#   gene_sets.gmt           — tiny local GMT file (6 gene sets) for pathway_enrichment
#   traits.csv              — synthetic clinical traits table (participant × trait)
#
# Usage: bash resources_test_scripts/beyond_test_data.sh

set -eo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

OUT="resources_test/beyond_test_data"
mkdir -p "$OUT"

python3 - <<'PYEOF'
import os
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA

out = "resources_test/beyond_test_data"
rng = np.random.default_rng(42)

# ── Parameters ───────────────────────────────────────────────────────────────
N_DONORS     = 8
CELL_TYPES   = ["ExN", "InN", "Ast"]
SUBPOPS      = {ct: [f"{ct}.{i}" for i in range(1, 4)] for ct in CELL_TYPES}
CELLS_PER_DONOR_PER_SUBPOP = 2   # 2 × 3 subpops × 3 cell types × 8 donors = 144 cells
N_GENES      = 2000

gene_names   = [f"GENE{i:05d}" for i in range(N_GENES)]
donors       = [f"donor_{i:02d}" for i in range(1, N_DONORS + 1)]

# ── Simulate expression profiles per cell type ────────────────────────────────
# Each cell type has a distinct 100-gene signature
sig_size = 100
ct_signature = {}
for i, ct in enumerate(CELL_TYPES):
    ct_signature[ct] = list(range(i * sig_size, (i + 1) * sig_size))

rows = []
for ct in CELL_TYPES:
    for sp in SUBPOPS[ct]:
        sp_idx = SUBPOPS[ct].index(sp)
        for donor in donors:
            for _ in range(CELLS_PER_DONOR_PER_SUBPOP):
                expr = rng.negative_binomial(2, 0.5, N_GENES).astype("float32")
                expr[ct_signature[ct]] += rng.negative_binomial(8, 0.3, sig_size)
                sp_sig = list(range(300 + sp_idx * 10, 300 + (sp_idx + 1) * 10))
                expr[sp_sig] += rng.negative_binomial(5, 0.4, 10)
                rows.append({
                    "counts": expr,
                    "cell_type": ct,
                    "subpopulation": sp,
                    "participant_id": donor,
                    "sample_id": donor,
                })

# Build AnnData
n_cells = len(rows)
X = csr_matrix(np.vstack([r["counts"] for r in rows]))
obs = pd.DataFrame([{k: v for k, v in r.items() if k != "counts"} for r in rows])
obs.index = [f"cell_{i:05d}" for i in range(n_cells)]

for col in ["cell_type", "subpopulation", "participant_id", "sample_id"]:
    obs[col] = pd.Categorical(obs[col])

obs["n_counts"]        = np.asarray(X.sum(axis=1)).ravel().astype(int)
obs["n_genes"]         = np.asarray((X > 0).sum(axis=1)).ravel().astype(int)
obs["scrublet_score"]  = rng.uniform(0, 0.3, n_cells)
obs["doublet_called"]  = obs["scrublet_score"] > 0.25
obs["leiden"]          = pd.Categorical(obs["subpopulation"])
obs["celltypist_pred"] = pd.Categorical(obs["cell_type"])

var = pd.DataFrame(index=pd.Index(gene_names, name="gene_symbol"))
var["highly_variable"] = False
var.loc[var.index[:500], "highly_variable"] = True

adata = ad.AnnData(X=X, obs=obs, var=var)

# ── PCA using sklearn (avoids scanpy/numba/numpy-2 issues) ────────────────────
X_dense = np.asarray(X.todense(), dtype="float32")
# log-normalise
row_sums = X_dense.sum(axis=1, keepdims=True) + 1e-9
X_lognorm = np.log1p(X_dense / row_sums * 1e4)

pca = PCA(n_components=30, random_state=42)
X_pca = pca.fit_transform(X_lognorm[:, :500]).astype("float32")  # HVG subset
adata.obsm["X_pca"] = X_pca

# Harmony-like correction: small donor shift
donor_shift = {d: rng.normal(0, 0.3, 30) for d in donors}
X_integrated = X_pca.copy()
for i, donor in enumerate(obs["participant_id"]):
    X_integrated[i] += donor_shift[donor]
adata.obsm["X_pca_integrated"] = X_integrated.astype("float32")

# Simple 2-D UMAP proxy: project first 2 PCs + noise (no internet/numba needed)
X_umap = X_integrated[:, :2].copy()
X_umap += rng.normal(0, 0.05, X_umap.shape)
adata.obsm["X_umap"] = X_umap.astype("float32")

# ── Write atlas h5mu ──────────────────────────────────────────────────────────
mdata = mu.MuData({"rna": adata})
atlas_path = f"{out}/atlas.h5mu"
mdata.write_h5mu(atlas_path, compression="gzip")
print(f"Wrote {atlas_path}  ({n_cells} cells × {N_GENES} genes, {N_DONORS} donors)")

# ── Synthetic DESeq2-like CSVs (one per cell type) ────────────────────────────
for ct in CELL_TYPES:
    lfc    = rng.normal(0, 2, N_GENES)
    pvals  = rng.uniform(0, 1, N_GENES)
    rank   = np.argsort(pvals).argsort() + 1
    padj   = np.clip(pvals * N_GENES / rank, 0, 1)
    de = pd.DataFrame({
        "log2FoldChange": lfc,
        "pvalue":         pvals,
        "padj":           padj,
        "baseMean":       rng.uniform(10, 1000, N_GENES),
        "lfcSE":          rng.uniform(0.1, 0.5, N_GENES),
    }, index=pd.Index(gene_names, name="gene"))
    csv_path = f"{out}/de_{ct}.csv"
    de.to_csv(csv_path)
    nsig = (padj < 0.05).sum()
    print(f"Wrote {csv_path}  ({N_GENES} genes, {nsig} sig at padj<0.05)")

# ── Local GMT file (6 gene sets, entirely from synthetic genes) ────────────────
gmt_path = f"{out}/gene_sets.gmt"
with open(gmt_path, "w") as fh:
    fh.write("SIG_A\tna\t"  + "\t".join(gene_names[  0: 40]) + "\n")
    fh.write("SIG_B\tna\t"  + "\t".join(gene_names[ 40:100]) + "\n")
    fh.write("SIG_C\tna\t"  + "\t".join(gene_names[100:160]) + "\n")
    fh.write("SIG_ExN\tna\t" + "\t".join(gene_names[  0:100]) + "\n")
    fh.write("SIG_InN\tna\t" + "\t".join(gene_names[100:200]) + "\n")
    fh.write("SIG_Ast\tna\t" + "\t".join(gene_names[200:300]) + "\n")
print(f"Wrote {gmt_path}  (6 gene sets)")

# ── Clinical traits CSV (participant_id × traits) ─────────────────────────────
traits = pd.DataFrame({
    "participant_id": donors,
    "age":            rng.integers(60, 90, N_DONORS).astype(float),
    "sex":            rng.choice(["M", "F"], N_DONORS),
    "diagnosis":      rng.choice(["control", "AD"], N_DONORS, p=[0.5, 0.5]),
    "cohort":         rng.choice(["cohort_A", "cohort_B"], N_DONORS),
    "pmi":            rng.uniform(2, 24, N_DONORS).round(1),
})
traits_path = f"{out}/traits.csv"
traits.to_csv(traits_path, index=False)
print(f"Wrote {traits_path}  ({N_DONORS} donors × 5 traits)")

# ── Proportions output (simulates metadata/calculate_proportions output) ──────
counts_df = (
    adata.obs.groupby(["participant_id", "subpopulation"], observed=True)
    .size()
    .unstack(fill_value=0)
)
proportions_df = counts_df.div(counts_df.sum(axis=1), axis=0)
adata.uns["proportions"] = proportions_df.to_dict()

# Per-cell proportion vectors in obsm (each cell gets its participant's row)
participant_ids_arr = adata.obs["participant_id"].values
obsm_matrix = np.array(
    [
        proportions_df.loc[pid].values if pid in proportions_df.index
        else np.zeros(len(proportions_df.columns))
        for pid in participant_ids_arr
    ],
    dtype=np.float64,
)
adata.obsm["proportions"] = pd.DataFrame(
    obsm_matrix,
    index=adata.obs_names,
    columns=proportions_df.columns.astype(str),
)
prop_path = f"{out}/proportions_output.h5mu"
mu.MuData({"rna": adata}).write_h5mu(prop_path, compression="gzip")
print(
    f"Wrote {prop_path}  "
    f"({len(proportions_df)} participants × {len(proportions_df.columns)} subpopulations)"
)

# ── Pseudotime output (simulates trajectory/palantir output) ──────────────────
donor_order = sorted(proportions_df.index)
donor_pt = {
    d: float(i) / max(len(donor_order) - 1, 1)
    for i, d in enumerate(donor_order)
}
adata.obs["palantir_pseudotime"] = [
    float(np.clip(donor_pt[pid] + rng.normal(0, 0.02), 0, 1))
    for pid in adata.obs["participant_id"]
]
pt_path = f"{out}/pseudotime_output.h5mu"
mu.MuData({"rna": adata}).write_h5mu(pt_path, compression="gzip")
print(f"Wrote {pt_path}  ({len(donor_order)} donors with synthetic palantir_pseudotime)")

# ── Dynamics output (simulates trajectory/pseudotime_dynamics output) ─────────
n_bins = 100
grid = np.linspace(0, 1, n_bins)
subpop_list = list(proportions_df.columns)

def _sigmoid(x, x0, k=10.0):
    return 1.0 / (1.0 + np.exp(-k * (x - x0)))

dynamics = {}
for j, sp in enumerate(subpop_list):
    peak_t = 0.1 + j * (0.8 / max(len(subpop_list) - 1, 1))
    fitted = _sigmoid(grid, peak_t).tolist()
    dynamics[sp] = {
        "pseudotime_grid": grid.tolist(),
        "proportion_fitted": fitted,
        "peak_pseudotime": float(grid[int(np.argmax(fitted))]),
        "r_squared": 0.85,
        "p_value": 0.001,
    }
adata.uns["dynamics"] = dynamics
dyn_path = f"{out}/dynamics_output.h5mu"
mu.MuData({"rna": adata}).write_h5mu(dyn_path, compression="gzip")
print(f"Wrote {dyn_path}  ({len(subpop_list)} subpopulation dynamics curves)")
PYEOF

echo ""
echo "Done. Test data in $OUT:"
ls -lh "$OUT"
