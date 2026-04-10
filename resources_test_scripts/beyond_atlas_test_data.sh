#!/bin/bash
# Generates synthetic per-sample h5mu files and a reference h5mu for testing
# the beyond/atlas_building workflow.
#
# Produces: resources_test/beyond_atlas_test/
#   donor_01.h5mu … donor_04.h5mu  — 180 cells × 500 genes, raw counts, 3 cell types
#                                     var["gene_symbol"]: GENE00000…GENE00486 + 13 MT- genes
#   reference.h5mu                 — 90 cells × 500 genes, raw counts, obs["cell_type"]
#                                     same var["gene_symbol"] layout as donor files
#
# Usage: bash resources_test_scripts/beyond_atlas_test_data.sh

set -eo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

OUT="resources_test/beyond_atlas_test"
mkdir -p "$OUT"

python3 - <<'PYEOF'
import os
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu
from scipy.sparse import csr_matrix

out = "resources_test/beyond_atlas_test"
rng = np.random.default_rng(42)

N_GENES = 500
# Last 13 genes are mitochondrial (MT- prefix), matching the default mito regex
MT_GENES = [
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
    "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-CYB",
]
N_MT = len(MT_GENES)
gene_names  = [f"GENE{i:05d}" for i in range(N_GENES - N_MT)] + MT_GENES
gene_symbols = gene_names  # index and gene_symbol are the same here

CELL_TYPES = ["ExN", "InN", "Ast"]
# Cell-type signatures: genes 0-49 for ExN, 50-99 for InN, 100-149 for Ast
CT_SIGNAL = {
    "ExN": (0,   50),
    "InN": (50,  100),
    "Ast": (100, 150),
}

DONORS = [f"donor_{i:02d}" for i in range(1, 5)]   # 4 donors
CELLS_PER_CT = 60                                    # 60 cells per cell type per donor (≥50 needed for PCA n_comps=50)

def make_counts(n_cells, cell_type, rng):
    """Simulate negative-binomial raw counts with a cell-type signal and ~8% mito."""
    X = rng.negative_binomial(2, 0.5, (n_cells, N_GENES)).astype("float32")
    lo, hi = CT_SIGNAL[cell_type]
    X[:, lo:hi] += rng.negative_binomial(8, 0.3, (n_cells, hi - lo))
    # Boost MT gene counts so mito fraction ≈ 5–15% (well under max_fraction_mito=0.2)
    X[:, -N_MT:] += rng.negative_binomial(5, 0.4, (n_cells, N_MT))
    return X

# ── Per-donor files (donor_01 … donor_04) ─────────────────────────────────────
for donor in DONORS:
    rows_X   = []
    cell_ids = []
    cell_type_col = []

    for ct in CELL_TYPES:
        X_ct = make_counts(CELLS_PER_CT, ct, rng)
        rows_X.append(X_ct)
        for i in range(CELLS_PER_CT):
            cell_ids.append(f"{donor}_{ct}_{i:03d}")
            cell_type_col.append(ct)

    X_mat = csr_matrix(np.vstack(rows_X))
    n_cells = X_mat.shape[0]

    obs = pd.DataFrame({
        "cell_type": pd.Categorical(cell_type_col, categories=CELL_TYPES),
    }, index=pd.Index(cell_ids, name="obs_names"))

    var = pd.DataFrame(
        {"gene_symbol": gene_symbols},
        index=pd.Index(gene_names, name="var_names"),
    )

    adata = ad.AnnData(X=X_mat, obs=obs, var=var)
    mdata = mu.MuData({"rna": adata})

    path = f"{out}/{donor}.h5mu"
    mdata.write_h5mu(path, compression="gzip")
    print(f"Wrote {path}  ({n_cells} cells x {N_GENES} genes, cell types: {CELL_TYPES})")

# ── Reference file (90 cells, 30 per cell type) ───────────────────────────────
ref_rows_X   = []
ref_cell_ids = []
ref_ct_col   = []

N_REF_PER_CT = 30
for ct in CELL_TYPES:
    X_ct = make_counts(N_REF_PER_CT, ct, rng)
    ref_rows_X.append(X_ct)
    for i in range(N_REF_PER_CT):
        ref_cell_ids.append(f"ref_{ct}_{i:03d}")
        ref_ct_col.append(ct)

X_ref = csr_matrix(np.vstack(ref_rows_X))
obs_ref = pd.DataFrame({
    "cell_type": pd.Categorical(ref_ct_col, categories=CELL_TYPES),
}, index=pd.Index(ref_cell_ids, name="obs_names"))
var_ref = pd.DataFrame(
    {"gene_symbol": gene_symbols},
    index=pd.Index(gene_names, name="var_names"),
)

adata_ref = ad.AnnData(X=X_ref, obs=obs_ref, var=var_ref)
mdata_ref = mu.MuData({"rna": adata_ref})

ref_path = f"{out}/reference.h5mu"
mdata_ref.write_h5mu(ref_path, compression="gzip")
print(f"Wrote {ref_path}  ({X_ref.shape[0]} cells x {N_GENES} genes, 3 cell types)")
PYEOF

echo ""
echo "Done. Test data in $OUT:"
ls -lh "$OUT"
