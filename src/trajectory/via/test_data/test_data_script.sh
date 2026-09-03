#!/usr/bin/env bash
# Generate synthetic test data for trajectory/via.
# Run from repo root: bash src/trajectory/via/test_data/test_data_script.sh
set -eo pipefail

python3 - <<'PYEOF'
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu
import pathlib

out_dir = pathlib.Path("src/trajectory/via/test_data")
out_dir.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(0)
n_cells = 300
n_clusters = 6

# Linear trajectory embedded in 10-D (VIA needs >2 dims for graph construction)
pt = np.linspace(0, 1, n_cells)
# Primary component encodes pseudotime; remaining 9 are correlated noise
embedding = np.column_stack(
    [pt] + [pt * rng.uniform(0.1, 0.9) + rng.normal(0, 0.05, n_cells) for _ in range(9)]
).astype("float32")

cluster_ids = (pt * n_clusters).astype(int).clip(0, n_clusters - 1)
cluster_labels = [str(c) for c in cluster_ids]

obs = pd.DataFrame(
    {"leiden": cluster_labels},
    index=[f"c{i}" for i in range(n_cells)],
)
adata = ad.AnnData(
    X=rng.integers(0, 100, (n_cells, 5)).astype("float32"),
    obs=obs,
    var=pd.DataFrame(index=[f"g{j}" for j in range(5)]),
)
adata.obsm["X_phate"] = embedding

out = out_dir / "via_input.h5mu"
mu.MuData({"rna": adata}).write_h5mu(str(out))
print(f"Written: {out}  ({n_cells} cells, {n_clusters} clusters)")
PYEOF

echo ""
echo "Run manually:"
echo "  viash run src/trajectory/via/config.vsh.yaml -- \\"
echo "    --input  \$(pwd)/src/trajectory/via/test_data/via_input.h5mu \\"
echo "    --root_user 0 \\"
echo "    --output \$(pwd)/src/trajectory/via/test_data/via_output.h5mu"
