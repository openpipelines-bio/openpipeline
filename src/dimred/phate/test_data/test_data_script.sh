#!/usr/bin/env bash
# Generate synthetic test data for dimred/phate.
# Run from the repo root:
#   bash src/dimred/phate/test_data/test_data_script.sh

set -eo pipefail

python3 - <<'PYEOF'
import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData
import pathlib

out_dir = pathlib.Path("src/dimred/phate/test_data")
out_dir.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(42)
n_obs = 150
n_pcs = 20
n_genes = 200

# Simulate a trajectory: points spiralling in PCA space
t = np.linspace(0, 4 * np.pi, n_obs)
curve = np.column_stack([np.sin(t), np.cos(t), t / (4 * np.pi)])
X_pca = curve @ rng.standard_normal((3, n_pcs))
X_pca += rng.standard_normal((n_obs, n_pcs)) * 0.15

obs = pd.DataFrame(
    {"cell_type": np.repeat(["TypeA", "TypeB", "TypeC"], n_obs // 3)},
    index=[f"cell_{i}" for i in range(n_obs)],
)
var = pd.DataFrame(index=[f"gene_{g}" for g in range(n_genes)])
X = rng.integers(0, 1000, size=(n_obs, n_genes)).astype(float)

adata = AnnData(X=X, obs=obs, var=var)
adata.obsm["X_pca"] = X_pca.astype(np.float32)

mdata = MuData({"rna": adata})
out_file = out_dir / "phate_input.h5mu"
mdata.write_h5mu(str(out_file))
print(f"Written: {out_file}  ({n_obs} cells, {n_pcs} PCs)")
PYEOF

echo ""
echo "Run the component manually:"
echo ""
echo "  viash run src/dimred/phate/config.vsh.yaml -- \\"
echo "    --input  src/dimred/phate/test_data/phate_input.h5mu \\"
echo "    --output src/dimred/phate/test_data/phate_output.h5mu"
