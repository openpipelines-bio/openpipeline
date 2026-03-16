#!/usr/bin/env bash
# Generate synthetic test data for metadata/calculate_proportions.
# Run from the repo root:
#   bash src/metadata/calculate_proportions/test_data/test_data_script.sh

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 - <<'PYEOF'
import sys
import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData
import pathlib

out_dir = pathlib.Path("src/metadata/calculate_proportions/test_data")
out_dir.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(42)

participants = ["donor_A", "donor_B", "donor_C"]
subpopulations = ["ExN_1", "ExN_2", "InN_1", "Astro_1", "Micro_1"]

# Unequal cell counts per group to make proportions more interesting
cell_counts = {
    ("donor_A", "ExN_1"):  25,
    ("donor_A", "ExN_2"):  15,
    ("donor_A", "InN_1"):  20,
    ("donor_A", "Astro_1"): 30,
    ("donor_A", "Micro_1"): 10,
    ("donor_B", "ExN_1"):  40,
    ("donor_B", "ExN_2"):  10,
    ("donor_B", "InN_1"):  15,
    ("donor_B", "Astro_1"): 20,
    ("donor_B", "Micro_1"):  5,
    ("donor_C", "ExN_1"):  10,
    ("donor_C", "ExN_2"):  30,
    ("donor_C", "InN_1"):  25,
    ("donor_C", "Astro_1"): 15,
    ("donor_C", "Micro_1"): 20,
}

obs_rows = []
for (pid, subpop), n in cell_counts.items():
    for _ in range(n):
        obs_rows.append({"participant_id": pid, "subpopulation": subpop})

obs = pd.DataFrame(obs_rows)
obs.index = [f"cell_{k}" for k in range(len(obs))]

n_cells = len(obs)
n_genes = 100
X = rng.integers(0, 500, size=(n_cells, n_genes)).astype(float)
var = pd.DataFrame(index=[f"gene_{g}" for g in range(n_genes)])

adata = AnnData(X=X, obs=obs, var=var)
mdata = MuData({"rna": adata})

out_file = out_dir / "proportions_input.h5mu"
mdata.write_h5mu(str(out_file))
print(f"Written: {out_file}  ({n_cells} cells, {len(participants)} donors, {len(subpopulations)} subpopulations)")
PYEOF

echo ""
echo "Run the component manually:"
echo ""
echo "  viash run src/metadata/calculate_proportions/config.vsh.yaml -- \\"
echo "    --input  src/metadata/calculate_proportions/test_data/proportions_input.h5mu \\"
echo "    --output src/metadata/calculate_proportions/test_data/proportions_output.h5mu"
