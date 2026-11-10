"""Simulated 12-donor snRNA-seq atlas with a planted severity gradient, for the BEYOND tests.

Each donor has a latent severity s in [0, 1]. Subpopulation composition, a set of
disease genes and the traits `amyloid` / `braak` follow s; `age`, `sex`, `pmi` and
`cohort` do not. Counts are negative binomial with a per-donor gene effect.

Produces resources_test/beyond_simulated_test_data/:
  atlas.h5mu                          864 cells x 2000 genes, 12 donors, 3 cell types,
                                      9 subpopulations
  traits.csv                          participant x trait

Usage: python3 resources_test_scripts/beyond_trajectory_simulated_test_data.py
"""

from pathlib import Path

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

out = Path(__file__).resolve().parents[1] / "resources_test/beyond_simulated_test_data"
out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)

N_DONORS = 12
CELL_TYPES = ["ExN", "InN", "Ast"]
SUBPOPS = {ct: [f"{ct}.{i}" for i in range(1, 4)] for ct in CELL_TYPES}
N_GENES = 2000
CELLS_PER_DONOR_PER_TYPE = 24
NB_DISPERSION = 5.0

gene_names = [f"GENE{i:05d}" for i in range(N_GENES)]
donors = [f"donor_{i:02d}" for i in range(1, N_DONORS + 1)]
severity = dict(zip(donors, np.linspace(0.0, 1.0, N_DONORS)))

# -- Gene programmes -----------------------------------------------------------
base_mean = rng.gamma(shape=0.8, scale=6.0, size=N_GENES) + 0.1
# 100 cell-type markers (4x) per type, 20 subpopulation markers (3x) per subpopulation
ct_markers = {ct: np.arange(i * 100, (i + 1) * 100) for i, ct in enumerate(CELL_TYPES)}
sp_markers = {}
offset = 300
for ct in CELL_TYPES:
    for sp in SUBPOPS[ct]:
        sp_markers[sp] = np.arange(offset, offset + 20)
        offset += 20
disease_up = np.arange(1000, 1060)
disease_down = np.arange(1060, 1120)
donor_effect = {d: np.exp(rng.normal(0, 0.12, N_GENES)) for d in donors}


def composition(s):
    """Subpopulation fractions at severity s: one up, one down, one flat."""
    weights = np.array([0.2 + 0.6 * s, 0.8 - 0.6 * s, 0.4])
    return weights / weights.sum()


# -- Simulate cells ------------------------------------------------------------
counts_rows, obs_rows = [], []
for donor in donors:
    s = severity[donor]
    for ct in CELL_TYPES:
        n_per_sp = rng.multinomial(CELLS_PER_DONOR_PER_TYPE, composition(s))
        for sp, n_cells in zip(SUBPOPS[ct], n_per_sp):
            for _ in range(int(n_cells)):
                mu_gene = base_mean.copy()
                mu_gene[ct_markers[ct]] *= 4.0
                mu_gene[sp_markers[sp]] *= 3.0
                mu_gene[disease_up] *= 1.0 + 8.0 * s
                mu_gene[disease_down] *= 1.0 / (1.0 + 8.0 * s)
                mu_gene = mu_gene * donor_effect[donor]
                mu_gene = mu_gene * rng.lognormal(0.0, 0.2)  # library size
                p = NB_DISPERSION / (NB_DISPERSION + mu_gene)
                counts_rows.append(
                    rng.negative_binomial(NB_DISPERSION, p).astype("float32")
                )
                obs_rows.append(
                    {
                        "cell_type": ct,
                        "subpopulation": sp,
                        "participant_id": donor,
                        "sample_id": donor,
                    }
                )

n_cells = len(counts_rows)
X = csr_matrix(np.vstack(counts_rows))
obs = pd.DataFrame(obs_rows)
obs.index = [f"cell_{i:05d}" for i in range(n_cells)]
for col in ["cell_type", "subpopulation", "participant_id", "sample_id"]:
    obs[col] = pd.Categorical(obs[col])

obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel().astype(int)
obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel().astype(int)
obs["leiden"] = pd.Categorical(obs["subpopulation"])
obs["celltypist_pred"] = pd.Categorical(obs["cell_type"])

var = pd.DataFrame(index=pd.Index(gene_names, name="gene_symbol"))
var["highly_variable"] = False
var.loc[var.index[:600], "highly_variable"] = True

adata = ad.AnnData(X=X, obs=obs, var=var)
mdata = mu.MuData({"rna": adata})
mdata.write_h5mu(out / "atlas.h5mu", compression="gzip")
print(f"Wrote atlas.h5mu  ({n_cells} cells x {N_GENES} genes, {N_DONORS} donors)")

sev_vec = np.array([severity[d] for d in donors])

# -- Donor traits: amyloid / braak follow severity, the rest are null -------------
traits = pd.DataFrame(
    {
        "participant_id": donors,
        "amyloid": (sev_vec * 10 + rng.normal(0, 0.5, N_DONORS)).round(2),
        "braak": np.clip(np.round(sev_vec * 6 + rng.normal(0, 0.4, N_DONORS)), 0, 6),
        "diagnosis": np.where(sev_vec > 0.5, "AD", "control"),
        "age": rng.integers(60, 90, N_DONORS).astype(float),
        "sex": rng.choice(["M", "F"], N_DONORS),
        "pmi": rng.uniform(2, 24, N_DONORS).round(1),
        "cohort": rng.choice(["cohort_A", "cohort_B"], N_DONORS),
    }
)
traits.to_csv(out / "traits.csv", index=False)

counts_df = (
    obs.groupby(["participant_id", "subpopulation"], observed=True)
    .size()
    .unstack(fill_value=0)
)
proportions = counts_df.div(counts_df.sum(axis=1), axis=0)
corr = proportions.corrwith(pd.Series(severity), axis=0)
print("Proportion vs severity: " + ", ".join(f"{k}={v:+.2f}" for k, v in corr.items()))
