#!/bin/bash
# Generates a simulated multi-donor snRNA-seq atlas for testing the BEYOND components.
#
# This is the companion of `beyond_trajectory_test_data.sh`, which builds a fixture from
# a real cohort. The two answer different questions and neither replaces the other:
#
#   - the real fixture (PsychAD RADC, 152 donors) says whether the workflow survives real
#     annotation hierarchies, real sparsity and real donor metadata;
#   - this simulation says whether the workflow *recovers* what is in the data, because
#     here we know what is in it.
#
# The real cohort cannot answer the second question. Scanning its 152 donors over 3
# annotation levels x 5 donor traits finds no compositional association surviving BH
# correction (closest: AD_status x OPC at class level, q = 0.053); the BEYOND paper's own
# result rests on 437 donors. A fixture with no signal can only test plumbing.
#
# The simulation is structured, not noise: every quantity the BEYOND workflow is supposed
# to recover is put into the data on purpose, so a broken component changes the result.
#
#   - Each donor has a latent severity s in [0, 1]. Donors are ordered by s, which is what
#     PHATE + Palantir are expected to recover as a trajectory.
#   - Subpopulation composition is a function of s: within each cell type one subpopulation
#     expands with severity, one contracts, one is flat. "Proportions" throughout means, per
#     donor, the fraction of that donor's cells that fall in each subpopulation (rows sum
#     to 1) - that is the matrix `stats/calculate_label_proportions` produces and the input
#     of PHATE, the dynamics fit, the communities and the association tests.
#   - Counts come from a negative-binomial model with per-gene means, per-cell library
#     sizes, cell-type and subpopulation marker genes, a per-donor gene-level effect
#     (applied at the count level, so integration has something real to remove) and a
#     disease effect on a fixed set of genes, scaled by s.
#   - The DE tables are computed from the simulated counts (donor pseudobulk, Welch t-test
#     of high- vs low-severity donors, BH correction) rather than drawn at random, so the
#     genes they call are the genes the simulation perturbed, and they overlap the
#     DISEASE_UP / DISEASE_DOWN sets in the GMT file.
#   - The traits table carries `amyloid` and `braak`, both functions of s, plus traits that
#     are independent of it - so an association test has both a true and a null answer.
#
# Produces: resources_test/beyond_simulated_test_data/
#   atlas.h5mu      - ~864 cells x 2000 genes; 12 donors; 3 cell types; 9 subpopulations;
#                     obsm["X_pca"], obsm["X_pca_integrated"]
#   de_ExN.csv      - DE results per cell type, computed from the counts above
#   de_InN.csv
#   de_Ast.csv
#   gene_sets.gmt   - local GMT (no internet needed), incl. the true disease gene sets
#   traits.csv      - donor-level traits (participant x trait)
#
# Usage: bash resources_test_scripts/beyond_trajectory_test_data.sh

set -eo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

OUT="resources_test/beyond_simulated_test_data"
mkdir -p "$OUT"

python3 - <<'PYEOF'
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu
from scipy.sparse import csr_matrix
from scipy.stats import t as t_dist
from sklearn.decomposition import PCA

out = "resources_test/beyond_simulated_test_data"
rng = np.random.default_rng(42)

# -- Parameters ---------------------------------------------------------------
N_DONORS   = 12
CELL_TYPES = ["ExN", "InN", "Ast"]
SUBPOPS    = {ct: [f"{ct}.{i}" for i in range(1, 4)] for ct in CELL_TYPES}
N_GENES    = 2000
CELLS_PER_DONOR_PER_TYPE = 24   # split over the 3 subpopulations by composition
NB_DISPERSION = 5.0             # negative-binomial size parameter

gene_names = [f"GENE{i:05d}" for i in range(N_GENES)]
donors     = [f"donor_{i:02d}" for i in range(1, N_DONORS + 1)]

# Latent severity per donor: the trajectory the workflow should recover
severity = dict(zip(donors, np.linspace(0.0, 1.0, N_DONORS)))

# -- Gene programmes -----------------------------------------------------------
# Baseline expression level per gene (long-tailed, as in real data)
base_mean = rng.gamma(shape=0.8, scale=6.0, size=N_GENES) + 0.1

# Cell-type marker genes: 100 per type, 4x up in that type
ct_markers = {ct: np.arange(i * 100, (i + 1) * 100) for i, ct in enumerate(CELL_TYPES)}
# Subpopulation marker genes: 20 per subpopulation, 3x up in that subpopulation
sp_markers = {}
offset = 300
for ct in CELL_TYPES:
    for sp in SUBPOPS[ct]:
        sp_markers[sp] = np.arange(offset, offset + 20)
        offset += 20

# Disease genes: up- and down-regulated with severity, in every cell type
disease_up   = np.arange(1000, 1060)
disease_down = np.arange(1060, 1120)

# Per-donor gene-level effect, applied to the counts (not to the embedding)
donor_effect = {d: np.exp(rng.normal(0, 0.12, N_GENES)) for d in donors}


def composition(cell_type, s):
    """Subpopulation fractions within a cell type at severity s: one up, one down, one flat."""
    weights = np.array([0.2 + 0.6 * s, 0.8 - 0.6 * s, 0.4])
    return weights / weights.sum()


# -- Simulate cells ------------------------------------------------------------
counts_rows, obs_rows = [], []
for donor in donors:
    s = severity[donor]
    for ct in CELL_TYPES:
        fractions = composition(ct, s)
        # Multinomial over the subpopulations of this cell type, so composition
        # varies with severity and cell counts still vary between donors.
        n_per_sp = rng.multinomial(CELLS_PER_DONOR_PER_TYPE, fractions)
        for sp, n_cells in zip(SUBPOPS[ct], n_per_sp):
            for _ in range(int(n_cells)):
                mu_gene = base_mean.copy()
                mu_gene[ct_markers[ct]] *= 4.0
                mu_gene[sp_markers[sp]] *= 3.0
                mu_gene[disease_up]   *= 1.0 + 8.0 * s
                mu_gene[disease_down] *= 1.0 / (1.0 + 8.0 * s)
                mu_gene = mu_gene * donor_effect[donor]
                # Library size varies per cell
                mu_gene = mu_gene * rng.lognormal(0.0, 0.2)
                # Negative binomial with fixed dispersion
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

obs["n_counts"]        = np.asarray(X.sum(axis=1)).ravel().astype(int)
obs["n_genes"]         = np.asarray((X > 0).sum(axis=1)).ravel().astype(int)
obs["leiden"]          = pd.Categorical(obs["subpopulation"])
obs["celltypist_pred"] = pd.Categorical(obs["cell_type"])

var = pd.DataFrame(index=pd.Index(gene_names, name="gene_symbol"))
var["highly_variable"] = False
var.loc[var.index[:600], "highly_variable"] = True

adata = ad.AnnData(X=X, obs=obs, var=var)

# -- Embeddings ----------------------------------------------------------------
# PCA with sklearn on the log-normalised highly variable genes. sklearn rather than
# scanpy so this script needs only numpy/scipy/sklearn/anndata/mudata; the result is
# the same decomposition scanpy would compute.
X_dense   = np.asarray(X.todense(), dtype="float32")
X_lognorm = np.log1p(X_dense / (X_dense.sum(axis=1, keepdims=True) + 1e-9) * 1e4)
hvg       = np.flatnonzero(var["highly_variable"].to_numpy())

X_pca = PCA(n_components=30, random_state=42).fit_transform(X_lognorm[:, hvg])
adata.obsm["X_pca"] = X_pca.astype("float32")

# The donor effect is in the counts, so it is in X_pca. "Integrated" means it has been
# removed: centre each donor's cells on the global mean (a stand-in for Harmony).
X_integrated = X_pca.copy()
participants = obs["participant_id"].to_numpy()
for donor in donors:
    mask = participants == donor
    X_integrated[mask] += X_pca.mean(axis=0) - X_pca[mask].mean(axis=0)
adata.obsm["X_pca_integrated"] = X_integrated.astype("float32")

# No .obsm["X_umap"]: nothing in the BEYOND components reads it, and a fake 2-D
# projection would only look like a real UMAP.

mdata = mu.MuData({"rna": adata})
atlas_path = f"{out}/atlas.h5mu"
mdata.write_h5mu(atlas_path, compression="gzip")
print(f"Wrote {atlas_path}  ({n_cells} cells x {N_GENES} genes, {N_DONORS} donors)")

# -- DE tables, computed from the simulated counts -----------------------------
# Donor pseudobulk per cell type, then a linear regression of log2 CPM on donor severity
# across all donors (which uses the whole gradient, not a high/low split). The genes this
# calls are the genes the simulation perturbed. `log2FoldChange` is the slope, i.e. the
# log2 change between the least and the most affected donor.
sev_vec = np.array([severity[d] for d in donors])
x_centered = sev_vec - sev_vec.mean()
dof = N_DONORS - 2

for ct in CELL_TYPES:
    ct_mask = (obs["cell_type"] == ct).to_numpy()
    pseudobulk = {}
    for donor in donors:
        mask = ct_mask & (participants == donor)
        summed = X_dense[mask].sum(axis=0)
        cpm = summed / max(summed.sum(), 1.0) * 1e6
        pseudobulk[donor] = np.log2(cpm + 1.0)
    pb = pd.DataFrame(pseudobulk, index=gene_names).T

    Y = pb.to_numpy()
    Y_centered = Y - Y.mean(axis=0)
    slope = (x_centered @ Y_centered) / (x_centered @ x_centered)
    fitted = np.outer(x_centered, slope)
    resid_var = ((Y_centered - fitted) ** 2).sum(axis=0) / dof
    se = np.sqrt(resid_var / (x_centered @ x_centered)) + 1e-12
    stat = slope / se
    pvals = 2 * t_dist.sf(np.abs(stat), dof)
    pvals = np.nan_to_num(pvals, nan=1.0)
    lfc = slope

    order = np.argsort(pvals)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, N_GENES + 1)
    padj = np.minimum.accumulate(
        (pvals * N_GENES / ranks)[order][::-1]
    )[::-1]
    padj_full = np.empty(N_GENES)
    padj_full[order] = np.clip(padj, 0, 1)

    de = pd.DataFrame(
        {
            "baseMean": pb.mean(axis=0).to_numpy(),
            "log2FoldChange": lfc,
            "lfcSE": se,
            "stat": np.nan_to_num(stat),
            "pvalue": pvals,
            "padj": padj_full,
        },
        index=pd.Index(gene_names, name="gene"),
    )
    csv_path = f"{out}/de_{ct}.csv"
    de.to_csv(csv_path)
    n_sig = int((padj_full < 0.05).sum())
    n_sig_disease = int(
        (padj_full[np.concatenate([disease_up, disease_down])] < 0.05).sum()
    )
    print(
        f"Wrote {csv_path}  ({N_GENES} genes, {n_sig} sig at padj<0.05, "
        f"{n_sig_disease}/120 of them in the simulated disease sets, "
        f"min padj {padj_full.min():.2g})"
    )

# -- Local GMT file ------------------------------------------------------------
gmt_path = f"{out}/gene_sets.gmt"
with open(gmt_path, "w") as fh:
    fh.write("DISEASE_UP\tna\t"   + "\t".join(gene_names[i] for i in disease_up) + "\n")
    fh.write("DISEASE_DOWN\tna\t" + "\t".join(gene_names[i] for i in disease_down) + "\n")
    for ct in CELL_TYPES:
        fh.write(
            f"MARKERS_{ct}\tna\t"
            + "\t".join(gene_names[i] for i in ct_markers[ct])
            + "\n"
        )
    fh.write("BACKGROUND\tna\t" + "\t".join(gene_names[1500:1600]) + "\n")
print(f"Wrote {gmt_path}  (6 gene sets, 2 of them the true disease sets)")

# -- Donor traits --------------------------------------------------------------
sev = np.array([severity[d] for d in donors])
traits = pd.DataFrame(
    {
        "participant_id": donors,
        # traits that follow the latent severity: an association test must find these
        "amyloid": (sev * 10 + rng.normal(0, 0.5, N_DONORS)).round(2),
        "braak":   np.clip(np.round(sev * 6 + rng.normal(0, 0.4, N_DONORS)), 0, 6),
        "diagnosis": np.where(sev > 0.5, "AD", "control"),
        # traits independent of it: an association test must not find these
        "age": rng.integers(60, 90, N_DONORS).astype(float),
        "sex": rng.choice(["M", "F"], N_DONORS),
        "pmi": rng.uniform(2, 24, N_DONORS).round(1),
        "cohort": rng.choice(["cohort_A", "cohort_B"], N_DONORS),
    }
)
traits_path = f"{out}/traits.csv"
traits.to_csv(traits_path, index=False)
print(f"Wrote {traits_path}  ({N_DONORS} donors x {traits.shape[1] - 1} traits)")

# -- Report the composition signal that was built in --------------------------
counts_df = (
    adata.obs.groupby(["participant_id", "subpopulation"], observed=True)
    .size()
    .unstack(fill_value=0)
)
proportions = counts_df.div(counts_df.sum(axis=1), axis=0)
corr = proportions.corrwith(pd.Series(severity), axis=0)
print(
    "Subpopulation proportion vs severity correlation: "
    + ", ".join(f"{sp}={corr[sp]:+.2f}" for sp in proportions.columns)
)
PYEOF

echo ""
echo "Done. Test data in $OUT:"
ls -lh "$OUT"
