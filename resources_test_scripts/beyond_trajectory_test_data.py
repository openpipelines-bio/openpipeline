"""BEYOND test fixture from the PsychAD RADC cohort (152 donors, DLPFC snRNA-seq).

Source: PsychAD RADC_Cohort, 693682 nuclei x 34176 genes, CZ CELLxGENE Discover
  https://cellxgene.cziscience.com/collections/84ce6837-548d-4a1f-919f-0bc0d9a3952f
  doi:10.1101/2024.10.31.24316513. Licence CC-BY 4.0 (attribution required).

Produces resources_test/beyond_test_data/:
  atlas.h5mu                          ~200 nuclei per donor, stratified by subtype;
                                      ~300-gene panel (class markers + top expressed)
  proportions.csv                     donor x subtype, within-class prevalence, all nuclei
  traits.csv                          donor x trait

The 6.28 GB source is downloaded once to resources_test/beyond_test_data/; set
BEYOND_RADC_H5AD to use a copy already on disk.

Usage: python3 resources_test_scripts/beyond_trajectory_test_data.py
"""

import os
import urllib.request
from pathlib import Path

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import scanpy as sc

RADC_URL = "https://datasets.cellxgene.cziscience.com/54293783-669c-410e-919d-474960f8761b.h5ad"
SEED = 0
CELLS_PER_DONOR = 200
N_MARKERS_PER_CLASS = 25
N_TOP_EXPRESSED = 100

out = Path(__file__).resolve().parents[1] / "resources_test/beyond_test_data"
out.mkdir(parents=True, exist_ok=True)
src = Path(os.environ.get("BEYOND_RADC_H5AD", out / "temp_RADC_Cohort.h5ad"))
if not src.exists():
    print(f"Downloading PsychAD RADC_Cohort (6.28 GB) to {src}")
    src.parent.mkdir(parents=True, exist_ok=True)
    # urlretrieve raises on a short read; the rename keeps a partial file from being
    # taken for the source on the next run
    part = src.with_name(src.name + ".part")
    urllib.request.urlretrieve(RADC_URL, part)
    part.rename(src)

rng = np.random.default_rng(SEED)
full = ad.read_h5ad(src, backed="r")
gene_names = full.var["feature_name"].astype(str).to_numpy()
obs = full.obs[
    [
        "donor_id",
        "class",
        "subclass",
        "subtype",
        "sex",
        "genetic_ancestry",
        "development_stage",
        "AD_status",
        "Parkinson_disease",
        "Vascular_status",
        "DLBD_status",
        "Schizophrenia",
        "ASCVD_status",
    ]
].reset_index(drop=True)

# -- proportions: within-class prevalence over all nuclei -------------------------
# Normalised within class, as in BEYOND 2. Cell-type analysis/3.create.proportion.matrix.R
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
    index="donor_id",
    columns="subtype",
    values="prevalence",
    fill_value=0,
    aggfunc="sum",
    observed=True,
)
proportions.index.name = "participant_id"
proportions.to_csv(out / "proportions.csv")

# -- donor traits -------------------------------------------------------------------
traits = obs.groupby("donor_id", observed=True).agg(lambda s: s.iloc[0])
traits = traits[
    [
        "AD_status",
        "Parkinson_disease",
        "Vascular_status",
        "DLBD_status",
        "Schizophrenia",
        "ASCVD_status",
        "sex",
        "genetic_ancestry",
    ]
].copy()
age = obs.groupby("donor_id", observed=True)["development_stage"].first().astype(str)
traits["age"] = age.str.extract(r"(\d+)").astype(float).to_numpy()
traits.index.name = "participant_id"
traits.to_csv(out / "traits.csv")

# -- stratified cell subset -----------------------------------------------------------
keep = []
for _, group in obs.groupby("donor_id", observed=True):
    if len(group) <= CELLS_PER_DONOR:
        keep.append(group.index.to_numpy())
        continue
    frac = CELLS_PER_DONOR / len(group)
    for _, sub in group.groupby("subtype", observed=True):
        k = max(1, int(round(len(sub) * frac)))
        keep.append(
            rng.choice(sub.index.to_numpy(), size=min(k, len(sub)), replace=False)
        )
selected = np.sort(np.concatenate(keep))
subset = full[selected].to_memory()
subset.obs = subset.obs[["donor_id", "class", "subclass", "subtype"]]
print(f"Subset: {subset.n_obs} nuclei, {subset.obs.donor_id.nunique()} donors")

# -- class markers: mean log2 CPM in class minus mean outside, expressed genes only ---
logcpm = subset.copy()
sc.pp.normalize_total(logcpm, target_sum=1e6)
sc.pp.log1p(logcpm, base=2)
class_sums = sc.get.aggregate(logcpm, by="class", func="sum")
total_sum = np.asarray(logcpm.X.sum(axis=0)).ravel()
overall_mean = total_sum / logcpm.n_obs
n_per_class = logcpm.obs["class"].value_counts()

markers, used = {}, set()
for cls in pd.unique(logcpm.obs["class"]):
    cls_sum = np.asarray(class_sums[cls].layers["sum"]).ravel()
    n_cls = n_per_class[cls]
    cls_mean = cls_sum / n_cls
    specificity = cls_mean - (total_sum - cls_sum) / max(1, logcpm.n_obs - n_cls)
    expressed = cls_mean > 0.5
    markers[cls] = []
    for j in np.argsort(-np.where(expressed, specificity, -np.inf)):
        if not expressed[j] or len(markers[cls]) == N_MARKERS_PER_CLASS:
            break
        if gene_names[j] not in used:
            markers[cls].append(j)
            used.add(gene_names[j])
    print(
        f"  {cls:7s} {len(markers[cls])} markers, e.g. {list(gene_names[markers[cls][:4]])}"
    )

# -- gene panel: all markers + top expressed genes, one index per gene symbol ---------
panel = {}
for j in [j for idx in markers.values() for j in idx]:
    panel.setdefault(gene_names[j], j)
n_panel = len(panel) + N_TOP_EXPRESSED
for j in np.argsort(-overall_mean):
    if len(panel) == n_panel:
        break
    panel.setdefault(gene_names[j], j)
panel = np.sort(list(panel.values()))
print(f"Gene panel: {len(panel)} genes")

full.file.close()

# -- atlas.h5mu ---------------------------------------------------------------------------
atlas_obs = subset.obs.rename(
    columns={
        "donor_id": "participant_id",
        "class": "cell_class",
        "subtype": "subpopulation",
    }
)[["participant_id", "cell_class", "subclass", "subpopulation"]]
atlas_obs = atlas_obs.astype(str).astype("category")
atlas_obs.index = pd.Index(
    [f"cell_{i:06d}" for i in range(len(atlas_obs))], name="cell_id"
)
atlas = ad.AnnData(
    X=subset.X[:, panel],
    obs=atlas_obs,
    var=pd.DataFrame(index=pd.Index(gene_names[panel], name="gene_symbol")),
)
atlas.layers["counts"] = atlas.X.copy()
mu.MuData({"rna": atlas}).write_h5mu(out / "atlas.h5mu", compression="gzip")
print(f"Wrote atlas.h5mu  ({atlas.n_obs} nuclei x {atlas.n_vars} genes)")
