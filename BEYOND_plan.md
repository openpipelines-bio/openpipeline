# BEYOND Implementation Plan for OpenPipeline

Reference: [BEYOND_DLPFC](https://github.com/naomihabiblab/BEYOND_DLPFC) — Habib lab methodology for
uncovering trajectories of brain aging via coordinated cellular community dynamics in snRNA-seq data.

---

## Full Pipeline: FASTQ → BEYOND Outputs

```
┌─────────────────────────────────────────────────────────────────────────────────────────────┐
│  INPUT: per-library scRNA-seq FASTQ files + reference genome                                │
└────────────────────────────────┬────────────────────────────────────────────────────────────┘
                                 │
                    ┌────────────▼────────────┐
                    │  mapping/               │   [EXISTING]
                    │  cellranger_multi       │
                    └────────────┬────────────┘
                                 │  output folder (raw + filtered matrices)
                    ┌────────────▼────────────┐
                    │ convert/                │   [EXISTING]
                    │  from_cellranger_       │
                    │  multi_to_h5mu          │
                    │  (use raw matrix)       │
                    └────────────┬────────────┘
                                 │  raw h5mu per sample (unfiltered, for cellbender)
                    ┌────────────▼────────────┐
                    │ correction/             │   [EXISTING]
                    │  cellbender_remove_     │
                    │  background             │
                    └────────────┬────────────┘
                                 │  denoised h5mu
                    ┌────────────▼────────────┐
                    │ qc/calculate_qc_metrics │   [EXISTING]
                    └────────────┬────────────┘
                                 │  h5mu + QC obs columns
                    ┌────────────▼────────────┐
                    │ filter/filter_with_     │   [EXISTING]
                    │  counts                 │
                    │ filter/do_filter        │
                    └────────────┬────────────┘
                                 │  filtered h5mu
                    ┌────────────▼────────────┐
                    │ filter/filter_with_     │   [EXISTING]
                    │  scrublet               │   (doublet score + call)
                    └────────────┬────────────┘
                                 │  h5mu + doublet scores
                    ┌────────────▼────────────┐
                    │ dataflow/concatenate_   │   [EXISTING]
                    │  h5mu                   │   (all libraries → one object)
                    └────────────┬────────────┘
                                 │  multi-donor h5mu
                    ┌────────────▼────────────┐
                    │ transform/normalize_    │   [EXISTING]
                    │  total + log1p + scale  │
                    └────────────┬────────────┘
                                 │
                    ┌────────────▼────────────┐
                    │ feature_annotation/     │   [EXISTING]
                    │  highly_variable_       │
                    │  features_scanpy        │
                    └────────────┬────────────┘
                                 │
                    ┌────────────▼────────────┐
                    │ dimred/pca              │   [EXISTING]
                    └────────────┬────────────┘
                                 │
                    ┌────────────▼────────────┐
                    │ integrate/harmonypy     │   [EXISTING]
                    │                         │   (multi-donor batch correction)
                    └────────────┬────────────┘
                                 │  batch-corrected latent space (obsm["X_pca_integrated"])
                    ┌────────────▼────────────┐
                    │ neighbors/find_neighbors│   [EXISTING]
                    └────────────┬────────────┘
                                 │
                    ┌────────────▼────────────┐
                    │ cluster/leiden          │   [EXISTING]
                    │  (initial cell types)   │
                    └────────────┬────────────┘
                                 │
                    ┌────────────▼────────────┐
                    │ annotate/celltypist     │   [EXISTING]
                    │                         │   (broad cell-type labels)
                    └────────────┬────────────┘
                                 │  h5mu with obs["celltypist_pred"] → rename to obs["cell_type"]
                                 │
        ┌────────────────────────┴────────────────────────┐
        │          SPLIT by cell type (per-type objects)  │
        │          dataflow/split_h5mu                    │   [EXISTING]
        └──┬──────────────────────────────────────────────┘
           │  (one h5mu per cell type, processed in parallel)
           │
  ┌────────▼────────────────────────────────────────────────────────────┐
  │  PER-CELL-TYPE BLOCK (Excitatory neurons, Inhibitory, Astrocytes…)  │
  │                                                                     │
  │  ┌─────────────────────────────┐                                    │
  │  │ dimred/pca (cell-type PCA)  │  [EXISTING]                        │
  │  └──────────────┬──────────────┘                                    │
  │                 │                                                   │
  │  ┌──────────────▼──────────────┐                                    │
  │  │ neighbors/find_neighbors    │  [EXISTING]                        │
  │  └──────────────┬──────────────┘                                    │
  │                 │                                                   │
  │  ┌──────────────▼──────────────┐                                    │
  │  │ cluster/leiden              │  [EXISTING]  (subpopulations)      │
  │  └──────────────┬──────────────┘                                    │
  │                 │                                                   │
  │  ┌──────────────▼──────────────┐                                    │
  │  │ filter/do_filter            │  [EXISTING]  (remove low-quality   │
  │  │ (remove doublets + low-QC)  │               subclusters)         │
  │  └──────────────┬──────────────┘                                    │
  │                 │                                                   │
  │  ┌──────────────▼──────────────┐                                    │
  │  │ dimred/umap                 │  [EXISTING]                        │
  │  └──────────────┬──────────────┘                                    │
  │                 │                                                   │
  │  ┌──────────────▼──────────────┐                                    │
  │  │ differential_expression/    │  [EXISTING]                        │
  │  │  create_pseudobulk          │                                    │
  │  └──────────────┬──────────────┘                                    │
  │                 │                                                   │
  │  ┌──────────────▼──────────────┐                                    │
  │  │ differential_expression/    │  [EXISTING]                        │
  │  │  deseq2                     │  (subpopulation marker genes)      │
  │  └──────────────┬──────────────┘                                    │
  │                 │                                                   │
  │  ┌──────────────▼──────────────┐                                    │
  │  │ interpret/pathway_enrichment│  [NEW #7]  GSEA/ORA per            │
  │  │                             │            subpopulation           │
  │  └──────────────┬──────────────┘                                    │
  └─────────────────┴───────────────────────────────────────────────────┘
                    │  annotated subpopulations (all cell types merged)
                    │
        ┌───────────▼─────────────────────────┐
        │ dataflow/concatenate_h5mu           │  [EXISTING]
        │  (merge cell-type objects → atlas)  │
        └───────────┬─────────────────────────┘
                    │  full atlas h5mu with obs["subpopulation"]
                    │
        ┌───────────▼─────────────────────────┐
        │ metadata/calculate_proportions      │  [NEW #1]
        │  obs["participant_id"] ×            │
        │  obs["subpopulation"] → proportion  │
        │  matrix stored in uns/obsm          │
        └───────────┬─────────────────────────┘
                    │  h5mu with participant × subpopulation proportions
                    │
        ┌───────────▼─────────────────────────┐
        │ dimred/phate                        │  [NEW #2]
        │  input: proportion matrix           │
        │  output: obsm["X_phate"]            │
        │  → cellular landscape               │
        └───────────┬─────────────────────────┘
                    │  h5mu with PHATE landscape
                    │
                    ├────────────────────────────────────┐
                    │                                    │
        ┌───────────▼─────────────────────────┐ ┌────────▼───────────────────────────┐
        │ trajectory/palantir                 │ │ trajectory/via                     │
        │  input: X_phate + start cell        │ │  input: X_phate + root cluster     │
        │  output:                            │ │  output:                           │
        │    obs["palantir_pseudotime"]       │ │    obs["via_pseudotime"]           │
        │    obsm["fate_probabilities"]       │ │    obsm["via_graph"]               │
        │  [NEW #3]                           │ │  [NEW #4]  (validation/alt)        │
        └───────────┬─────────────────────────┘ └────────────────────────────────────┘
                    │  pseudotime + fate probabilities
                    │
        ┌───────────▼─────────────────────────┐
        │ trajectory/pseudotime_dynamics      │  [NEW #5]
        │  input: proportion matrix +         │
        │         pseudotime                  │
        │  output: GAM fits per subpopulation │
        │    uns["dynamics"]["peak_time"]     │
        │    uns["dynamics"]["gam_params"]    │
        └───────────┬─────────────────────────┘
                    │  dynamics curves per subpopulation
                    │
        ┌───────────▼─────────────────────────┐
        │ cluster/cellular_communities        │  [NEW #6]
        │  input: proportion matrix +         │
        │         dynamics curves             │
        │  method:                            │
        │    1. co-occurrence (Pearson corr   │
        │       of proportion vectors)        │
        │    2. dynamics similarity (DTW /    │
        │       corr of GAM curves)           │
        │    3. combined similarity matrix    │
        │    4. hierarchical / spectral       │
        │       clustering → communities      │
        │  output:                            │
        │    uns["cellular_communities"]      │
        │    obs["community_id"]              │
        └───────────┬─────────────────────────┘
                    │  cellular communities
                    │
        ┌───────────▼─────────────────────────┐
        │ stats/trait_associations            │  [NEW #8]
        │  input: proportion matrix +         │
        │         community assignments +     │
        │         clinical/trait table (CSV)  │
        │  method: linear mixed models        │
        │  output:                            │
        │    uns["trait_associations"]        │
        │    (subpopulation/community ×       │
        │     trait: β, p, FDR)               │
        └───────────┬─────────────────────────┘
                    │
┌───────────────────▼───────────────────────────────────────────────────┐
│  OUTPUTS                                                              │
│                                                                       │
│  atlas.h5mu                 Full annotated single-cell atlas          │
│    obs["cell_type"]         Broad cell type                           │
│    obs["subpopulation"]     Fine-grained subcluster ID                │
│    obs["community_id"]      BEYOND cellular community                 │
│    obs["palantir_pseudotime"] Disease trajectory pseudotime           │
│    obsm["X_phate"]          Cellular landscape coordinates            │
│    obsm["fate_probabilities"] Branch fate probabilities               │
│    uns["dynamics"]          GAM-fitted proportion dynamics            │
│    uns["cellular_communities"] Community membership + metadata        │
│    uns["trait_associations"]  Per-subpop/community trait stats        │
│    uns["pathway_enrichment"]  Per-subpop GSEA/ORA results             │
│                                                                       │
│  proportions.h5mu           Participant × subpopulation matrix        │
│  de_results/                Per-subpopulation DESeq2 tables           │
│  pathway_results/           Per-subpopulation enrichment tables       │
│  trait_associations.csv     Community × trait association table       │
└───────────────────────────────────────────────────────────────────────┘
```

---

## Gap Analysis

### Already available (existing components)

| BEYOND step | OpenPipeline component |
|---|---|
| Read alignment + count matrix | `mapping/cellranger_multi` |
| CellRanger multi → h5mu conversion | `convert/from_cellranger_multi_to_h5mu` |
| Ambient RNA removal | `correction/cellbender_remove_background` |
| snRNA-seq QC metrics | `qc/calculate_qc_metrics` |
| Cell/feature filtering | `filter/filter_with_counts`, `filter/do_filter` |
| Doublet detection | `filter/filter_with_scrublet` |
| Multi-library concatenation | `dataflow/concatenate_h5mu` |
| Normalization + log-transform | `transform/normalize_total`, `transform/log1p`, `transform/scale` |
| Highly variable features | `feature_annotation/highly_variable_features_scanpy` |
| PCA | `dimred/pca` |
| Batch correction (multi-donor) | `integrate/harmonypy` |
| Neighbor graph | `neighbors/find_neighbors` |
| Leiden clustering | `cluster/leiden` |
| Broad cell type annotation | `annotate/celltypist` |
| UMAP visualization | `dimred/umap` |
| Split by cell type | `dataflow/split_h5mu` |
| Pseudobulk aggregation | `differential_expression/create_pseudobulk` |
| Differential expression | `differential_expression/deseq2` |
| Format conversion | `convert/` namespace |
| Label transfer (replication cohort) | `labels_transfer/knn`, `labels_transfer/xgboost` |

### Partially available (alternative exists, not exact match)

| BEYOND step | Gap | Closest existing |
|---|---|---|
| Cellular community biology | BEYOND uses co-occurrence + dynamics; not ligand-receptor | `interpret/lianapy` |
| RNA velocity context | scVelo present; BEYOND uses Palantir/VIA for pseudotime | `velocity/scvelo` |

### Component compatibility notes

| Component | Issue | Resolution |
|---|---|---|
| `integrate/harmony` | `status: disabled` in config | Use `integrate/harmonypy` (Python port, identical I/O: obsm["X_pca"] → obsm["X_pca_integrated"]; `--obs_covariates` is required) |
| `mapping/cellranger_multi` | Outputs a folder, not h5mu | Follow with `convert/from_cellranger_multi_to_h5mu` using the raw matrix; cellbender then runs on the resulting raw h5mu |
| `annotate/celltypist` | Output stored as `obs["celltypist_pred"]`, not `obs["cell_type"]` | Add `metadata/move_obsm_to_obs` or direct rename step before `dataflow/split_h5mu`; also requires a pretrained `--model` (.pkl) or a `--reference` h5mu with `--reference_obs_target` |

### Not available — must be built (8 new components)

| # | Component | Namespace | Method | Language |
|---|---|---|---|---|
| 1 | `calculate_proportions` | `metadata/` | Group-by aggregation (participant × subpopulation) | Python |
| 2 | `phate` | `dimred/` | PHATE embedding (KrishnaswamyLab/phate) | Python |
| 3 | `palantir` | `trajectory/` | Pseudotime + fate probabilities (palantir package) | Python |
| 4 | `via` | `trajectory/` | Trajectory inference (pyVIA package) | Python | **BLOCKED** — pyVIA 0.2.4 has hard incompatibilities with scipy>=1.12 (required by anndata) and numpy>=2. See [blocking issues](#trajectory-via-blocked). |
| 5 | `pseudotime_dynamics` | `trajectory/` | GAM fits of proportions along pseudotime (pygam / mgcv) | Python/R |
| 6 | `cellular_communities` | `cluster/` | Co-occurrence + dynamics similarity → community clustering | Python |
| 7 | `pathway_enrichment` | `interpret/` | GSEA / ORA (gseapy, MSigDB/GO/KEGG) | Python |
| 8 | `trait_associations` | `stats/` *(new namespace)* | Linear mixed models per subpop × clinical trait | R |

---

## Implementation Plan

### Phase 1 — Proportion matrix (prerequisite for everything BEYOND)

**`metadata/calculate_proportions`** — ~1 day

- Input: h5mu with `obs["participant_id"]` and `obs["subpopulation"]`
- Output: h5mu with participant × subpopulation proportion matrix in `uns["proportions"]`
  and a condensed AnnData in `obsm["proportions"]`
- Simple pandas group-by; model after existing `metadata/` components

---

### Phase 2 — PHATE landscape (`dimred/phate`)  — ~1 day

- Input: h5mu, `--obsm_input` key (e.g. `X_pca` or `proportions`)
- Output: h5mu with `obsm["X_phate"]`
- Docker: Python image + `pip install phate scprep`
- Mirror structure of `dimred/umap/config.vsh.yaml` and `dimred/umap/script.py`
- Key params: `n_components`, `knn`, `decay`, `t` (auto or fixed)

---

### Phase 3 — Trajectory inference

**`trajectory/palantir`** — ~3 days

- Input: h5mu, `--obsm_key` (X_phate recommended), `--start_cell` or `--start_cell_cluster`
- Output: h5mu with `obs["palantir_pseudotime"]`, `obsm["palantir_fate_probabilities"]`,
  `uns["palantir_waypoints"]`
- Package: `palantir` (pip)
- Key params: `num_waypoints`, `terminal_states` (obs column or explicit cell IDs),
  `n_jobs`, `scale_components`
- Note: use `reticulate`-free pure Python; write output back with anndata

**`trajectory/via`** — ~2 days

- Input: h5mu, `--obsm_key`, `--root_user` (cluster label or cell ID)
- Output: h5mu with `obs["via_pseudotime"]`, `uns["via_graph"]`
- Package: `pyVIA` (pip)
- Serve as validation/alternative to Palantir; both can run in parallel in workflow

**`trajectory/pseudotime_dynamics`** — ~3 days

- Input: h5mu with `obs["palantir_pseudotime"]` and `uns["proportions"]`
- Output: h5mu with `uns["dynamics"]`:
  - per-subpopulation GAM curve (smoothed proportion vs pseudotime)
  - `peak_pseudotime` per subpopulation
  - fit statistics (R², p-value)
- Method: `pygam` (Python) — GeneralizedAdditiveModel with n_splines, lam
- Key params: `--pseudotime_key`, `--n_splines`, `--n_pseudotime_bins`, `--min_cells`

---

### Phase 4 — BEYOND core (`cluster/cellular_communities`) — ~5 days

The novel BEYOND method itself:

1. **Co-occurrence matrix**: Pearson correlation of participant proportion vectors
   across subpopulations → subpopulation × subpopulation similarity
2. **Dynamics similarity matrix**: correlation (or DTW) of GAM-fitted dynamics curves
   → subpopulation × subpopulation temporal similarity
3. **Combined similarity**: weighted sum of both matrices (tunable `--alpha`)
4. **Community detection**: hierarchical clustering (Ward linkage) or spectral clustering
   on combined similarity → community labels

- Input: h5mu with `uns["proportions"]` + `uns["dynamics"]`
- Output: h5mu with `uns["cellular_communities"]` (subpopulation → community mapping),
  `obs["community_id"]`
- Key params: `--n_communities` (or auto via silhouette), `--alpha` (co-occurrence weight),
  `--method` (hierarchical / spectral)

---

### Phase 5 — Supporting analyses

**`interpret/pathway_enrichment`** — ~2 days

- Input: DE result table (CSV from deseq2), gene set GMT file or database name
- Output: h5mu `uns["pathway_enrichment"]` + CSV result table
- Package: `gseapy` (supports Enrichr, prerank GSEA, ORA)
- Key params: `--gene_set` (MSigDB_Hallmark, GO_Biological_Process, KEGG…),
  `--method` (gsea / enrichr / ora), `--fc_column`, `--pval_column`

**`stats/trait_associations`** — ~3 days  *(new namespace `stats/`)*

- Input: h5mu with proportion matrix + community assignments, clinical trait CSV
  (participant_id + trait columns)
- Output: association table CSV + `uns["trait_associations"]`
- Method: `lme4` linear mixed model in R (random effect = cohort/batch)
- Key params: `--trait_columns`, `--covariates`, `--random_effects`, `--fdr_method`

---

### Phase 6 — Nextflow workflows — ~3 days

New workflow directory: `src/workflows/beyond/`

```
src/workflows/beyond/
  config.vsh.yaml
  1_preprocessing.nf          cellranger_multi → from_cellranger_multi_to_h5mu → cellbender → qc → filter → scrublet
  2_cell_type_atlas.nf        normalize → HVF → pca → harmonypy → neighbors → leiden
                              → celltypist (+ rename obs["celltypist_pred"] → obs["cell_type"]) → split
                              → per-type leiden → concatenate → calculate_proportions
  3_landscape.nf              phate
  4_trajectories.nf           palantir → pseudotime_dynamics  (+ via in parallel for validation)
  5_communities.nf            cellular_communities
  6_associations.nf           pathway_enrichment + trait_associations
  beyond.nf                   master workflow chaining 1→2→3→4→5→6
```

---

## Build Order & Effort

| # | Component | Effort | Depends on | Priority |
|---|---|---|---|---|
| 1 | `metadata/calculate_proportions` | 1 day | atlas h5mu | Critical |
| 2 | `dimred/phate` | 1 day | proportions | Critical |
| 3 | `trajectory/palantir` | 3 days | phate | Critical |
| 5 | `trajectory/pseudotime_dynamics` | 3 days | palantir | Critical |
| 6 | `cluster/cellular_communities` | 5 days | pseudotime_dynamics | Critical |
| 7 | `interpret/pathway_enrichment` | 2 days | deseq2 output | High |
| 4 | `trajectory/via` | 2 days | phate | ~~Medium~~ **BLOCKED** |
| 8 | `stats/trait_associations` | 3 days | communities | Medium |
| — | Workflows | 3 days | all components | Last |
| **Total** | | **~23 days** | | |

Components 1–6 are the **critical path**. Components 7, 8, and workflows can be
parallelized once the critical path is unblocked. Component 4 is blocked (see below).

---

## trajectory/via — BLOCKED {#trajectory-via-blocked}

**Status**: Skipped — pyVIA 0.2.4 is incompatible with the anndata base stack.

| Issue | Detail |
|---|---|
| `numpy.bool8` removed in NumPy 2.x | `nptyping` (pyVIA dep) uses `np.bool8`; requires `numpy<2` |
| `scipy>=1.12` required by anndata | pyVIA's `csr_matrix(zip(*edges))` breaks with newer scipy API when edge lists are empty after Jaccard pruning |
| Markov NaN crash | Even after patching `get_sparse_from_igraph`, isolated nodes in pruned graph → zero row-sum → NaN transition probabilities → `_simulate_markov` crashes |

**Resolution options (future)**:
1. Check whether the original BEYOND_DLPFC R code uses pyVIA or an R-native tool (slingshot / monocle3); if R-native, rewrite as `r_script`
2. Wait for a pyVIA release that supports scipy>=1.12
3. Pin `scipy<1.11` in a separate Docker image (breaks anndata base stack sharing)

Palantir (`trajectory/palantir`) fully covers the pseudotime requirement for the BEYOND critical path.

---

#   Architecture to run the "full" pipeline

The right split:

1. ingestion/cellranger_multi (EXISTS) — FASTQs → per-sample h5mu
2. ingestion/cellranger_postprocessing (EXISTS) — cellbender + initial filter
3. beyond/atlas_building (NEW, done ↑) — per-sample h5mu → atlas with obs["subpopulation"]
4. beyond/trajectory_analysis (NEW, to build) — atlas → BEYOND outputs

Key reasons: atlas_building is useful for any snRNA-seq project (not BEYOND-specific); trajectory_analysis can be applied to any existing atlas; Nextflow -resume works at sub-workflow boundaries; cellranger version bumps don't touch trajectory code.

