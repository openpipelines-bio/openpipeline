# OpenPipeline — Repository Summary

## Project Overview

**Name:** OpenPipeline (`openpipeline`)
**Organization:** openpipelines-bio
**Version:** 4.0.3 (4.0.4-dev in progress)
**License:** MIT

Extensible single-cell analysis pipelines for reproducible and large-scale single-cell processing built on **Viash 0.9.x** and **Nextflow DSL2**. Covers the full analysis lifecycle: demultiplexing → read mapping → QC → normalization → dimensionality reduction → clustering → integration → annotation.

**Links:**
- Website: https://openpipelines.bio
- Docs: https://openpipelines.bio/fundamentals
- Viash Hub: https://www.viash-hub.com/packages/openpipeline
- Docker Registry: ghcr.io/openpipelines-bio/openpipeline

---

## Directory Structure

```
openpipeline/
├── src/                            # All components and workflows (33 namespaces)
│   └── workflows/                  # Nextflow workflow definitions
├── target/                         # (Generated) Compiled Nextflow modules
├── resources_test_scripts/         # Test data download scripts (22 scripts, S3-backed)
├── .github/workflows/              # GitHub Actions CI/CD (5 workflows)
├── _viash.yaml                     # Root Viash config (v0.9.4, ghcr.io registry)
├── nextflow.config                 # Nextflow template configuration
├── main.nf                         # Placeholder Nextflow entry point
├── CHANGELOG.md                    # Full version history (~91 KB)
├── tasks.md                        # Developer workflow guide
└── ruff.toml / .pylintrc / .lintr  # Code quality configuration
```

---

## Namespaces and Components

33 namespaces with ~120 components total.

### Data Ingestion & Mapping

| Namespace | Key Components | Purpose |
|-----------|----------------|---------|
| `mapping` | cellranger_count, cellranger_multi, star_align, htseq_count, bd_rhapsody (12 total) | Read alignment and count matrix generation |
| `demux` | bcl2fastq, bcl_convert, cellranger_mkfastq (4 total) | BCL → FASTQ demultiplexing |
| `genetic_demux` | souporcell, vireo, demuxlet, freemuxlet, bcftools, cellsnp (10 total) | Genetic demultiplexing |
| `reference` | build_cellranger_reference, build_star_reference, make_reference (6 total) | Reference genome construction |
| `correction` | cellbender_remove_background | Ambient RNA correction |

### Format Conversion

| Namespace | Key Components | Purpose |
|-----------|----------------|---------|
| `convert` | from_10xh5_to_h5mu, from_10xmtx_to_h5mu, from_bdrhap_to_h5mu, from_h5mu_to_h5ad, from/to_seurat, from/to_tiledb (13 total) | Inter-format conversion (h5mu, h5ad, Seurat, TileDB) |
| `dataflow` | concatenate_h5mu, split_h5mu, split_modalities, merge | Data merging and splitting |
| `compression` | compress_h5mu, tar_extract | File compression |
| `tiledb` | move_mudata_obs/obsm/obsp_to_tiledb | TileDB-SOMA integration |

### QC & Filtering

| Namespace | Key Components | Purpose |
|-----------|----------------|---------|
| `qc` | calculate_qc_metrics, calculate_atac_qc_metrics, fastqc, multiqc | QC metric calculation |
| `filter` | filter_with_counts, filter_with_scrublet, do_filter, subset_h5mu, intersect_obs (11 total) | Cell/feature filtering and doublet detection |
| `process_10xh5` | filter_10xh5 | 10X HDF5 processing |

### Normalization & Transformation

| Namespace | Key Components | Purpose |
|-----------|----------------|---------|
| `transform` | normalize_total, log1p, scale, clr, regress_out, tfidf, bpcells_regress_out (9 total) | Normalization and scaling |
| `feature_annotation` | highly_variable_features_scanpy, score_genes_scanpy, score_genes_cell_cycle_scanpy, align_query_reference | Feature selection and gene scoring |

### Dimensionality Reduction & Clustering

| Namespace | Components | Purpose |
|-----------|-----------|---------|
| `dimred` | pca, lsi, umap, tsne, densmap | Dimensionality reduction |
| `neighbors` | find_neighbors, bbknn | Neighbor graph |
| `cluster` | leiden | Graph-based clustering |

### Integration & Batch Correction

| Namespace | Components | Purpose |
|-----------|-----------|---------|
| `integrate` | harmony, harmonypy, scanorama, scvi, scarches, totalvi, totalvi_scarches | Batch correction methods |

### Annotation

| Namespace | Components | Purpose |
|-----------|-----------|---------|
| `annotate` | celltypist, singler, popv, scanvi, onclass, random_forest_annotation, svm_annotation | Cell type annotation |
| `labels_transfer` | knn, xgboost, api | Label transfer |

### Downstream Analysis

| Namespace | Components | Purpose |
|-----------|-----------|---------|
| `velocity` | scvelo, velocyto, velocyto_to_h5mu | RNA velocity |
| `differential_expression` | create_pseudobulk, deseq2 | Pseudobulk DE analysis |
| `interpret` | lianapy | Cell-cell interaction analysis |
| `query` | cellxgene_census, tiledb_soma_healthcheck | External data querying |

### Supporting

| Namespace | Purpose |
|-----------|---------|
| `metadata` | add_id, join_csv, grep_annotation_column, move_obsm_to_obs |
| `download` | download_file, sync_test_resources |
| `report` | mermaid (workflow diagrams) |
| `base` | Shared base requirements |
| `utils` | Shared test helpers and utilities |
| `authors` | Author profile definitions |

---

## Workflows (Nextflow)

31 workflows in `src/workflows/`, organized by analysis stage:

| Category | Workflows |
|----------|-----------|
| **Ingestion** | demux, cellranger_mapping, cellranger_multi, cellranger_postprocessing, bd_rhapsody, conversion, make_reference |
| **Single-sample** | rna_singlesample, prot_singlesample, gdo_singlesample |
| **Multi-sample** | rna_multisample, prot_multisample, log_normalize |
| **Multiomics** | process_samples, process_batches, neighbors_leiden_umap, dimensionality_reduction, split_h5mu, split_modalities |
| **Integration** | harmony_leiden, scvi_leiden, totalvi_leiden, bbknn_leiden, scanorama_leiden, totalvi_scarches_leiden |
| **Annotation** | celltypist, scanvi_scarches, harmony_knn, scvi_knn |
| **QC** | qc |
| **DE** | pseudobulk_deseq2 |

---

## Data Formats

**Primary format:** MuData (H5MU — HDF5-based multi-modal)

| Format | Use |
|--------|-----|
| H5MU | Primary internal format (multi-modal) |
| H5AD | Single-modality AnnData (legacy/interop) |
| Seurat RDS | R interoperability |
| TileDB-SOMA | Columnar storage for large datasets |

**Input sources:** 10X Genomics (GEX/VDJ/ATAC/ADT), BD Rhapsody, FASTQ, BCL, BAM/SAM

**Modalities supported:** RNA, ATAC, Antibody Capture (ADT/Protein), VDJ, AIRR

---

## Testing

Each component contains:
- `config.vsh.yaml` — declares `test_resources`
- `test.sh` or `test.py` — test script
- `test_data/test_data_script.sh` — generates test inputs

Test data is stored on S3 (`s3://openpipelines-data`) and downloaded via scripts in `resources_test_scripts/`.

**Commands:**
```bash
viash test src/<ns>/<comp>/config.vsh.yaml          # single component
viash ns test --parallel -q <namespace>             # entire namespace
viash ns build --setup cb                           # build all
```

---

## CI/CD (GitHub Actions)

| Workflow | Trigger | Purpose |
|----------|---------|---------|
| `viash-test.yml` | PR / push to main | Lint (ruff, styler, lintr), pre-commit checks |
| `main-build.yml` | Workflow dispatch | Build components, push Docker images |
| `release-build.yml` | Release creation | Tag and deploy release artifacts |
| `integration-test.yml` | Workflow dispatch | Cross-dependency integration tests |
| `create-documentation-pr.yml` | Schedule / dispatch | Generate documentation PRs |

Docker images pushed to `ghcr.io/openpipelines-bio/openpipeline`, using BioContainers base images (quay.io/biocontainers). Each image tag includes the build string, e.g., `tool:2.0.3--h5b5514e_1`.

---

## Resource Labels

Tiered resource allocation via Nextflow labels (`src/workflows/utils/labels.config`):

- `lowmem` / `midmem` / `highmem` / `veryhighmem`
- `lowcpu` / `midcpu` / `highcpu` / `veryhighcpu`
- `lowdisk` / `middisk` / `highdisk` / `veryhighdisk`

---

## Languages

| Language | Role |
|----------|------|
| Python | Primary implementation (scanpy, scvi-tools, numpy, scipy) |
| R | Supplementary analysis (Seurat, DESeq2, SingleR) |
| Bash/Shell | Wrapper scripts, pipeline orchestration |
| Nextflow (DSL2) | Workflow definitions |
| Groovy | Nextflow config files |

---

## Key Viash Conventions

- Viash version: **0.9.4**
- Use `engines:` / `runners:` (not deprecated `platforms:`)
- Use `argument_groups:` to organize arguments; `boolean_true` for flags
- Script rules: 2-space indent, `set -eo pipefail` after `## VIASH END`
- Unset boolean params before use: `[[ "$par_flag" == "false" ]] && unset par_flag`
- Use `meta_cpus` / `meta_memory_gb` (never `par_cpus` / `par_threads`)
- Build CLI args as arrays: `cmd_args=(...) && tool "${cmd_args[@]}"`
- Always put `--` between config path and component arguments in `viash run`

---

## Recent Activity (as of 2026-03-16)

- Cell Ranger 10 / AIRR rearrangement support
- `calculate_qc_metrics` rework (memory optimization)
- Highly variable features masking
- `split_h5mu` obsp subset fix
- docker/login-action bump (3 → 4)
