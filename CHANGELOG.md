# openpipeline 0.3.1

## NEW FUNCTIONALITY

* `bin/port_from_czbiohub_utilities.sh`: Added helper script to import components and pipelines from `czbiohub/utilities`

Imported components from `czbiohub/utilities`:

* `demux/cellranger_mkfastq`: Demultiplex raw sequencing data.

* `mapping/cellranger_count`: Align fastq files using Cell Ranger count.

* `mapping/cellranger_count_split`: Split 10x Cell Ranger output directory into separate output fields.

Imported workflows from `czbiohub/utilities`:

* `workflows/1_ingestion/cellranger`: Use Cell Ranger to preprocess 10x data.

* `workflows/1_ingestion/cellranger_demux`: Use cellranger demux to demultiplex sequencing BCL output to FASTQ.

* `workflows/1_ingestion/cellranger_mapping`: Use cellranger count to align 10x fastq files to a reference.


## MINOR CHANGES

* Fix `interactive/run_cirrocumulus` script raising `NotImplementedError` caused by using `MutData.var_names_make_unique()` 
on each modality instead of on the whole `MuData` object.

* Fix `transform/normalize_total` and `interactive/run_cirrocumulus` component build missing a hdf5 dependency.

* `interactive/run_cellxgene`: Updated container to ubuntu:focal because it contains python3.6 but cellxgene dropped python3.6 support.

* `download/download_file`: Switch base container from `ubuntu:22.04` to `bash:5.1` because `ubuntu:22.04`
  was causing some edge case issues.

* `mapping/bd_rhapsody_wta`: Set `--parallel` to true by default.

* `mapping/bd_rhapsody_wta`: Translate Bash script into Python.

* `download/sync_test_resources`: Add `--dryrun` and `--quiet` arguments.

## BUG FIXES

* `mapping/bd_rhapsody`: Use a smaller test dataset to reduce test time and make sure that the Github Action runners do not run out of disk space.

* `download/sync_test_resources`: Disable the use of the Amazon EC2 instance metadata service to make script work on Github Actions runners.

* `convert/from_h5mu_to_seurat`: Fix unit test requiring Seurat by using native R functions to test the Seurat object instead.

# openpipeline 0.3.0

* Add `tx_processing` pipeline with following components:
  - `filter_with_counts`
  - `filter_with_scrublet`
  - `filter_with_hvg`
  - `do_filter`
  - `normalize_total`
  - `regress_out`
  - `log1p`
  - `pca`
  - `find_neighbors`
  - `leiden`
  - `umap`

# openpipeline 0.2.0

## NEW FUNCTIONALITY

* Added `from_10x_to_h5ad` and `download_10x_dataset` components.

## MINOR CHANGES
* Workflow `bd_rhapsody_wta`: Minor change to workflow to allow for easy processing of multiple samples with a tsv.

* Component `bd_rhapsody_wta`: Added more parameters, `--parallel` and `--timestamps`.

* Added `pbmc_1k_protein_v3` as a test resource.

* Translate `bd_rhapsody_extracth5ad` from R into Python script.

* `bd_rhapsody_wta`: Remove temporary directory after execution.


# openpipeline 0.1.0

* Initial release containing only a `bd_rhapsody_wta` pipeline and corresponding components.
