# openpipeline 0.5.0

## NEW FUNCTIONALITY

* `graph/bbknn`: BBKNN network generation.

* `workflows/integration/multiomics`: implement pipeline for processing multiple multiomics samples.

* `transform/scaling`: Scale data to unit variance and zero mean.

* `integrate/harmony` and `integrate/harmonypy`: Run a Harmony integration analysis (R-based and Python-based, respectively).

* `integrate/concat`: Added 'move' option to `--other_axis_mode`, which allows merging .obs and .var by only keeping elements of the matrices which are the same in each of the samples, moving the conflicting values to .varm or .obsm.

* `reference/make_reference`: Download a transcriptomics reference and preprocess it (adding ERCC spikeins and filtering with a regex).

* `reference/build_bdrhap_reference`: Compile a reference into a STAR index in the format expected by BD Rhapsody.

## MAJOR CHANGES

* Multiple components: update to anndata 0.8 with mudata 0.2.0. This means that the format of the .h5mu files have changed.

* `process_rna/singlesample`: move transformation counts into layers instead of overwriting .X.

## MINOR CHANGES

* `velocity/velocyto`: Allow configuring memory and parallellisation.

## BUG FIXES

* Several components: revert using slim versions of containers because they do not provide the tools to run nextflow with trace capabilities.

* `integrate/concat`: Fix an issue where joining boolean values caused `TypeError`.

* `workflows/process_rna/multisample`, `workflows/process_rna/singlesample` and `workflows/integration/multimodal_integration`: use nextflow trace reporting when running integration tests.

## DEPRECATED

# openpipeline 0.4.1

## BUG FIXES

* `workflows/ingestion/bd_rhapsody_wta`: use ':' as a seperator for multiple input files and fix integration tests.

## MINOR CHANGES

* Several components: pin mudata and scanpy dependencies so that anndata version <0.8.0 is used.

# openpipeline 0.4.0

## NEW FUNCTIONALITY

* `convert/from_bdrhap_to_h5mu`: Merge one or more BD rhapsody outputs into an h5mu file.

* `split/split_modalities`: Split the modalities from a single .h5mu multimodal sample into seperate .h5mu files. 

* `integrate/concat`: Combine data from multiple samples together.

## MINOR CHANGES

* `mapping/bd_rhapsody_wta`: Update to BD Rhapsody 1.10.1.

* `mapping/bd_rhapsody_wta`: Add parameters for overriding the minimum RAM & cores. Add `--dryrun` parameter.

* Switch to Viash 0.5.14.

* `convert/from_bdrhap_to_h5mu`: Update to BD Rhapsody 1.10.1.

* `resources_test/bdrhap_5kjrt`: Add subsampled BD rhapsody datasets to test pipeline with.

* `resources_test/bdrhap_ref_gencodev40_chr1`: Add subsampled reference to test BD rhapsody pipeline with.

* `integrate/merge`: Merge several unimodal .h5mu files into one multimodal .h5mu file.

* Updated several python docker images to slim version.

* `mapping/cellranger_count_split`: update container from ubuntu focal to ubuntu jammy

* `download/sync_test_resources`: update AWS cli tools from 2.7.11 to 2.7.12 by updating docker image

* `download/download_file`: now uses bash container instead of python.

* `mapping/bd_rhapsody_wta`: Use squashed docker image in which log4j issues are resolved.

## BUG FIXES

* `workflows/utils/WorkflowHelper.nf`: Renamed `utils.nf` to `WorkflowHelper.nf`.

* `workflows/utils/WorkflowHelper.nf`: Fix error message when required parameter is not specified.

* `workflows/utils/WorkflowHelper.nf`: Added helper functions:
  - `readConfig`: Read a Viash config from a yaml file.
  - `viashChannel`: Create a channel from the Viash config and the params object.
  - `helpMessage`: Print a help message and exit.

* `mapping/bd_rhapsody_wta`: Update picard to 2.27.3.

## DEPRECATED

* `convert/from_bdrhap_to_h5ad`: Deprecated in favour for `convert/from_bdrhap_to_h5mu`.

* `convert/from_10xh5_to_h5ad`: Deprecated in favour for `convert/from_10xh5_to_h5mu`.

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

* `mapping/bd_rhapsody_wta`: Set `--parallel` to true by default.

* `mapping/bd_rhapsody_wta`: Translate Bash script into Python.

* `download/sync_test_resources`: Add `--dryrun`, `--quiet`, and `--delete` arguments.

* `convert/from_h5mu_to_seurat`: Use `eddelbuettel/r2u:22.04` docker container in order to speed up builds by downloading precompiled R packages.

* `mapping/cellranger_count`: Use 5Gb for testing (to adhere to github CI runner memory constraints).

* `convert/from_bdrhap_to_h5ad`: change test data to output from `mapping/bd_rhapsody_wta` after reducing the BD Rhapsody test data size.

* Various `config.vsh.yaml`s: Renamed `values:` to `choices:`.

* `download/download_file` and `transfer/publish`: Switch base container from `bash:5.1` to `python:3.10`.

* `mapping/bd_rhapsody_wta`: Make sure procps is installed.

## BUG FIXES

* `mapping/bd_rhapsody_wta`: Use a smaller test dataset to reduce test time and make sure that the Github Action runners do not run out of disk space.

* `download/sync_test_resources`: Disable the use of the Amazon EC2 instance metadata service to make script work on Github Actions runners.

* `convert/from_h5mu_to_seurat`: Fix unit test requiring Seurat by using native R functions to test the Seurat object instead.

* `mapping/cellranger_count` and `bcl_demus/cellranger_mkfastq`: cellranger uses the `--parameter=value` formatting instead of `--parameter value` to set command line arguments.

* `mapping/cellranger_count`: `--nosecondary` is no longer always applied.

* `mapping/bd_rhapsody_wta`: Added workaround for bug in Viash 0.5.12 where triple single quotes are incorrectly escaped (viash-io/viash#139).

## DEPRECATED

* `bcl_demux/cellranger_mkfastq`: Duplicate of `demux/cellranger_mkfastq`.

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
