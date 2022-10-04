# openpipeline 0.5.1

## BREAKING CHANGES

* `reference/make_reference`: Input files changed from `type: string` to `type: file` to allow Nextflow to cache the input files fetched from URL.

* several components (except `from_h5ad_to_5hmu`): the `--modality` arguments no longer accept multiple values.

* Remove outdated `resources_test_scripts`.

* `convert/from_h5mu_to_seurat`: Disabled because MuDataSeurat is currently broken, see [https://github.com/PMBio/MuDataSeurat/issues/9](PMBio/MuDataSeurat#9).

* `integrate/concat`: Renamed --sample_names to --input_id and moved the ability to add sample id and to join the sample ids with the observation names to `metadata/add_id`

## NEW FUNCTIONALITY

* `metadata/add_id`: Add an id to a column in .obs. Also allows joining the id to the .obs_names.

* `workflows/ingestion/make_reference`: A generic component to build a transcriptomics reference into one of many formats.

* `integrate/add_metadata`: Add a csv containing metadata to the .obs or .var field of a mudata file.

* `DataFlowHelper.nf`: Added `passthroughMap`. Usage:

  ```groovy
  include { passthroughMap as pmap } from "./DataFlowHelper.nf"
  
  workflow {
    Channel.fromList([["id", [input: "foo"], "passthrough"]])
      | pmap{ id, data ->
        [id, data + [arg: 10]]
      }
  }
  ```
  Note that in the example above, using a regular `map` would result in an exception being thrown,
  that is, "Invalid method invocation `call` with arguments".

  A synonymous of doing this with a regular `map()` would be:
  ```groovy
  workflow {
    Channel.fromList([["id", [input: "foo"], "passthrough"]])
      | map{ tup ->
        def (id, data) = tup
        [id, data + [arg: 10]] + tup.drop(2)
      }
  }
  ```

* `correction/cellbender_remove_background`: Eliminating technical artifacts from high-throughput single-cell RNA sequencing data.

## MAJOR CHANGES

* `workflows/utils/DataFlowHelper.nf`: Added helper functions `setWorkflowArguments()` and `getWorkflowArguments()` to split the data field of a channel event into a hashmap. Example usage:
  ```groovy
  | setWorkflowArguments(
    pca: [ "input": "input", "obsm_output": "obsm_pca" ]
    integration: [ "obs_covariates": "obs_covariates", "obsm_input": "obsm_pca" ]
  )
  | getWorkflowArguments("pca")
  | pca
  | getWorkflowArguments("integration")
  | integration
  ```

* `mapping/cellranger_count`: Allow passing both directories as well as individual fastq.gz files as inputs.

* `convert/from_10xh5_to_h5mu`: Allow reading in QC metrics, use gene ids as `.obs_names` instead of gene symbols.

## MINOR CHANGES

* `dimred/umap`: Streamline UMAP parameters by adding `--obsm_output` parameter to allow choosing the output `.obsm` slot.

* `workflows/multiomics/integration`: Added arguments for tuning the various output slots of the integration pipeline, namely `--obsm_pca`, `--obsm_integrated`, `--uns_neighbors`, `--obsp_neighbor_distances`, `--obsp_neighbor_connectivities`, `--obs_cluster`, `--obsm_umap`.

* Switch to Viash 0.6.1.

## BUG FIXES

* `convert/from_bd_to_10x_molecular_barcode_tags`: Replaced UTF8 characters with ASCII. OpenJDK 17 or lower might throw the following exception when trying to read a UTF8 file: `java.nio.charset.MalformedInputException: Input length = 1`.

* `integrate/concat`: Overriding sample name in .obs no longer raises `AttributeError`.

# openpipeline 0.5.0

Major redesign of the integration and multiomic workflows. Current list of workflows:

* `ingestion/bd_rhapsody`: A generic pipeline for running BD Rhapsody WTA or Targeted mapping, with support for AbSeq, VDJ and/or SMK.

* `ingestion/cellranger_mapping`: A pipeline for running Cell Ranger mapping.

* `ingestion/demux`: A generic pipeline for running bcl2fastq, bcl-convert or Cell Ranger mkfastq.

* `multiomics/rna_singlesample`: Processing unimodal single-sample RNA transcriptomics data.

* `multiomics/rna_multisample`: Processing unimodal multi-sample RNA transcriptomics data.

* `multiomics/integration`: A pipeline for demultiplexing multimodal multi-sample RNA transcriptomics data.

* `multiomics/full_pipeline`: A pipeline to analyse multiple multiomics samples.

## BREAKING CHANGES

* Many components: Renamed `.var["gene_ids"]` and `.var["feature_types"]` to `.var["gene_id"]` and `.var["feature_type"]`.

## DEPRECATED

* `convert/from_10xh5_to_h5ad` and `convert/from_bdrhap_to_h5ad`: Removed h5ad based components.

* `mapping/bd_rhapsody_wta` and `workflows/ingestion/bd_rhapsody_wta`: Deprecated in favour for more generic `mapping/bd_rhapsody` and `workflows/ingestion/bd_rhapsody` pipelines.

* `convert/from_csv_to_h5mu`: Disable until it is needed again.

* `integrate/concat`: Deprecated `"concat"` option for `--other_axis_mode`.

## NEW COMPONENTS

* `graph/bbknn`: Batch balanced KNN.

* `transform/scaling`: Scale data to unit variance and zero mean.

* `mapping/bd_rhapsody`: Added generic component for running the BD Rhapsody WTA or Targeted analysis, with support for AbSeq, VDJ and/or SMK.

* `integrate/harmony` and `integrate/harmonypy`: Run a Harmony integration analysis (R-based and Python-based, respectively).

* `integrate/scanorama`: Use Scanorama to integrate different experiments.

* `reference/make_reference`: Download a transcriptomics reference and preprocess it (adding ERCC spikeins and filtering with a regex).

* `reference/build_bdrhap_reference`: Compile a reference into a STAR index in the format expected by BD Rhapsody.

## NEW WORKFLOWS

* `workflows/ingestion/bd_rhapsody`: Added generic workflow for running the BD Rhapsody WTA or Targeted analysis, with support for AbSeq, VDJ and/or SMK.

* `workflows/multiomics/full_pipeline`: Implement pipeline for processing multiple multiomics samples.

## NEW FUNCTIONALITY

* `convert/from_bdrhap_to_h5mu`: Added support for being able to deal with WTA, Targeted, SMK, AbSeq and VDJ data.

* `integrate/concat`: Added `"move"` option to `--other_axis_mode`, which allows merging `.obs` and `.var` by only keeping elements of the matrices which are the same in each of the samples, moving the conflicting values to `.varm` or `.obsm`.

## MAJOR CHANGES

* Multiple components: Update to anndata 0.8 with mudata 0.2.0. This means that the format of the `.h5mu` files have changed.

* `multiomics/rna_singlesample`: Move transformation counts into layers instead of overwriting `.X`.

* Updated to Viash 0.6.0.

## MINOR CHANGES

* `velocity/velocyto`: Allow configuring memory and parallellisation.

* `cluster/leiden`: Add `--obsp_connectivities` parameter to allow choosing the output slot.

* `workflows/multiomics/rna_singlesample`, `workflows/multiomics/rna_multisample` and `workflows/multiomics/integration`: Allow choosing the output paths.

* `neighbors/bbknn` and `neighbors/find_neighbors`: Add parameters for choosing the input/output slots.

* `dimred/pca` and `dimred/umap`: Add parameters for choosing the input/output slots.

* `integrate/concat`: Optimize concat performance by adding multiprocessing and refactoring functions.

* `workflows/multimodal_integration`: Add `obs_covariates` argument to pipeline.

## BUG FIXES

* Several components: Revert using slim versions of containers because they do not provide the tools to run nextflow with trace capabilities.

* `integrate/concat`: Fix an issue where joining boolean values caused `TypeError`.

* `workflows/multiomics/rna_multisample`, `workflows/multiomics/rna_singlesample` and `workflows/multiomics/integration`: Use nextflow trace reporting when running integration tests.


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
