# openpipeline 0.2.1

## Minor changes

* Translate `bd_rhapsody_extracth5ad` from R into Python script.

* `bd_rhapsody_wta`: Remove temporary directory after execution.


# openpipeline 0.2.0

## Functionality

* Added `convert_10x_to_h5ad` and `download_10x_dataset` components.

## Minor changes
* Workflow `bd_rhapsody_wta`: Minor change to workflow to allow for easy processing of multiple samples with a tsv.

* Component `bd_rhapsody_wta`: Added more parameters, `--parallel` and `--timestamps`.

* Added `pbmc_1k_protein_v3` as a test resource.


# openpipeline 0.1.0

* Initial release containing only a `bd_rhapsody_wta` pipeline and corresponding components.
