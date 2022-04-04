# Conversion

## Goal
Pipeline to convert the different formats of raw count matrices to the input format required by the pipelines. 

## Input

- ### ___input_type___

- ### ___{input_name}.{extension}___
Takes a single or multiple input files of the format described below and the layer variable where the final data will be stored.

#### __input_type == "10x_h5"__ (extension = h5)

The 10x h5 format, supports both the legacy single modality version from cellranger 2 and the multi-modality version from cellranger 3 and upwards. 

- __modality__ (e.g. "Gene Expression", "Antibody Capture")
    
The feature type to be extracted to the final h5ad.

#### __input_type == "mtx format"__ (extenstion = none)

The standard mtx format.

#### __input_type == "csv"__ (extension = tsv,csv,dge)

Count matrices in different matrix forms.

####  __fcs format__ (extension = fcs)

#### __seurat format__ (extension = rds)


## Output

- ### {input_name}-converted.h5ad

    - X
    - {raw}/X

    The raw count data that was in the original files.
