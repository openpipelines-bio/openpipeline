# Openpipeline - data definition

All data will be stored and exchanged between the different components as h5mu. 

The h5mu will contain the standard field naming for the different data modalities as defined as _rna_, _prot_, _velo_, _atac_, _vdj_ and/or _crispr_.

Metadata will be added to the master h5mu object's obs slot. 

Filtering of the data will be performed using a _mask_ observation. 


## Object overview

data.h5mu

/obs: clinical/technical metadata over all cells

# - _mask_ (Boolean): The cells to be retained in the downstream analysis. ==> Filtering as alternative (keep the mask as a backdrop)

/mod: The different data modalities with their corresponding analysis.

- _rna_ (h5ad) (Gene level selection)
- _velo_ (h5ad)
- _prot_ (h5ad)
- _atac_ (h5ad)
- _crispr_ (h5ad)
- _vdj_ (h5ad) 

## Data usage

Every component will use the different modalities required from the original h5mu data and store the data on the appropriate modality. E.g. performing an analysis only on the _rna_ modality will store the data back to the _rna_ slot. When performing analysis using multiple modalities the data will be stored on the master h5mu object. 

## Final data export

The final h5mu object will be exported into the corresponding h5ad, seurat, and sce format collapsing the multi-modal if required. 
