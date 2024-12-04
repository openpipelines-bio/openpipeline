import mudata as mu
import anndata as ad
import sys
from pathlib import Path
import shutil

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modalities": ["rna", "prot"],
    "output": "output.h5mu"
}
meta = {
    
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import compress_h5mu

logger = setup_logger()

def main():
    modality_names = par['modalities']

    if len(modality_names) < 2:
        raise ValueError("Please provide two more more modalities.")
    
    obs_names = {}
    for mod_name in par['modalities']:
        try:
            modality = mu.read_h5ad(filename=par['input'], mod=mod_name)
        except KeyError:
            raise ValueError(f"Modality {mod_name} does not exist for file {par['input']}.")

        obs_names[mod_name] = modality.obs_names.copy()
        del modality
    
    intersected_index = None
    for mod_name, mod_index in obs_names.items():
        if intersected_index is None:
            intersected_index = mod_index
            continue
        intersected_index = intersected_index.intersection(mod_index)
    

    output_file = Path(par['output'])
    output_file_uncompressed = output_file.with_name(output_file.stem + "_uncompressed.h5mu")
    output_file_uncompressed.touch()

    mdata = mu.MuData({modality: ad.AnnData() for modality in modality_names})
    mdata.write(output_file_uncompressed, compression=par['output_compression'])
    
    for mod_name in modality_names:
        modality = mu.read_h5ad(filename=par['input'], mod=mod_name)
        intersected_modality = modality[intersected_index]
        mu.write_h5ad(output_file_uncompressed, data=intersected_modality, mod=mod_name)

    if par['output_compression']:
        compress_h5mu(output_file_uncompressed, output_file, compression=par['output_compression'])
        output_file_uncompressed.unlink()
    else:
        shutil.move(output_file_uncompressed, output_file)
    
if __name__ == "__main__":
    main()