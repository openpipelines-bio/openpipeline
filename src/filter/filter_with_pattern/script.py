import mudata as mu
import sys
import pandas as pd

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "rna",
    "output": "output.h5mu",
    "pattern": ["MIR\\d+", "AL\\d+", "LINC\\d+", "AC\\d+", "AP\\d+"],
    "var_name_filter": "filter_with_pattern",
    "var_gene_names": "gene_symbol",
    "do_subset": True,
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
### VIASH END


sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading input data from %s, modality %s", par["input"], par["modality"])
modality_data = mu.read_h5ad(par["input"], mod=par["modality"])

logger.info("\tUnfiltered data: %s", modality_data)

gene_names = (
    modality_data.var[par["var_gene_names"]].values
    if par["var_gene_names"]
    else modality_data.var_names.values
)

pattern_string = "|".join(par["pattern"])

# Filter genes based on regex pattern
genes_to_keep = list(~pd.Series(gene_names).str.contains(pattern_string, na=False))

logger.info(
    f"Number of genes to keep: {sum(genes_to_keep)} out of {len(genes_to_keep)}"
)
modality_data.var[par["var_name_filter"]] = genes_to_keep

if par["do_subset"]:
    modality_data = modality_data[:, genes_to_keep]
    logger.info("\tFiltered data: %s", modality_data)

logger.info(
    "Writing output data to %s with compression %s",
    par["output"],
    par["output_compression"],
)
write_h5ad_to_h5mu_with_compression(
    par["output"],
    par["input"],
    par["modality"],
    modality_data,
    par["output_compression"],
)


logger.info("Finished")
