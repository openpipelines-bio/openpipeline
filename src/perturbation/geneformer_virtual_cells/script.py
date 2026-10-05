import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

## VIASH START
par = {
    "input": "tokenized.h5mu",
    "modality": "rna",
    "obsm_input_uncropped": "geneformer_tokens_uncropped",
    "uns_input": "geneformer_tokenize",
    "var_tokens": "geneformer_token",
    "var_gene_ids": "geneformer_ensembl_id",
    "var_gene_names": None,
    "top_n_genes": None,
    "output": "output.h5mu",
    "obsm_output": "geneformer_tokens",
    "obs_output_cell_id": "perturbation_cell_id",
    "obs_output_gene_id": "perturbation_gene_id",
    "obs_output_gene_name": "perturbation_gene_name",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402

logger = setup_logger()


def check_inputs(adata):
    for key, slot in [
        (par["obsm_input_uncropped"], adata.obsm),
        (par["uns_input"], adata.uns),
    ]:
        if key not in slot:
            raise ValueError(
                f"'{key}' not found. Tokenize the cells with transform/geneformer_tokenize "
                "and --obsm_output_uncropped first."
            )
    for column in [par["var_tokens"], par["var_gene_ids"], par["var_gene_names"]]:
        if column and column not in adata.var.columns:
            raise ValueError(
                f"Column '{column}' not found in .var. Available: {list(adata.var.columns)}"
            )


def n_model_genes(settings):
    """Number of gene tokens in the model input, <cls> and <eos> excluded."""
    size = int(settings["model_input_size"])
    return size - 2 if bool(settings["special_token"]) else size


def model_input(ranked, settings):
    genes = ranked[: n_model_genes(settings)]
    if bool(settings["special_token"]):
        return np.concatenate([[settings["cls_token"]], genes, [settings["eos_token"]]])
    return genes


def token_to_gene(var_tokens):
    """Row in .var for every token; genes sharing a token (collapsed Ensembl
    ids) are named after the first one in .var."""
    rows = np.flatnonzero(var_tokens >= 0)
    tokens, first = np.unique(var_tokens[rows], return_index=True)
    return dict(zip(tokens.tolist(), rows[first].tolist()))


def main():
    logger.info("Reading modality '%s' from %s", par["modality"], par["input"])
    adata = mu.read_h5ad(par["input"], mod=par["modality"])
    check_inputs(adata)

    settings = adata.uns[par["uns_input"]]
    pad = settings["pad_token"]
    size = int(settings["model_input_size"])
    uncropped = np.asarray(adata.obsm[par["obsm_input_uncropped"]])

    var_tokens = adata.var[par["var_tokens"]].to_numpy().astype(np.int64)
    gene_of_token = token_to_gene(var_tokens)
    gene_ids = adata.var[par["var_gene_ids"]].astype(str).to_numpy()
    gene_names = (
        adata.var[par["var_gene_names"]].astype(str).to_numpy()
        if par["var_gene_names"]
        else adata.var_names.astype(str).to_numpy()
    )

    token_rows, cell_ids, knocked_ids, knocked_names = [], [], [], []
    for row, cell_id in enumerate(adata.obs_names):
        ranked = uncropped[row]
        ranked = ranked[ranked != pad]
        # knock out the genes the model sees: the gene tokens of the model
        # input, in rank order. The knockout removes the token from the
        # uncropped list, so the next gene moves into the model input.
        knocked = ranked[: n_model_genes(settings)][: par["top_n_genes"]]
        if knocked.size == 0:
            raise ValueError(f"Cell '{cell_id}' has no gene tokens.")
        unknown = [int(t) for t in knocked if int(t) not in gene_of_token]
        if unknown:
            raise ValueError(
                f"Tokens {unknown[:5]} of cell '{cell_id}' belong to no gene in "
                f".var['{par['var_tokens']}']; were the tokens made from this data?"
            )

        genes = np.array([gene_of_token[int(t)] for t in knocked])
        for position in range(knocked.size):
            token_rows.append(model_input(np.delete(ranked, position), settings))
        cell_ids.extend([cell_id] * genes.size)
        knocked_ids.extend(gene_ids[genes])
        knocked_names.extend(gene_names[genes])
        logger.info(
            "%s: %i ranked genes, %i knockouts", cell_id, ranked.size, genes.size
        )

    n_virtual = len(token_rows)
    tokens = np.full((n_virtual, size), pad, dtype=uncropped.dtype)
    for row, sequence in enumerate(token_rows):
        tokens[row, : sequence.size] = sequence

    obs = pd.DataFrame(
        {
            par["obs_output_cell_id"]: cell_ids,
            par["obs_output_gene_id"]: knocked_ids,
            par["obs_output_gene_name"]: knocked_names,
        },
        index=pd.Index([f"{c}@{g}" for c, g in zip(cell_ids, knocked_ids)]),
    )
    virtual = ad.AnnData(
        X=csr_matrix((n_virtual, 0), dtype=np.float32),
        obs=obs,
        obsm={par["obsm_output"]: tokens},
        uns={par["uns_input"]: dict(settings)},
    )
    logger.info(
        "Writing %i virtual cells from %i cells to %s",
        n_virtual,
        adata.n_obs,
        par["output"],
    )
    mu.MuData({par["modality"]: virtual}).write_h5mu(
        par["output"], compression=par["output_compression"]
    )


if __name__ == "__main__":
    main()
