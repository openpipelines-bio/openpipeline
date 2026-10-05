import os
import sys
import tempfile
from pathlib import Path

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
from geneformer import TranscriptomeTokenizer
from scipy.sparse import csr_matrix, issparse

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "input_layer": None,
    "var_gene_ids": None,
    "sanitize_ensembl_ids": True,
    "allow_non_integer": False,
    "gene_map": None,
    "var_gene_names": None,
    "gene_map_symbol_column": "gene_symbol",
    "gene_map_ensembl_column": "ensembl_id",
    "model_version": "V2",
    "output": "output.h5mu",
    "obsm_output": "geneformer_tokens",
    "obsm_output_uncropped": None,
    "var_output_gene_ids": "geneformer_ensembl_id",
    "var_output_tokens": "geneformer_token",
    "uns_output": "geneformer_tokenize",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils", "cpus": 4, "temp_dir": None}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402
from compress_h5mu import write_h5ad_to_h5mu_with_compression  # noqa: E402
from set_var_index import strip_version_number  # noqa: E402

logger = setup_logger()

# Column names of the AnnData staged for the tokenizer. The tokenizer reads
# "ensembl_id" from .var and "n_counts" from .obs; the row number is carried
# through as a custom attribute so the tokens can be put back in input order
# (the tokenizer drops cells without any gene in its dictionary).
ROW_ATTRIBUTE = "geneformer_row"


def get_counts(adata):
    if par["input_layer"]:
        if par["input_layer"] not in adata.layers:
            raise ValueError(
                f"Layer '{par['input_layer']}' not found in .layers. "
                f"Available: {list(adata.layers)}"
            )
        counts = adata.layers[par["input_layer"]]
    else:
        counts = adata.X
    counts = counts if issparse(counts) else csr_matrix(counts)
    return csr_matrix(counts)


def check_raw_counts(counts):
    if counts.sum() == 0:
        raise ValueError(
            "The count matrix is all zeros. Geneformer ranks genes by their raw "
            "counts; point --input_layer at the raw counts."
        )
    sample = counts.data[:1_000_000]
    if not np.allclose(sample, np.round(sample)):
        if not par["allow_non_integer"]:
            raise ValueError(
                "The count matrix holds non-integer values, which usually means "
                "it was normalized, while Geneformer needs raw counts. Point "
                "--input_layer at the raw counts, or pass --allow_non_integer if "
                "these really are counts."
            )
        logger.warning("Non-integer counts accepted because of --allow_non_integer")


def get_ensembl_ids(adata):
    if par["gene_map"]:
        if par["var_gene_names"]:
            if par["var_gene_names"] not in adata.var.columns:
                raise ValueError(
                    f"Column '{par['var_gene_names']}' not found in .var. "
                    f"Available: {list(adata.var.columns)}"
                )
            symbols = adata.var[par["var_gene_names"]].astype(str)
        else:
            symbols = adata.var_names.to_series().astype(str)

        gene_map = pd.read_csv(par["gene_map"])
        columns = [par["gene_map_symbol_column"], par["gene_map_ensembl_column"]]
        missing = [col for col in columns if col not in gene_map.columns]
        if missing:
            raise ValueError(
                f"Columns {missing} not found in {par['gene_map']}. "
                f"Available: {list(gene_map.columns)}"
            )
        lookup = (
            gene_map[columns]
            .dropna()
            .drop_duplicates(subset=par["gene_map_symbol_column"])
            .set_index(par["gene_map_symbol_column"])[par["gene_map_ensembl_column"]]
        )
        ensembl_ids = symbols.map(lookup)
        logger.info(
            "%i of %i gene symbols map to an Ensembl id via %s",
            ensembl_ids.notna().sum(),
            adata.n_vars,
            par["gene_map"],
        )
    elif par["var_gene_ids"]:
        if par["var_gene_ids"] not in adata.var.columns:
            raise ValueError(
                f"Column '{par['var_gene_ids']}' not found in .var. "
                f"Available: {list(adata.var.columns)}"
            )
        ensembl_ids = adata.var[par["var_gene_ids"]]
    else:
        ensembl_ids = adata.var_names.to_series()

    ensembl_ids = ensembl_ids.astype(object).where(ensembl_ids.notna(), "")
    ensembl_ids = pd.Series(ensembl_ids.astype(str).to_numpy(), index=adata.var_names)
    if par["sanitize_ensembl_ids"]:
        ensembl_ids = strip_version_number(ensembl_ids)
    return ensembl_ids.str.upper()


def crop_and_pad(cells, rows, n_obs, width, pad_token):
    matrix = np.full((n_obs, width), pad_token, dtype=np.int64)
    for row, tokens in zip(rows, cells):
        matrix[row, : len(tokens)] = tokens
    return matrix


def format_model_input(tokens, tokenizer):
    # Same cropping as TranscriptomeTokenizer.create_dataset: keep the highest
    # ranked genes and, for V2, wrap them in <cls> ... <eos>.
    if tokenizer.special_token:
        cropped = tokens[: tokenizer.model_input_size - 2]
        return np.concatenate(
            [
                [tokenizer.gene_token_dict["<cls>"]],
                cropped,
                [tokenizer.gene_token_dict["<eos>"]],
            ]
        )
    return tokens[: tokenizer.model_input_size]


def main():
    logger.info("Reading modality '%s' from %s", par["modality"], par["input"])
    adata = mu.read_h5ad(par["input"], mod=par["modality"])

    counts = get_counts(adata)
    check_raw_counts(counts)
    ensembl_ids = get_ensembl_ids(adata)

    tokenizer = TranscriptomeTokenizer(
        custom_attr_name_dict={ROW_ATTRIBUTE: ROW_ATTRIBUTE},
        nproc=int(meta["cpus"] or 1),
        model_version=par["model_version"],
    )
    pad_token = tokenizer.gene_token_dict["<pad>"]
    logger.info(
        "Tokenizer settings: model_version=%s, model_input_size=%i, special_token=%s",
        par["model_version"],
        tokenizer.model_input_size,
        tokenizer.special_token,
    )

    staged = ad.AnnData(
        X=counts,
        obs=pd.DataFrame(
            {
                "n_counts": np.asarray(counts.sum(axis=1)).ravel(),
                ROW_ATTRIBUTE: np.arange(adata.n_obs),
            },
            index=pd.Index(np.arange(adata.n_obs).astype(str)),
        ),
        var=pd.DataFrame(
            {"ensembl_id": ensembl_ids.to_numpy()},
            index=pd.Index(np.arange(adata.n_vars).astype(str)),
        ),
    )

    with tempfile.TemporaryDirectory(dir=meta.get("temp_dir")) as temp_dir:
        staged.write_h5ad(os.path.join(temp_dir, "staged.h5ad"))
        del staged
        ranked_cells, cell_metadata, _ = tokenizer.tokenize_files(
            Path(temp_dir), file_format="h5ad"
        )
    rows = np.asarray(cell_metadata[ROW_ATTRIBUTE], dtype=int)

    n_empty = adata.n_obs - len(rows)
    if n_empty:
        empty = np.setdiff1d(np.arange(adata.n_obs), rows)
        raise ValueError(
            f"{n_empty} cells express no gene of the Geneformer token dictionary, "
            f"e.g. {list(adata.obs_names[empty[:5]])}. Remove them before "
            "tokenizing, they cannot be embedded."
        )

    model_input = [format_model_input(tokens, tokenizer) for tokens in ranked_cells]
    max_token = max(tokenizer.gene_token_dict.values())
    dtype = np.min_scalar_type(max_token)
    adata.obsm[par["obsm_output"]] = crop_and_pad(
        model_input, rows, adata.n_obs, tokenizer.model_input_size, pad_token
    ).astype(dtype)
    lengths = np.array([len(tokens) for tokens in model_input])
    logger.info(
        "Tokenized %i cells, %i to %i tokens per cell (median %i)",
        adata.n_obs,
        lengths.min(),
        lengths.max(),
        np.median(lengths),
    )

    if par["obsm_output_uncropped"]:
        width = max(len(tokens) for tokens in ranked_cells)
        adata.obsm[par["obsm_output_uncropped"]] = crop_and_pad(
            ranked_cells, rows, adata.n_obs, width, pad_token
        ).astype(dtype)
        logger.info(
            "Stored uncropped ranked tokens in .obsm['%s'], width %i",
            par["obsm_output_uncropped"],
            width,
        )

    # Same gene id resolution as the tokenizer, so that a gene in .var can be
    # matched to its token without the geneformer dictionaries.
    collapsed = ensembl_ids.map(tokenizer.gene_mapping_dict)
    tokens = collapsed.map(tokenizer.gene_token_dict)
    adata.var[par["var_output_gene_ids"]] = collapsed.fillna("").to_numpy()
    adata.var[par["var_output_tokens"]] = tokens.fillna(-1).astype(np.int64).to_numpy()
    logger.info(
        "%i of %i genes are in the Geneformer token dictionary",
        (tokens.notna()).sum(),
        adata.n_vars,
    )

    adata.uns[par["uns_output"]] = {
        "model_version": par["model_version"],
        "model_input_size": int(tokenizer.model_input_size),
        "special_token": bool(tokenizer.special_token),
        "pad_token": int(pad_token),
        "cls_token": int(tokenizer.gene_token_dict.get("<cls>", -1)),
        "eos_token": int(tokenizer.gene_token_dict.get("<eos>", -1)),
    }

    logger.info("Writing output to %s", par["output"])
    write_h5ad_to_h5mu_with_compression(
        output_file=par["output"],
        h5mu=par["input"],
        modality_name=par["modality"],
        modality_data=adata,
        output_compression=par["output_compression"],
    )


if __name__ == "__main__":
    main()
