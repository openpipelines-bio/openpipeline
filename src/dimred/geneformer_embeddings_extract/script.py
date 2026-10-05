import os
import sys
import tempfile

import mudata as mu
import numpy as np
import torch
from datasets import Dataset
from geneformer import EmbExtractor

## VIASH START
par = {
    "input": "tokenized.h5mu",
    "modality": "rna",
    "obsm_input": "geneformer_tokens",
    "uns_input": "geneformer_tokenize",
    "model": "resources_test/geneformer/Geneformer-V1-10M",
    "model_version": "V1",
    "emb_mode": "cell",
    "emb_layer": 0,
    "forward_batch_size": 32,
    "output": "output.h5mu",
    "obsm_output": "X_geneformer",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils", "cpus": 4, "temp_dir": None}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402
from compress_h5mu import write_h5ad_to_h5mu_with_compression  # noqa: E402

logger = setup_logger()

# Dataset column used to put the embeddings back in input order:
# EmbExtractor sorts the cells by length before the forward pass.
ROW_LABEL = "cell_index"


def get_tokens(adata):
    if par["obsm_input"] not in adata.obsm:
        raise ValueError(
            f"'{par['obsm_input']}' not found in .obsm. Tokenize the cells with "
            f"transform/geneformer_tokenize first. Available: {list(adata.obsm)}"
        )
    if par["uns_input"] not in adata.uns:
        raise ValueError(
            f"'{par['uns_input']}' not found in .uns. Tokenize the cells with "
            f"transform/geneformer_tokenize first. Available: {list(adata.uns)}"
        )
    settings = adata.uns[par["uns_input"]]
    if settings["model_version"] != par["model_version"]:
        raise ValueError(
            f"The cells were tokenized for model version {settings['model_version']}, "
            f"but --model_version is {par['model_version']}. Tokenize again with "
            "the model version of --model."
        )

    tokens = np.asarray(adata.obsm[par["obsm_input"]])
    lengths = (tokens != settings["pad_token"]).sum(axis=1)
    if np.any(lengths == 0):
        empty = adata.obs_names[lengths == 0][:5].tolist()
        raise ValueError(
            f"{int((lengths == 0).sum())} cells have no tokens, e.g. {empty}"
        )
    return tokens, lengths


def main():
    if par["emb_mode"] == "cls" and par["model_version"] == "V1":
        raise ValueError(
            "--emb_mode cls needs a <cls> token, which V1 models do not have"
        )
    if not os.path.isdir(par["model"]):
        raise ValueError(f"--model must be a directory, got {par['model']}")

    nproc = int(meta["cpus"] or 1)
    if torch.cuda.is_available():
        logger.info(
            "Running the forward pass on cuda (%s)", torch.cuda.get_device_name(0)
        )
    else:
        torch.set_num_threads(nproc)
        logger.warning(
            "No CUDA device visible, running the forward pass on the CPU with %i "
            "threads. This is orders of magnitude slower than on a GPU.",
            nproc,
        )

    logger.info("Reading modality '%s' from %s", par["modality"], par["input"])
    adata = mu.read_h5ad(par["input"], mod=par["modality"])
    tokens, lengths = get_tokens(adata)
    logger.info(
        "%i cells, %i to %i tokens per cell", adata.n_obs, lengths.min(), lengths.max()
    )

    dataset = Dataset.from_dict(
        {
            "input_ids": [row[:n].tolist() for row, n in zip(tokens, lengths)],
            "length": lengths.tolist(),
            ROW_LABEL: list(range(adata.n_obs)),
        }
    )

    extractor = EmbExtractor(
        model_type="Pretrained",
        num_classes=0,
        emb_mode=par["emb_mode"],
        max_ncells=None,
        emb_layer=par["emb_layer"],
        emb_label=[ROW_LABEL],
        forward_batch_size=par["forward_batch_size"],
        model_version=par["model_version"],
        nproc=nproc,
    )

    with tempfile.TemporaryDirectory(dir=meta.get("temp_dir")) as temp_dir:
        dataset_path = os.path.join(temp_dir, "tokens.dataset")
        dataset.save_to_disk(dataset_path)
        embeddings = extractor.extract_embs(
            par["model"], dataset_path, temp_dir, "embeddings"
        )

    embeddings = embeddings.sort_values(ROW_LABEL)
    if not np.array_equal(embeddings[ROW_LABEL].to_numpy(), np.arange(adata.n_obs)):
        raise RuntimeError("The model did not return one embedding per cell")
    embedding_matrix = embeddings.drop(columns=[ROW_LABEL]).to_numpy(dtype=np.float32)

    adata.obsm[par["obsm_output"]] = embedding_matrix
    logger.info(
        "Stored %i x %i embeddings in .obsm['%s']",
        *embedding_matrix.shape,
        par["obsm_output"],
    )

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
