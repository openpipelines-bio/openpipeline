import re
import sys
from subprocess import CalledProcessError

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

## VIASH START
meta = {
    "executable": "./target/executable/perturbation/geneformer_virtual_cells/geneformer_virtual_cells",
    "resources_dir": "src/utils",
    "config": "./src/perturbation/geneformer_virtual_cells/config.vsh.yaml",
}
## VIASH END

PAD, CLS, EOS = 0, 2, 3

# gene_4 shares token 40 with gene_3; gene_5 has no token
TOKENS = [10, 20, 30, 40, 40, -1]
SYMBOLS = ["A", "B", "C", "D", "D2", "E"]


def settings(special_token, size):
    return {
        "model_version": "V2" if special_token else "V1",
        "model_input_size": size,
        "special_token": special_token,
        "pad_token": PAD,
        "cls_token": CLS if special_token else -1,
        "eos_token": EOS if special_token else -1,
    }


@pytest.fixture
def make_input(random_h5mu_path):
    def wrapper(special_token=False, size=3):
        # the rank order in the tokens is set by hand and differs from the
        # count order, as normalisation by the gene medians can make it differ:
        # A has the most counts of the genes with a token but ranks last, E has
        # the most counts of all but no token
        counts = np.array(
            [
                [9, 5, 0, 7, 2, 100],  # cell_a
                [0, 4, 6, 0, 0, 1],  # cell_b
            ]
        )
        uncropped = np.array(
            [
                [20, 40, 10, PAD],  # cell_a ranked: B, D(+D2), A
                [30, 20, PAD, PAD],  # cell_b ranked: C, B
            ],
            dtype=np.uint16,
        )
        rna = ad.AnnData(
            X=csr_matrix(counts.astype(np.float32)),
            obs=pd.DataFrame(index=["cell_a", "cell_b"]),
            var=pd.DataFrame(
                {
                    "gene_symbol": SYMBOLS,
                    "geneformer_token": TOKENS,
                    "geneformer_ensembl_id": [
                        "ENSG0",
                        "ENSG1",
                        "ENSG2",
                        "ENSG3",
                        "ENSG3",
                        "",
                    ],
                },
                index=[f"gene_{i}" for i in range(6)],
            ),
        )
        rna.obsm["geneformer_tokens_uncropped"] = uncropped
        rna.uns["geneformer_tokenize"] = settings(special_token, size)
        path = random_h5mu_path()
        mu.MuData({"rna": rna}).write_h5mu(path)
        return path

    return wrapper


def run(run_component, input_path, output_path, *extra):
    run_component(["--input", input_path, "--output", output_path, *extra])
    return mu.read_h5mu(output_path).mod["rna"]


def tokens_of(virtual, name):
    return virtual[name].obsm["geneformer_tokens"][0]


def test_knockouts_are_the_model_input_and_refill_the_crop(
    run_component, make_input, random_h5mu_path
):
    virtual = run(
        run_component,
        make_input(special_token=False, size=2),
        random_h5mu_path(),
        "--var_gene_names",
        "gene_symbol",
    )

    # model input of cell_a is B, D: A (below the crop) and E (no token) are
    # not knocked out, whatever their counts. D2 shares D's token, the gene is
    # named after the first one in .var
    cell_a = virtual[virtual.obs["perturbation_cell_id"] == "cell_a"]
    assert list(cell_a.obs["perturbation_gene_name"]) == ["B", "D"]
    assert list(cell_a.obs["perturbation_gene_id"]) == ["ENSG1", "ENSG3"]

    # knocking out a gene of the model input pulls A into it
    assert list(tokens_of(virtual, "cell_a@ENSG1")) == [40, 10]
    assert list(tokens_of(virtual, "cell_a@ENSG3")) == [20, 10]

    cell_b = virtual[virtual.obs["perturbation_cell_id"] == "cell_b"]
    assert list(cell_b.obs["perturbation_gene_name"]) == ["C", "B"]
    assert list(tokens_of(virtual, "cell_b@ENSG2")) == [20, PAD]
    assert list(tokens_of(virtual, "cell_b@ENSG1")) == [30, PAD]
    assert virtual.n_obs == 4


def test_model_input_size_bounds_the_knockouts(
    run_component, make_input, random_h5mu_path
):
    virtual = run(run_component, make_input(size=3), random_h5mu_path())
    cell_a = virtual[virtual.obs["perturbation_cell_id"] == "cell_a"]
    # all three ranked genes fit in the model input now
    assert list(cell_a.obs["perturbation_gene_id"]) == ["ENSG1", "ENSG3", "ENSG0"]
    assert list(tokens_of(virtual, "cell_a@ENSG0")) == [20, 40, PAD]


def test_special_tokens_and_top_n(run_component, make_input, random_h5mu_path):
    virtual = run(
        run_component,
        make_input(special_token=True, size=4),
        random_h5mu_path(),
        "--top_n_genes",
        "1",
        "--obsm_output",
        "my_tokens",
        "--output_compression",
        "gzip",
    )
    assert list(virtual.obs_names) == ["cell_a@ENSG1", "cell_b@ENSG2"]
    # ranked B, D, A minus B, cropped to 4 - 2 and wrapped in <cls> ... <eos>
    assert list(virtual.obsm["my_tokens"][0]) == [CLS, 40, 10, EOS]
    assert list(virtual.obsm["my_tokens"][1]) == [CLS, 20, EOS, PAD]
    # gene names default to the .var index
    assert list(virtual.obs["perturbation_gene_name"]) == ["gene_1", "gene_2"]
    assert virtual.uns["geneformer_tokenize"]["model_input_size"] == 4


def test_tokens_not_matching_var_error(run_component, make_input, random_h5mu_path):
    path = make_input()
    mdata = mu.read_h5mu(path)
    # token 99 belongs to no gene in .var
    mdata.mod["rna"].obsm["geneformer_tokens_uncropped"][1] = [99, 20, PAD, PAD]
    mdata.write_h5mu(path)
    with pytest.raises(CalledProcessError) as err:
        run(run_component, path, random_h5mu_path())
    assert re.search(
        r"Tokens \[99\] of cell 'cell_b'", err.value.stdout.decode("utf-8")
    )


def test_untokenized_input_errors(run_component, make_input, random_h5mu_path):
    with pytest.raises(CalledProcessError) as err:
        run(
            run_component,
            make_input(),
            random_h5mu_path(),
            "--obsm_input_uncropped",
            "missing",
        )
    assert re.search(r"--obsm_output_uncropped first", err.value.stdout.decode("utf-8"))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
