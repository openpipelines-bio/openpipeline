import sys
import pytest
from mudata import read_h5mu

## VIASH START
meta = {
    "config": "src/velocity/scvelo/config.vsh.yaml",
    "executable": "./target/executable/velocity/scvelo/scvelo",
    "resources_dir": "resources_test/rna_velocity/velocyto_processed/",
}
## VIASH END

input_h5mu = f"{meta['resources_dir']}/velocyto.h5mu"


def test_scvelo(run_component, tmp_path):
    output_dir = tmp_path / "foo"
    output_h5mu = tmp_path / "output.h5mu"
    run_component(
        [
            "--input",
            input_h5mu,
            "--modality",
            "velocyto",
            "--output",
            str(output_dir),
            "--output_h5mu",
            output_h5mu,
            "--layer_spliced",
            "velo_spliced",
            "--layer_unspliced",
            "velo_unspliced",
            "--layer_ambiguous",
            "velo_ambiguous",
            "--output_compression",
            "gzip",
        ]
    )

    assert output_dir.is_dir()
    assert (output_dir / "scvelo_proportions.pdf").is_file()
    assert (output_dir / "scvelo_embedding.pdf").is_file()
    assert (output_dir / "scvelo_graph.pdf").is_file()
    assert (output_dir / "proportions.txt").is_file()
    assert output_h5mu.is_file()

    output_data = read_h5mu(output_h5mu)
    items_per_slot = {
        "uns": (
            "recover_dynamics",
            "velocity_params",
            "velocity_graph",
            "velocity_graph_neg",
        ),
        "varm": ("loss",),
        "obsm": ("velocity_pca",),
        "layers": (
            "Ms",
            "Mu",
            "fit_t",
            "fit_tau",
            "fit_tau_",
            "velocity",
            "velocity_u",
        ),
    }
    for dict_slot, dict_items in items_per_slot.items():
        for dict_item in dict_items:
            assert (
                getattr(output_data.mod["velocyto"], dict_slot).get(dict_item)
                is not None
            ), f"Expected {dict_item} to be present in {dict_slot}"
            del getattr(output_data.mod["velocyto"], dict_slot)[dict_item]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
