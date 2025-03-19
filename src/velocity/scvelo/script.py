import sys
import mudata
import anndata
import tempfile
import shutil
from contextlib import redirect_stdout
from pathlib import Path
import matplotlib as mpl
import scvelo

## VIASH START
from collections import defaultdict


def none_factory():
    return None


par = defaultdict(
    none_factory,
    {
        "input": "resources_test/rna_velocity/velocyto_processed/velocyto.h5mu",
        "modality": "velocyto",
        "output": "./foo",
        "output_h5mu": "output.h5mu",
        "log_transform": True,
        "n_neighbors": 30,
        "layer_spliced": "velo_spliced",
        "layer_unspliced": "velo_unspliced",
        "layer_ambiguous": "velo_ambiguous",
    },
)

meta = {"resources_dir": "src/utils", "temp_dir": "/tmp/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import compress_h5mu

logger = setup_logger()

mpl.rcParams["savefig.dpi"] = 150


# Script must be wrapped into a main function because scvelo spawn subprocesses
# and this fails when the functions are not wrapped.
def main():
    # Create output directory
    output_dir = Path(par["output"])
    output_dir.mkdir(parents=True, exist_ok=True)
    scvelo.settings.figdir = str(output_dir)

    # Load the input data
    adata_in = mudata.read_h5ad(par["input"], mod=par["modality"])

    # Create a copy of the data as input
    # as many scvelo functions do not take input layer arguments
    layers_mapping = {
        "spliced": par["layer_spliced"],
        "unspliced": par["layer_unspliced"],
        "ambiguous": par["layer_ambiguous"],
    }
    layer_data = {
        default: (adata_in.layers.get(arg_val) if arg_val else adata_in.layers[default])
        for default, arg_val in layers_mapping.items()
    }
    adata = anndata.AnnData(
        X=adata_in.X
        if not par["counts_layer"]
        else adata_in.layers[par["counts_layer"]],
        layers=layer_data,
    )

    # Calculate the sample name
    sample_name = par["output"].removesuffix(".h5mu")
    sample_name = Path(sample_name).name

    # Read the input data

    # Save spliced vs unspliced proportions to file
    with (output_dir / "proportions.txt").open("w") as target:
        with redirect_stdout(target):
            scvelo.utils.show_proportions(adata)

    # Plot piecharts of spliced vs unspliced proportions
    scvelo.pl.proportions(adata, save=True, show=False)

    # Perform preprocessing
    scvelo.pp.filter_and_normalize(
        adata,
        min_counts=par["min_counts"],
        min_counts_u=par["min_counts_u"],
        min_cells=par["min_cells"],
        min_cells_u=par["min_cells_u"],
        min_shared_counts=par["min_shared_counts"],
        min_shared_cells=par["min_shared_cells"],
        n_top_genes=par["n_top_genes"],
        log=par["log_transform"],
    )

    # Fitting
    scvelo.pp.moments(
        adata, n_pcs=par["n_principal_components"], n_neighbors=par["n_neighbors"]
    )

    # Second step in velocyto calculations
    # Velocity calculation and visualization
    # From the scvelo manual:
    # The solution to the full dynamical model is obtained by setting mode='dynamical',
    # which requires to run scv.tl.recover_dynamics(adata) beforehand
    scvelo.tl.recover_dynamics(adata)
    scvelo.tl.velocity(adata, mode="dynamical")
    scvelo.tl.velocity_graph(adata)
    scvelo.pl.velocity_graph(
        adata, save=str(output_dir / "scvelo_graph.pdf"), show=False
    )

    # Plotting
    # TODO: add more here.
    scvelo.pl.velocity_embedding_stream(
        adata, save=str(output_dir / "scvelo_embedding.pdf"), show=False
    )

    # Copy over slots to output
    for slot in ("obs", "var"):
        setattr(
            adata_in,
            slot,
            getattr(adata_in, slot)
            .assign(**getattr(adata, slot).to_dict())
            .convert_dtypes(),
        )
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
        setattr(
            adata_in,
            dict_slot,
            dict(
                getattr(adata_in, dict_slot),
                **{key_: getattr(adata, dict_slot)[key_] for key_ in dict_items},
            ),
        )
    with tempfile.NamedTemporaryFile(
        suffix=".h5mu", delete_on_close=False
    ) as temp_h5mu:
        shutil.copyfile(par["input"], temp_h5mu.name)
        # Create output
        mudata.write_h5ad(temp_h5mu.name, mod=par["modality"], data=adata_in)
        compression = par["output_compression"]

        if compression:
            compress_h5mu(temp_h5mu.name, par["output_h5mu"], compression=compression)
        else:
            shutil.move(temp_h5mu.name, par["output_h5mu"])


if __name__ == "__main__":
    main()
