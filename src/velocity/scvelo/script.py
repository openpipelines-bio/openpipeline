import sys
import mudata
from contextlib import redirect_stdout
from pathlib import Path
import matplotlib as mpl

# Backwards compatibility for numpy 2.0
import numpy
numpy_module = sys.modules['numpy']
numpy_module.float_ = numpy.float64
sys.modules['numpy'] = numpy_module

# Backwards compatibility for scirpy
import scipy
scipy_module = sys.modules['scipy']
scipy_module.sparse._base._spbase.A = property(lambda self: self.toarray())

sys.modules['scipy'] = scipy_module

import scvelo

## VIASH START
from collections import defaultdict

def none_factory():
    return None

par = defaultdict(none_factory, {
    'input': './resources_test/rna_velocity/velocyto_processed/cellranger_tiny.loom',
    'output': './foo',
    'log_transform': True,
    'n_neighbors': 30
})
## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

mpl.rcParams['savefig.dpi']=150

# Script must be wrapped into a main function because scvelo spawn subprocesses
# and this fails when the functions are not wrapped.
def main():
    # Create output directory
    output_dir = Path(par['output'])
    output_dir.mkdir(parents=True, exist_ok=True)
    scvelo.settings.figdir = str(output_dir)


    # Calculate the sample name
    sample_name = par["output"].removesuffix(".loom")
    sample_name = Path(sample_name).name

    # Read the input data
    adata = scvelo.read(par['input'])

    # Save spliced vs unspliced proportions to file
    with (output_dir / "proportions.txt").open('w') as target:
        with redirect_stdout(target):
            scvelo.utils.show_proportions(adata)

    # Plot piecharts of spliced vs unspliced proportions
    scvelo.pl.proportions(adata, save=True, show=False)

    # Perform preprocessing
    scvelo.pp.filter_and_normalize(adata,
                                   min_counts=par["min_counts"],
                                   min_counts_u=par["min_counts_u"],
                                   min_cells=par["min_cells"],
                                   min_cells_u=par["min_cells_u"],
                                   min_shared_counts=par["min_shared_counts"],
                                   min_shared_cells=par["min_shared_cells"],
                                   n_top_genes=par["n_top_genes"],
                                   log=par["log_transform"])

    # Fitting
    scvelo.pp.moments(adata,
                      n_pcs=par["n_principal_components"],
                      n_neighbors=par["n_neighbors"])


    # Second step in velocyto calculations
    # Velocity calculation and visualization
    # From the scvelo manual:
    # The solution to the full dynamical model is obtained by setting mode='dynamical',
    # which requires to run scv.tl.recover_dynamics(adata) beforehand
    scvelo.tl.recover_dynamics(adata)
    scvelo.tl.velocity(adata, mode="dynamical")
    scvelo.tl.velocity_graph(adata)
    scvelo.pl.velocity_graph(adata, save=str(output_dir / "scvelo_graph.pdf"), show=False)

    # Plotting
    # TODO: add more here.
    scvelo.pl.velocity_embedding_stream(adata, save=str(output_dir / "scvelo_embedding.pdf"), show=False)

    # Create output
    ouput_data = mudata.MuData({'rna_velocity': adata})
    ouput_data.write_h5mu(output_dir / f"{sample_name}.h5mu", compression=par["output_compression"])

if __name__ == "__main__":
    main()