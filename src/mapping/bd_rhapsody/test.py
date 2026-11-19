import sys
from pathlib import Path
from uuid import uuid4

import mudata as md
import pytest

## VIASH START
meta = {
    "executable": "target/executable/mapping/bd_rhapsody/bd_rhapsody",
    "resources_dir": "resources_test/",
    "config": "src/mapping/bd_rhapsody/config.vsh.yaml",
    "cpus": 2,
    "memory_gb": 10,
}
## VIASH END

resources_dir = Path(meta["resources_dir"])
raw_dir = resources_dir / "raw"
reference_archive = resources_dir / "reference_bd_rhapsody.tar.gz"
abseq_reference = raw_dir / "BDAbSeq_ImmuneDiscoveryPanel.fasta"
reads = [
    raw_dir / f"{library}_S1_L432_{read}_001_subset.fastq.gz"
    for library in ["12WTA", "12ABC"]
    for read in ["R1", "R2"]
]


@pytest.fixture
def random_path(tmp_path):
    def wrapper(extension=None):
        extension = f".{extension}" if extension else ""
        return tmp_path / f"{uuid4()}{extension}"

    return wrapper


def test_wta_abseq(run_component, random_path):
    output_dir = random_path()
    outputs = {
        "output_seurat": random_path(extension="rds"),
        "output_mudata": random_path(extension="h5mu"),
        "metrics_summary": random_path(extension="csv"),
        "pipeline_report": random_path(extension="html"),
        "dimred_tsne": random_path(extension="csv"),
        "dimred_umap": random_path(extension="csv"),
    }

    run_component(
        [
            "--reads",
            ";".join(str(read) for read in reads),
            "--reference_archive",
            reference_archive,
            "--abseq_reference",
            abseq_reference,
            "--exact_cell_count",
            "4900",
            "--output_dir",
            output_dir,
        ]
        + [arg for key, path in outputs.items() for arg in (f"--{key}", path)]
    )

    for file_name in [
        "sample_Bioproduct_Stats.csv",
        "sample_Metrics_Summary.csv",
        "sample_Pipeline_Report.html",
        "sample_RSEC_MolsPerCell_MEX.zip",
        "sample_RSEC_MolsPerCell_Unfiltered_MEX.zip",
        "sample_Seurat.rds",
        "sample.cellismo",
        "sample_biotypeCounts.csv",
    ]:
        assert (output_dir / file_name).is_file(), (
            f"Output directory should contain '{file_name}'."
        )

    for key, path in outputs.items():
        assert path.is_file(), f"Output file for --{key} should be created."

    data = md.read_h5mu(outputs["output_mudata"])
    assert data.n_obs == 4900, "Number of cells is incorrect."
    assert "rna" in data.mod, "RNA modality is missing."
    assert "prot" in data.mod, "Protein modality is missing."

    data_rna = data.mod["rna"]
    assert data_rna.n_vars > 1000, "Number of genes is incorrect."
    assert data_rna.X.sum(axis=1).min() > 0, "Number of reads per cell is incorrect."
    assert data_rna.var["Raw_Reads"].sum() > 100000, "Number of reads is incorrect."


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
