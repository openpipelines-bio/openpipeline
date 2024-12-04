import sys
import pytest
from pathlib import Path

## VIASH START
meta = {
    "name": "cellranger_count",
    "resources_dir": "resources_test"
}
## VIASH END

input = Path(meta["resources_dir"] + "/cellranger_tiny_fastq/cellranger_tiny_fastq/")
reference = meta["resources_dir"] + "/cellranger_tiny_fastq/cellranger_tiny_ref/"

def test_cellranger_count_with_folder(run_component, random_path):
    output = random_path() 
    run_component([
        "--input", input,
        "--reference", reference,
        "--output", output,
        "--lanes", "1",
    ])

    assert (output / "filtered_feature_bc_matrix.h5").is_file(), "No output was created."


def test_cellranger_count_with_fastq_files(run_component, random_path):
    output = random_path()
    run_component([
        "--input", input / "tinygex_S1_L001_R1_001.fastq.gz",
        "--input", input / "tinygex_S1_L001_R2_001.fastq.gz",
        "--reference", reference,
        "--output", output, 
    ])
    assert (output / "filtered_feature_bc_matrix.h5").is_file(), "No output was created."


@pytest.mark.parametrize("chemistry", ["auto", "SC3Pv2"])
def test_cellranger_chemistry(run_component, random_path, chemistry):
    output = random_path()
    run_component([
        "--input", input / "tinygex_S1_L001_R1_001.fastq.gz",
        "--input", input / "tinygex_S1_L001_R2_001.fastq.gz",
        "--reference", reference,
        "--output", output, 
        "--chemistry", chemistry,
    ])
    assert (output / "filtered_feature_bc_matrix.h5").is_file(), "No output was created."

def test_cellranger_no_bam(run_component, random_path):
    output = random_path()
    run_component([
        "--input", input / "tinygex_S1_L001_R1_001.fastq.gz",
        "--input", input / "tinygex_S1_L001_R2_001.fastq.gz",
        "--reference", reference,
        "--output", output, 
        "--generate_bam", "false",
    ])
    assert (output / "filtered_feature_bc_matrix.h5").is_file(), "No output was created."

def test_cellranger_no_secondary_analysis(run_component, random_path):
    output = random_path()
    run_component([
        "--input", input / "tinygex_S1_L001_R1_001.fastq.gz",
        "--input", input / "tinygex_S1_L001_R2_001.fastq.gz",
        "--reference", reference,
        "--output", output, 
        "--secondary_analysis", "false",
    ])
    assert (output / "filtered_feature_bc_matrix.h5").is_file(), "No output was created."

def test_cellranger_no_secondary_analysis(run_component, random_path):
    output = random_path()
    run_component([
        "--input", input / "tinygex_S1_L001_R1_001.fastq.gz",
        "--input", input / "tinygex_S1_L001_R2_001.fastq.gz",
        "--reference", reference,
        "--output", output, 
        "--secondary_analysis", "false",
    ])
    assert (output / "filtered_feature_bc_matrix.h5").is_file(), "No output was created."

def test_cellranger_exclude_introns(run_component, random_path):
    output = random_path()
    run_component([
        "--input", input / "tinygex_S1_L001_R1_001.fastq.gz",
        "--input", input / "tinygex_S1_L001_R2_001.fastq.gz",
        "--reference", reference,
        "--output", output, 
        "--include_introns", "false",
    ])
    assert (output / "filtered_feature_bc_matrix.h5").is_file(), "No output was created."

def test_cellranger_trim_reads(run_component, random_path):
    output = random_path()
    run_component([
        "--input", input / "tinygex_S1_L001_R1_001.fastq.gz",
        "--input", input / "tinygex_S1_L001_R2_001.fastq.gz",
        "--reference", reference,
        "--output", output, 
        "--r1_length", "100",
        "--r2_length", "100",
    ])
    assert (output / "filtered_feature_bc_matrix.h5").is_file(), "No output was created."

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))

