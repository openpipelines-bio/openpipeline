import subprocess
from pathlib import Path

## VIASH START
meta = {
    "functionality_name": "cellranger_count",
    "resources_dir": "resources_test"
}
## VIASH END

# find common input files
resources_dir = Path(meta["resources_dir"])
input_dir = resources_dir / "cellranger_tiny_fastq" / "cellranger_tiny_fastq"
reference_index = resources_dir / "cellranger_tiny_fastq" / "cellranger_tiny_ref_v2_7_10_a/"
reference_gtf = resources_dir / "cellranger_tiny_fastq" / "cellranger_tiny_ref" / "genes" / "genes.gtf.gz"

def test_two_samples():
    input_id = ["mysample1", "mysample2"]
    input_r1 = [
        input_dir / "tinygex_S1_L001_R1_001.fastq.gz",
        input_dir / "tinygex_S1_L002_R1_001.fastq.gz"
    ]
    input_r2 = [
        input_dir / "tinygex_S1_L001_R2_001.fastq.gz",
        input_dir / "tinygex_S1_L002_R2_001.fastq.gz"
    ]
    output = Path("test_output")

    cmd_pars = [
        meta["executable"],
        "--input_id", ';'.join(input_id),
        "--input_r1", ';'.join([str(r1) for r1 in input_r1]),
        "--input_r2", ';'.join([str(r2) for r2 in input_r2]),
        "--reference_index", reference_index,
        "--reference_gtf", reference_gtf,
        "--output", output,
        "---cpus", "8",
        "--outSAMattributes", "NH;HI;NM;MD"
    ]
    subprocess.run([str(x) for x in cmd_pars], check=True)

    expected_files = [
        "Log.final.out",
        "Aligned.out.bam",
        "Aligned.sorted.out.bam",
        "htseq-count.txt"
    ]
    for iid in input_id:
        for expected_file in expected_files:
            path = output / "per" / iid / expected_file
            assert path.exists(), f"Required file '{path}' is missing"

def test_one_sample():
    input_id = ["mysample", "mysample"]
    input_r1 = [
        input_dir / "tinygex_S1_L001_R1_001.fastq.gz",
        input_dir / "tinygex_S1_L002_R1_001.fastq.gz"
    ]
    input_r2 = [
        input_dir / "tinygex_S1_L001_R2_001.fastq.gz",
        input_dir / "tinygex_S1_L002_R2_001.fastq.gz"
    ]
    output = Path("test_output")

    cmd_pars = [
        meta["executable"],
        "--input_id", ';'.join(input_id),
        "--input_r1", ';'.join([str(r1) for r1 in input_r1]),
        "--input_r2", ';'.join([str(r2) for r2 in input_r2]),
        "--reference_index", reference_index,
        "--reference_gtf", reference_gtf,
        "--output", output,
        "---cpus", "8"
    ]
    subprocess.run([str(x) for x in cmd_pars], check=True)

    assert (output / "feature_info.tsv").exists()
    expected_files = [
        "Log.final.out",
        "Aligned.out.bam",
        "Aligned.sorted.out.bam",
        "htseq-count.txt"
    ]
    for iid in input_id:
        for expected_file in expected_files:
            path = output / "per" / iid / expected_file
            assert path.exists(), f"Required file '{path}' is missing"

if __name__ == '__main__':
    test_two_samples()
    test_one_sample()
