import subprocess
from pathlib import Path
import mudata as md

## VIASH START
meta = {
  "functionality_name": "bd_rhapsody2",
  "executable": "target/docker/mapping/bd_rhapsody2/bd_rhapsody2",
  "resources_dir": "src/mapping/bd_rhapsody2",
  "cpus": 8,
  "memory_mb": 4096,
}

reference_small_gtf = "resources_test/reference_gencodev41_chr1/reference.gtf"
reference_small_fa = "resources_test/reference_gencodev41_chr1/reference.fa"
bdabseq_panel_fa = "resources_test/bdrhap_5kjrt/raw/BDAbSeq_ImmuneDiscoveryPanel.fasta"
reference_file = "resources_test/reference_gencodev41_chr1/Rhap_reference.tar.gz"
abc_reads = "resources_test/bdrhap_5kjrt/raw/12ABC_S1_L432_R1_001_subset.fastq.gz;resources_test/bdrhap_5kjrt/raw/12ABC_S1_L432_R2_001_subset.fastq.gz"
wta_reads = "resources_test/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz;resources_test/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz"
## VIASH END

# wta_reads = f"{meta["resources_dir"]}/12WTA_S1_L432_R1_001_subset.fastq.gz;{meta["resources_dir"]}/12WTA_S1_L432_R2_001_subset.fastq.gz"
# abc_reads = f"{meta["resources_dir"]}/raw/12ABC_S1_L432_R1_001_subset.fastq.gz;{meta["resources_dir"]}/raw/12ABC_S1_L432_R2_001_subset.fastq.gz"
# reference_file = f"{meta["resources_dir"]}/Rhap_reference_small.tar.gz"
# bdabseq_panel_fa = f"{meta["resources_dir"]}/raw/BDAbSeq_ImmuneDiscoveryPanel.fasta"

# Run executable
print(f">> Run {meta['functionality_name']}", flush=True)
output_dir = Path("output")
subprocess.run([
  meta['executable'],
  "--reads=resources_test/bdrhap_5kjrt/raw/12ABC_S1_L432_R1_001_subset.fastq.gz",
  "--reads=resources_test/bdrhap_5kjrt/raw/12ABC_S1_L432_R2_001_subset.fastq.gz",
  "--reads=resources_test/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz",
  "--reads=resources_test/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz",
  f"--reference_archive=resources_test/reference_gencodev41_chr1/Rhap_reference.tar.gz",
  f"--abseq_reference=resources_test/bdrhap_5kjrt/raw/BDAbSeq_ImmuneDiscoveryPanel.fasta",
  "--output_dir=output",
  "--exact_cell_count=4900",
  f"---cpus=2",
  f"---memory=10gb",
  "--output_seurat=seurat.rds",
  "--output_mudata=mudata.h5mu",
  "--metrics_summary=metrics_summary.csv",
  "--pipeline_report=pipeline_report.html",
])


# Check if output exists
print(">> Check if output exists", flush=True)
assert (output_dir / "sample_Bioproduct_Stats.csv").exists()
assert (output_dir / "sample_Metrics_Summary.csv").exists()
assert (output_dir / "sample_Pipeline_Report.html").exists()
assert (output_dir / "sample_RSEC_MolsPerCell_MEX.zip").exists()
assert (output_dir / "sample_RSEC_MolsPerCell_Unfiltered_MEX.zip").exists()
assert (output_dir / "sample_Seurat.rds").exists()
assert (output_dir / "sample.h5mu").exists()

# check individual outputs
assert Path("seurat.rds").exists()
assert Path("mudata.h5mu").exists()
assert Path("metrics_summary.csv").exists()
assert Path("pipeline_report.html").exists()

print(">> Check contents of output", flush=True)
data = md.read_h5mu("mudata.h5mu")

assert data.n_obs == 4900, "Number of cells is incorrect"
assert "rna" in data.mod, "RNA data is missing"
assert "prot" in data.mod, "RNA data is missing"


data_rna = data.mod["rna"]
assert data_rna.n_vars == 2231, "Number of genes is incorrect"

# row sum should be greater than 950
assert data_rna.X.sum(axis=1).min() > 0, "Number of reads per cell is incorrect"
assert data_rna.var.Raw_Reads.sum() == 157808, "Number of reads is incorrect"

# TODO: add VDJ, SMK, ATAC, and targeted RNA to test


#########################################################################################

print("> Test successful", flush=True)
