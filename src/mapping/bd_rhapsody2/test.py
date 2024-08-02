import subprocess
import gzip
from pathlib import Path
from typing import Tuple
import numpy as np
import random
import mudata as md

## VIASH START
meta = {
  "name": "bd_rhapsody2",
  "executable": "target/docker/mapping/bd_rhapsody_2/bd_rhapsody_2",
  "resources_dir": "src/mapping/bd_rhapsody2",
  "cpus": 8,
  "memory_mb": 4096,
}
## VIASH END

import sys
sys.path.append(meta["resources_dir"])

from rhapsody_cell_label import index_to_sequence

meta["executable"] = Path(meta["executable"])
meta["resources_dir"] = Path(meta["resources_dir"])

# #########################################################################################

# Generate samples
sample_path = Path(meta["resources_dir"]) / "BDAbSeq_ImmuneDiscoveryPanel.fasta"

sample_content = f"""\
cat > {meta["resources_dir"]}/BDAbSeq_ImmuneDiscoveryPanel.fasta <<EOF
CD11c:B-LY6|ITGAX|AHS0056|pAbO Catalog_940024
ATGCGTTGCGAGAGATATGCGTAGGTTGCTGATTGG
>CD14:MPHIP9|CD14|AHS0037|pAbO Catalog_940005
TGGCCCGTGGTAGCGCAATGTGAGATCGTAATAAGT
>CXCR5|CXCR5|AHS0039|pAbO Catalog_940042
AGGAAGGTCGATTGTATAACGCGGCATTGTAACGGC
>CD19:SJ25C1|CD19|AHS0030|pAbO Catalog_940004
TAGTAATGTGTTCGTAGCCGGTAATAATCTTCGTGG
>CD25:2A3|IL2RA|AHS0026|pAbO Catalog_940009
AGTTGTATGGGTTAGCCGAGAGTAGTGCGTATGATT
>CD27:M-T271|CD27|AHS0025|pAbO Catalog_940018
TGTCCGGTTTAGCGAATTGGGTTGAGTCACGTAGGT
>CD278|ICOS|AHS0012|pAbO Catalog_940043
ATAGTCCGCCGTAATCGTTGTGTCGCTGAAAGGGTT
>CD279:EH12-1|PDCD1|AHS0014|pAbO Catalog_940015
ATGGTAGTATCACGACGTAGTAGGGTAATTGGCAGT
>CD3:UCHT1|CD3E|AHS0231|pAbO Catalog_940307
AGCTAGGTGTTATCGGCAAGTTGTACGGTGAAGTCG
>GITR|TNFRSF18|AHS0104|pAbO Catalog_940096
TCTGTGTGTCGGGTTGAATCGTAGTGAGTTAGCGTG
>Tim3|HAVCR2|AHS0016|pAbO Catalog_940066
TAGGTAGTAGTCCCGTATATCCGATCCGTGTTGTTT
>CD4:SK3|CD4|AHS0032|pAbO Catalog_940001
TCGGTGTTATGAGTAGGTCGTCGTGCGGTTTGATGT
>CD45RA:HI100|PTPRC|AHS0009|pAbO Catalog_940011
AAGCGATTGCGAAGGGTTAGTCAGTACGTTATGTTG
>CD56:NCAM16.2|NCAM1|AHS0019|pAbO Catalog_940007
AGAGGTTGAGTCGTAATAATAATCGGAAGGCGTTGG
>CD62L:DREG-56|SELL|AHS0049|pAbO Catalog_940041
ATGGTAAATATGGGCGAATGCGGGTTGTGCTAAAGT
>CCR7|CCR7|AHS0273|pAbO Catalog_940394
AATGTGTGATCGGCAAAGGGTTCTCGGGTTAATATG
>CXCR6|CXCR6|AHS0148|pAbO Catalog_940234
GTGGTTGGTTATTCGGACGGTTCTATTGTGAGCGCT
>CD127|IL7R|AHS0028|pAbO Catalog_940012
AGTTATTAGGCTCGTAGGTATGTTTAGGTTATCGCG
>CD134:ACT35|TNFRSF4|AHS0013|pAbO Catalog_940060
GGTGTTGGTAAGACGGACGGAGTAGATATTCGAGGT
>CD28:L293|CD28|AHS0138|pAbO Catalog_940226
TTGTTGAGGATACGATGAAGCGGTTTAAGGGTGTGG
>CD272|BTLA|AHS0052|pAbO Catalog_940105
GTAGGTTGATAGTCGGCGATAGTGCGGTTGAAAGCT
>CD8:SK1|CD8A|AHS0228|pAbO Catalog_940305
AGGACATAGAGTAGGACGAGGTAGGCTTAAATTGCT
>HLA-DR|CD74|AHS0035|pAbO Catalog_940010
TGTTGGTTATTCGTTAGTGCATCCGTTTGGGCGTGG
>CD16:3G8|FCGR3A|AHS0053|pAbO Catalog_940006
TAAATCTAATCGCGGTAACATAACGGTGGGTAAGGT
>CD183|CXCR3|AHS0031|pAbO Catalog_940030
AAAGTGTTGGCGTTATGTGTTCGTTAGCGGTGTGGG
>CD196|CCR6|AHS0034|pAbO Catalog_940033
ACGTGTTATGGTGTTGTTCGAATTGTGGTAGTCAGT
>CD137|TNFRSF9|AHS0003|pAbO Catalog_940055
TGACAAGCAACGAGCGATACGAAAGGCGAAATTAGT
>CD161:HP-3G10|KLRB1|AHS0205|pAbO Catalog_940283
TTTAGGACGATTAGTTGTGCGGCATAGGAGGTGTTC
>IgM|IGHM|AHS0198|pAbO Catalog_940276
TTTGGAGGGTAGCTAGTTGCAGTTCGTGGTCGTTTC
>IgD|IGHD|AHS0058|pAbO Catalog_940026
TGAGGGATGTATAGCGAGAATTGCGACCGTAGACTT
EOF
"""

with open(sample_path, 'w') as file:
  file.write(sample_content)

# Generate index
print("> Generate index", flush=True)
cwl_file = meta["resources_dir"] / "bd_rhapsody_make_reference.cwl"
reference_small_gtf = meta["resources_dir"] / "reference.gtf"
reference_small_fa = meta["resources_dir"] / "reference.fa"
bdabseq_panel_fa = meta["resources_dir"] / "BDAbSeq_ImmuneDiscoveryPanel.fasta"

config_file = Path("reference_config.yml")
reference_file = Path("Rhap_reference.tar.gz")


subprocess.run([
    "cwl-runner",
    "--no-container",
    "--preserve-entire-environment",
    "--outdir",
    ".",
    str(cwl_file),
    "--Genome_fasta",
    str(reference_small_fa),
    "--Gtf", 
    str(reference_small_gtf),
    "--Extra_STAR_params",
    "--genomeSAindexNbases 4"
])

#########################################################################################
# Load reference in memory

subprocess.run(["gunzip", f"{reference_small_gtf}.gz"])
subprocess.run(["gunzip", f"{reference_small_fa}.gz"])

from Bio import SeqIO
import gffutils


# Load FASTA sequence
with open(str(reference_small_fa), "r") as handle:
  reference_fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
with open(str(bdabseq_panel_fa), "r") as handle:
  bdabseq_panel_fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

# create in memory db
reference_gtf_db = gffutils.create_db(
  str(reference_small_gtf),
  dbfn=":memory:",
  force=True,
  keep_order=True,
  merge_strategy="merge",
  sort_attribute_values=True,
  disable_infer_transcripts=True,
  disable_infer_genes=True
)

#############################################
# TODO: move helper functions to separate helper file


def generate_bd_read_metadata(
  instrument_id: str = "A00226",
  run_id: str = "970",
  flowcell_id: str = "H5FGVMXY",
  lane: int = 1,
  tile: int = 1101,
  x: int = 1000,
  y: int = 1000,
  illumina_flag: str = "1:N:0",
  sample_id: str = "CAGAGAGG",
) -> str:
  """
  Generate a FASTQ metadata line for a BD Rhapsody FASTQ file.

  Args:
    instrument_id: The instrument ID.
    run_id: The run ID.
    flowcell_id: The flowcell ID.
    lane: The lane number.
    tile: The tile number. Between 1101 and 1112 in the used example data.
    x: The x-coordinate. Between 1000 and 32967 in the used example data.
    y: The y-coordinate. Between 1000 and 37059 in the used example data.
    illumina_flag: The Illumina flag. Either 1:N:0 or 2:N:0 in the used example data.
    sample_id: The sample ID.
  """
  # format: @A00226:970:H5FGVDMXY:1:1101:2645:1000 2:N:0:CAGAGAGG
  return f"@{instrument_id}:{run_id}:{flowcell_id}:{lane}:{tile}:{x}:{y} {illumina_flag}:{sample_id}"


def generate_bd_wta_transcript(
  transcript_length: int = 42,
) -> str:
  """
  Generate a WTA transcript from a given GTF and FASTA file.
  """

  # Randomly select a gene
  gene = random.choice(list(reference_gtf_db.features_of_type("gene")))

  # Find all exons within the gene
  exons = list(reference_gtf_db.children(gene, featuretype="exon", order_by="start"))

  # Calculate total exon length
  total_exon_length = sum(exon.end - exon.start + 1 for exon in exons)

  # If total exon length is less than desired transcript length, use it as is
  max_transcript_length = min(total_exon_length, transcript_length)

  # Build the WTA transcript sequence
  sequence = ""
  for exon in exons:
    exon_seq = str(reference_fasta_dict[exon.seqid].seq[exon.start - 1 : exon.end])  
    sequence += exon_seq

    # Break if desired length is reached
    if len(sequence) >= max_transcript_length:
      sequence = sequence[:max_transcript_length]
      break
  
  # add padding if need be
  if len(sequence) < max_transcript_length:
    sequence += "N" * (max_transcript_length - len(sequence))

  return sequence


def generate_bd_wta_read(
  cell_index: int = 0,
  bead_version: str = "EnhV2",
  umi_length: int = 14,
  transcript_length: int = 42,
) -> Tuple[str, str]:
  """
  Generate a BD Rhapsody WTA read pair for a given cell index.

  Args:
    cell_index: The cell index to generate reads for.
    bead_version: The bead version to use for generating the cell label.
    umi_length: The length of the UMI to generate.
    transcript_length: The length of the transcript to generate

  Returns:
    A tuple of two strings, the first string being the R1 read and the second string being the R2 read.

  More info:
    
    See structure of reads:
    - https://bd-rhapsody-bioinfo-docs.genomics.bd.com/steps/top_steps.html
    - https://bd-rhapsody-bioinfo-docs.genomics.bd.com/steps/steps_cell_label.html
    - https://scomix.bd.com/hc/en-us/articles/360057714812-All-FAQ
    R1 is Cell Label + UMI + PolyT -> 60 bp
      actually, CLS1 + "GTGA" + CLS2 + "GACA" + CLS3 + UMI
    R2 is the actual read -> 42 bp

    Example R1
    CLS1       Link CLS2      Link CLS3       UMI
    AAAATCCTGT GTGA AACCAAAGT GACA GATAGAGGAG CGCATGTTTATAAC
  """
  
  # generate metadata
  per_row = np.floor((32967 - 1000) / 9)
  per_col = np.floor((37059 - 1000) / 9)
  
  assert cell_index >= 0 and cell_index < per_row * per_col, f"cell_index must be between 0 and {per_row} * {per_col}"
  x = 1000 + (cell_index % per_row) * 9
  y = 1000 + (cell_index // per_row) * 9
  meta_r1 = generate_bd_read_metadata(x=x, y=y, illumina_flag="1:N:0")
  meta_r2 = generate_bd_read_metadata(x=x, y=y, illumina_flag="2:N:0")

  # generate r1 (cls1 + link + cls2 + link + cls3 + umi)
  assert cell_index >= 0 and cell_index < 384 * 384 * 384
  cell_label = index_to_sequence(cell_index + 1, bead_version=bead_version)
  # sample random umi
  umi = "".join(random.choices("ACGT", k=umi_length))
  quality_r1 = "I" * (len(cell_label) + len(umi))
  r1 = f"{meta_r1}\n{cell_label}{umi}\n+\n{quality_r1}\n"

  # generate r2 by extracting sequence from fasta and gtf
  wta_transcript = generate_bd_wta_transcript(transcript_length=transcript_length)
  quality_r2 = "I" * transcript_length
  r2 = f"{meta_r2}\n{wta_transcript}\n+\n{quality_r2}\n"

  return r1, r2

def generate_bd_wta_fastq_files(
  num_cells: int = 100,
  num_reads_per_cell: int = 1000,
) -> Tuple[str, str]:
  """
  Generate BD Rhapsody WTA FASTQ files for a given number of cells and transcripts per cell.

  Args:
    num_cells: The number of cells to generate
    num_reads_per_cell: The number of reads to generate per cell

  Returns:
    A tuple of two strings, the first string being the R1 reads and the second string being the R2 reads.
  """
  r1_reads = ""
  r2_reads = ""
  for cell_index in range(num_cells):
    for _ in range(num_reads_per_cell):
      r1, r2 = generate_bd_wta_read(cell_index)
      r1_reads += r1
      r2_reads += r2

  return r1_reads, r2_reads

def generate_bd_abc_read(
  cell_index: int = 0,
  bead_version: str = "EnhV2",
  umi_length: int = 14,
  transcript_length: int = 72,
) -> Tuple[str, str]:
  """
  Generate a BD Rhapsody ABC read pair for a given cell index.

  Args:
    cell_index: The cell index to generate reads for.
    bead_version: The bead version to use for generating the cell label.
    umi_length: The length of the UMI to generate.
    transcript_length: The length of the transcript to generate

  Returns:
    A tuple of two strings, the first string being the R1 read and the second string being the R2 read.
  """
  # generate metadata
  per_row = np.floor((32967 - 1000) / 9)
  per_col = np.floor((37059 - 1000) / 9)
  
  assert cell_index >= 0 and cell_index < per_row * per_col, f"cell_index must be between 0 and {per_row} * {per_col}"
  x = 1000 + (cell_index % per_row) * 9
  y = 1000 + (cell_index // per_row) * 9
  meta_r1 = generate_bd_read_metadata(x=x, y=y, illumina_flag="1:N:0")
  meta_r2 = generate_bd_read_metadata(x=x, y=y, illumina_flag="2:N:0")

  # generate r1 (cls1 + link + cls2 + link + cls3 + umi)
  assert cell_index >= 0 and cell_index < 384 * 384 * 384
  cell_label = index_to_sequence(cell_index + 1, bead_version=bead_version)
  # sample random umi
  umi = "".join(random.choices("ACGT", k=umi_length))
  quality_r1 = "I" * (len(cell_label) + len(umi))
  r1 = f"{meta_r1}\n{cell_label}{umi}\n+\n{quality_r1}\n"

  # generate r2 by sampling sequence from bdabseq_panel_fa
  abseq_seq = str(random.choice(list(bdabseq_panel_fasta_dict.values())).seq)
  abc_prefix = "N" #+ "".join(random.choices("ACGT", k=12))
  abc_data = abseq_seq[:transcript_length - len(abc_prefix)]
  abc_suffix = "A" * (transcript_length - len(abc_prefix) - len(abc_data))

  abc_transcript = f"{abc_prefix}{abc_data}{abc_suffix}"

  quality_r2 = "#" + "I" * (len(abc_transcript) - 1)
  r2 = f"{meta_r2}\n{abc_transcript}\n+\n{quality_r2}\n"

  return r1, r2

def generate_bd_abc_fastq_files(
  num_cells: int = 100,
  num_reads_per_cell: int = 1000,
) -> Tuple[str, str]:
  """
  Generate BD Rhapsody ABC FASTQ files for a given number of cells and transcripts per cell.

  Args:
    num_cells: The number of cells to generate
    num_reads_per_cell: The number of reads to generate per cell

  Returns:
    A tuple of two strings, the first string being the R1 reads and the second string being the R2 reads.
  """
  r1_reads = ""
  r2_reads = ""
  for cell_index in range(num_cells):
    for _ in range(num_reads_per_cell):
      r1, r2 = generate_bd_abc_read(cell_index)
      r1_reads += r1
      r2_reads += r2

  return r1_reads, r2_reads


# Prepare WTA, ABC, and SMK test data
print("> Prepare WTA test data", flush=True)
wta_reads_r1_str, wta_reads_r2_str = generate_bd_wta_fastq_files(num_cells=100, num_reads_per_cell=1000)
with gzip.open("WTAreads_R1.fq.gz", "wt") as f:
  f.write(wta_reads_r1_str)
with gzip.open("WTAreads_R2.fq.gz", "wt") as f:
  f.write(wta_reads_r2_str)

print("> Prepare ABC test data", flush=True)
abc_reads_r1_str, abc_reads_r2_str = generate_bd_abc_fastq_files(num_cells=100, num_reads_per_cell=1000)
with gzip.open("ABCreads_R1.fq.gz", "wt") as f:
  f.write(abc_reads_r1_str)
with gzip.open("ABCreads_R2.fq.gz", "wt") as f:
  f.write(abc_reads_r2_str)

#########################################################################################

# Run executable
print(f">> Run {meta['name']}", flush=True)
output_dir = Path("output")
subprocess.run([
  meta['executable'],
  "--reads=WTAreads_R1.fq.gz;WTAreads_R2.fq.gz",
  "--reads=ABCreads_R1.fq.gz;ABCreads_R2.fq.gz",
  f"--reference_archive={reference_file}",
  f"--abseq_reference={bdabseq_panel_fa}",
  "--output_dir=output",
  "--exact_cell_count=100",
  f"---cpus={meta['cpus'] or 1}",
  f"---memory={meta['memory_mb'] or 2048}mb",
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

assert data.n_obs == 100, "Number of cells is incorrect"
assert "rna" in data.mod, "RNA data is missing"


data_rna = data.mod["rna"]
assert data_rna.n_vars == 1, "Number of genes is incorrect"

# row sum should be greater than 950
assert data_rna.X.sum(axis=1).min() > 950, "Number of reads per cell is incorrect"
assert data_rna.var.Raw_Reads.sum() == 100000, "Number of reads is incorrect"

# TODO: add VDJ, SMK, ATAC, and targeted RNA to test


#########################################################################################

print("> Test successful", flush=True)
