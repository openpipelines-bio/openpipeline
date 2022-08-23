import os
import re
import subprocess
import tempfile
import logging
from sys import stdout

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

## VIASH START
par = {
  'input': ['resources_test/bd_rhapsody_wta_test/raw/12SMK_S1_L432_R1_001.fastq.gz',
            'resources_test/bd_rhapsody_wta_test/raw/12SMK_S1_L432_R2_001.fastq.gz'],
  'output': 'output_dir',
  'subsample': None,
  'reference_genome': 'resources_test/bd_rhapsody_wta_test/raw/GRCh38_primary_assembly_genome_chr1.tar.gz',
  'transcriptome_annotation': 'resources_test/bd_rhapsody_wta_test/raw/gencode_v40_annotation_chr1.gtf',
  'exact_cell_count': None,
  'disable_putative_calling': False,
  'parallel': True,
  'timestamps': False,
  'abseq_reference': [],
  'supplemental_reference': []
}
meta = {
  'resources_dir': 'src/mapping/bd_rhapsody_wta',
  'temp_dir': os.getenv("VIASH_TEMP")
}
## VIASH END

if re.match("[^A-Za-z0-9]", par["sample_prefix"]):
  logger.warning("--sample_prefix should only consist of letters, numbers or hyphens. Replacing all '[^A-Za-z0-9]' with '-'.")
  par["sample_prefix"] = re.sub("[^A-Za-z0-9\\-]", "-", par["sample_prefix"])

def strip_margin(text):
  return re.sub('\n[ \t]*\|', '\n', text)

# if par_input is a directory, look for fastq files
if len(par["input"]) == 1 and os.path.isdir(par["input"][0]):
  par["input"] = [ os.path.join(dp, f) for dp, dn, filenames in os.walk(par["input"]) for f in filenames if re.match(r'.*\.fastq.gz', f) ]

# use absolute paths
par["input"] = [ os.path.abspath(f) for f in par["input"] ]
if par["reference"]:
  par["reference"] = os.path.abspath(par["reference"])
if par["abseq_reference"]:
  par["abseq_reference"] = [ os.path.abspath(f) for f in par["abseq_reference"] ]
par["output"] = os.path.abspath(par["output"])

# create output dir if not exists
if not os.path.exists(par["output"]):
  os.makedirs(par["output"])

# Create params file
config_file = os.path.join(par["output"], "config.yml")
endl = "\n"

content_list = [f"""#!/usr/bin/env cwl-runner

cwl:tool: rhapsody

# This is a YML file used to specify the inputs for a BD Genomics Targeted Rhapsody Analysis pipeline run. See the
# BD Genomics Analysis Setup User Guide (Doc ID: 47383) for more details.

## Reads (required) - Path to your read files in the FASTQ.GZ format. You may specify as many R1/R2 read pairs as you want.
Reads:
"""]

for file in par["input"]:
  content_list.append(strip_margin(f"""\
    | - class: File
    |   location: "{file}"
    |"""))

if par["reference"]:
  content_list.append(strip_margin(f"""\
    |
    |## Reference (optional) - Path to mRNA reference file for pre-designed, supplemental, or custom panel, in FASTA format.
    |Reference:
    | - class: File
    |   location: "{par["reference"]}"
    |"""))

if par["abseq_reference"]:
  content_list.append(strip_margin(f"""\
    |
    |## AbSeq_Reference (optional) - Path to the AbSeq reference file in FASTA format.  Only needed if BD AbSeq Ab-Oligos are used.
    |## For putative cell calling using an AbSeq dataset, please provide an AbSeq reference fasta file as the AbSeq_Reference.
    |AbSeq_Reference:
    |"""))
  for file in par["abseq_reference"]:
    content_list.append(strip_margin(f"""\
      | - class: File
      |   location: {file}
      |"""))

## Putative Cell Calling Settings
content_list.append(strip_margin(f"""\
  |
  |####################################
  |## Putative Cell Calling Settings ##
  |####################################
  |"""))

if par["putative_cell_call"]:
  content_list.append(strip_margin(f"""\
    |## Putative cell calling dataset (optional) - Specify the dataset to be used for putative cell calling: mRNA or AbSeq_Experimental.
    |## For putative cell calling using an AbSeq dataset, please provide an AbSeq_Reference fasta file above.
    |## By default, the mRNA data will be used for putative cell calling.
    |Putative_Cell_Call: {par["putative_cell_call"]}
    |"""))

if par["exact_cell_count"]:
  content_list.append(strip_margin(f"""\
    |## Exact cell count (optional) - Set a specific number (>=1) of cells as putative, based on those with the highest error-corrected read count
    |Exact_Cell_Count: {par["exact_cell_count"]}
    |"""))

if par["disable_putative_calling"]:
  content_list.append(strip_margin(f"""\
    |## Disable Refined Putative Cell Calling (optional) - Determine putative cells using only the basic algorithm (minimum second derivative along the cumulative reads curve).  The refined algorithm attempts to remove false positives and recover false negatives, but may not be ideal for certain complex mixtures of cell types.  Does not apply if Exact Cell Count is set.
    |## The values can be true or false. By default, the refined algorithm is used.
    |Basic_Algo_Only: {str(par["disable_putative_calling"]).lower()}
    |"""))

## Subsample Settings
content_list.append(strip_margin(f"""\
  |
  |########################
  |## Subsample Settings ##
  |########################
  |"""
))

if par["subsample"]:
  content_list.append(strip_margin(f"""\
    |## Subsample (optional) - A number >1 or fraction (0 < n < 1) to indicate the number or percentage of reads to subsample.
    |Subsample: {par["subsample"]}
    |"""))

if par["subsample_seed"]:
  content_list.append(strip_margin(f"""\
    |## Subsample seed (optional) - A seed for replicating a previous subsampled run.
    |Subsample_seed: {par["subsample_seed"]}
    |"""))


## Multiplex options
content_list.append(strip_margin(f"""\
  |
  |#######################
  |## Multiplex options ##
  |#######################
  |"""
))

if par["sample_tags_version"]:
  content_list.append(strip_margin(f"""\
    |## Sample Tags Version (optional) - Specify if multiplexed run: human, hs, mouse or mm
    |Sample_Tags_Version: {par["sample_tags_version"]}
    |"""))

if par["tag_names"]:
  content_list.append(strip_margin(f"""\
    |## Tag_Names (optional) - Specify the tag number followed by '-' and the desired sample name to appear in Sample_Tag_Metrics.csv
    |# Do not use the special characters: &, (), [], {{}},  <>, ?, |
    |Tag_Names: [{', '.join(par["tag_names"])}]
    |"""))

## VDJ options
content_list.append(strip_margin(f"""\
  |
  |#################
  |## VDJ options ##
  |#################
  |"""
))

if par["vdj_version"]:
  content_list.append(strip_margin(f"""\
    |## VDJ Version (optional) - Specify if VDJ run: human, mouse, humanBCR, humanTCR, mouseBCR, mouseTCR
    |VDJ_Version: {par["vdj_version"]}
    |"""))

## VDJ options
content_list.append(strip_margin(f"""\
  |
  |########################
  |## Additional Options ##
  |########################
  |"""
))

if par["sample_prefix"]:
  content_list.append(strip_margin(f"""\
    |## Run Name (optional) -  Specify a run name to use as the output file base name. Use only letters, numbers, or hyphens. Do not use special characters or spaces.
    |Run_Name: {par["sample_prefix"]}
    |"""))

## Write config to file
config_content = ''.join(content_list)

with open(config_file, "w") as f:
  f.write(config_content)

## Process parameters
proc_pars = ["--no-container", "--outdir", par["output"]]

if par["parallel"]:
  proc_pars.append("--parallel")

if par["timestamps"]:
  proc_pars.append("--timestamps")

# create cwl file (if need be)
orig_cwl_file=os.path.join(meta["resources_dir"], "rhapsody_targeted_1.10.1_nodocker.cwl")
if par["override_min_ram"] or par["override_min_cores"]:
  cwl_file = os.path.join(par["output"], "pipeline.cwl")

  # Read in the file
  with open(orig_cwl_file, 'r') as file :
    cwl_data = file.read()

  # Replace the target string
  if par["override_min_ram"]:
    cwl_data = re.sub('"ramMin": [^\n]*,\n', f'"ramMin": {par["override_min_ram"] * 1000},\n', cwl_data)
  if par["override_min_cores"]:
    cwl_data = re.sub('"coresMin": [^\n]*,\n', f'"coresMin": {par["override_min_cores"]},\n', cwl_data)

  # Write the file out again
  with open(cwl_file, 'w') as file:
    file.write(cwl_data)
else:
  cwl_file = orig_cwl_file

## Run pipeline
if not par["dryrun"]:
  with tempfile.TemporaryDirectory(prefix="cwl-bd_rhapsody_wta-", dir=meta["temp_dir"]) as temp_dir:
    cmd = ["cwl-runner"] + proc_pars + [cwl_file, os.path.basename(config_file)]

    env = dict(os.environ)
    env["TMPDIR"] = temp_dir

    logger.info("> " + ' '.join(cmd))

    p = subprocess.check_call(
      cmd,
      cwd=os.path.dirname(config_file),
      env=env
    )

  # # look for counts file
  # if not par["sample_prefix"]:
  #   par["sample_prefix"] = "sample"
  # counts_filename = par["sample_prefix"] + "_RSEC_MolsPerCell.csv"
  
  # if par["sample_tags_version"]:
  #   counts_filename = "Combined_" + counts_filename
  # counts_file = os.path.join(par["output"], counts_filename)
  
  # if not os.path.exists(counts_file):
  #   raise ValueError(f"Could not find output counts file '{counts_filename}'")

  # # look for metrics file
  # metrics_filename = par["sample_prefix"] + "_Metrics_Summary.csv"
  # metrics_file = os.path.join(par["output"], metrics_filename)
  # if not os.path.exists(metrics_file):
  #   raise ValueError(f"Could not find output metrics file '{metrics_filename}'")

