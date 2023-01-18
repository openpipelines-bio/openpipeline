import os
import re
import subprocess
import tempfile
import logging
from sys import stdout
from typing import Any
import pandas as pd
import gzip
import shutil

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

## VIASH START
par = {
    'sample_prefix': 'sample',
    'input': ['resources_test/bdrhap_5kjrt/raw/12ABC_S1_L432_R1_001_subset.fastq.gz',
                'resources_test/bdrhap_5kjrt/raw/12ABC_S1_L432_R2_001_subset.fastq.gz'],
    'mode': 'wta',
    'output': 'foo',
    'subsample': None,
    'reference': ['resources_test/reference_gencodev41_chr1/reference_bd_rhapsody.tar.gz'],
    'transcriptome_annotation': 'resources_test/reference_gencodev41_chr1/reference.gtf.gz',
    'exact_cell_count': None,
    'putative_cell_call': None,
    'disable_putative_calling': False,
    'subsample': None,
    'subsample_seed': None,
    'sample_tags_version': None,
    'tag_names': None,
    'vdj_version': None,
    'dryrun': False,
    'parallel': True,
    'timestamps': False,
    'abseq_reference': ["resources_test/bdrhap_5kjrt/raw/BDAbSeq_ImmuneDiscoveryPanel.fasta"],
    'supplemental_reference': []
}
meta = {
    'resources_dir': os.path.abspath('src/mapping/bd_rhapsody'),
    'temp_dir': os.getenv("VIASH_TEMP"),
    'memory_mb': None,
    'n_proc': None
}
## VIASH END

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def strip_margin(text: str) -> str:
    return re.sub('(\n?)[ \t]*\|', '\\1', text)

def process_params(par: dict[str, Any]) -> str:
    # check input parameters
    assert par["input"] is not None, "Pass at least one set of inputs to --input."
    if par["mode"] == "wta":
        assert len(par["reference"]) == 1, "When mode is \"wta\", --reference should be length 1"
        assert par["transcriptome_annotation"] is not None, "When mode is \"wta\", --transcriptome_annotation should be defined"
    elif par["mode"] == "targeted":
        assert par["transcriptome_annotation"] is None, "When mode is \"targeted\", --transcriptome_annotation should be undefined"
        assert par["supplemental_reference"] is None, "When mode is \"targeted\", --supplemental_reference should be undefined"

    # checking sample prefix
    if re.match("[^A-Za-z0-9]", par["sample_prefix"]):
        logger.warning("--sample_prefix should only consist of letters, numbers or hyphens. Replacing all '[^A-Za-z0-9]' with '-'.")
        par["sample_prefix"] = re.sub("[^A-Za-z0-9\\-]", "-", par["sample_prefix"])

    # if par_input is a directory, look for fastq files
    if len(par["input"]) == 1 and os.path.isdir(par["input"][0]):
        par["input"] = [ os.path.join(dp, f) for dp, dn, filenames in os.walk(par["input"]) for f in filenames if re.match(r'.*\.fastq.gz', f) ]

    # use absolute paths
    par["input"] = [ os.path.abspath(f) for f in par["input"] ]
    if par["reference"]:
        par["reference"] = [ os.path.abspath(f) for f in par["reference"] ]
    if par["transcriptome_annotation"]:
        par["transcriptome_annotation"] = os.path.abspath(par["transcriptome_annotation"])
    if par["abseq_reference"]:
        par["abseq_reference"] = [ os.path.abspath(f) for f in par["abseq_reference"] ]
    if par["supplemental_reference"]:
        par["supplemental_reference"] = [ os.path.abspath(f) for f in par["supplemental_reference"] ]
    par["output"] = os.path.abspath(par["output"])
    
    return par

def generate_config(par: dict[str, Any]) -> str:
    content_list = [strip_margin(f"""\
        |#!/usr/bin/env cwl-runner
        |
        |cwl:tool: rhapsody
        |
        |# This is a YML file used to specify the inputs for a BD Genomics {"WTA" if par["mode"] == "wta" else "Targeted" } Rhapsody Analysis pipeline run. See the
        |# BD Genomics Analysis Setup User Guide (Doc ID: 47383) for more details.
        |
        |## Reads (required) - Path to your read files in the FASTQ.GZ format. You may specify as many R1/R2 read pairs as you want.
        |Reads:
        |""")]

    for file in par["input"]:
        content_list.append(strip_margin(f"""\
            | - class: File
            |   location: "{file}"
            |"""))

    if par["reference"] and par["mode"] == "wta":
        content_list.append(strip_margin(f"""\
            |
            |## Reference_Genome (required) - Path to STAR index for tar.gz format. See Doc ID: 47383 for instructions to obtain pre-built STAR index file.
            |Reference_Genome:
            |   class: File
            |   location: "{par["reference"][0]}"
            |"""))

    if par["reference"] and par["mode"] == "targeted":
        content_list.append(strip_margin(f"""\
            |
            |## Reference (optional) - Path to mRNA reference file for pre-designed, supplemental, or custom panel, in FASTA format.
            |Reference:
            |"""))
        for file in par["reference"]:
            content_list.append(strip_margin(f"""\
                | - class: File
                |   location: {file}
                |"""))

    if par["transcriptome_annotation"]:
        content_list.append(strip_margin(f"""\
            |
            |## Transcriptome_Annotation (required) - Path to GTF annotation file
            |Transcriptome_Annotation:
            |   class: File
            |   location: "{par["transcriptome_annotation"]}"
            |"""))

    if par["abseq_reference"]:
        content_list.append(strip_margin(f"""\
            |
            |## AbSeq_Reference (optional) - Path to the AbSeq reference file in FASTA format.  Only needed if BD AbSeq Ab-Oligos are used.
            |AbSeq_Reference:
            |"""))
        for file in par["abseq_reference"]:
            content_list.append(strip_margin(f"""\
                | - class: File
                |   location: {file}
                |"""))

    if par["supplemental_reference"]:
        content_list.append(strip_margin(f"""\
            |
            |## Supplemental_Reference (optional) - Path to the supplemental reference file in FASTA format.  Only needed if there are additional transgene sequences used in the experiment.
            |Supplemental_Reference:
            |"""))
        for file in par["supplemental_reference"]:
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
    return ''.join(content_list)

def generate_cwl_file(par: dict[str, Any], meta: dict[str, Any]) -> str:
        # create cwl file (if need be)
    if par["mode"] == "wta":
        orig_cwl_file=os.path.join(meta["resources_dir"], "rhapsody_wta_1.10.1_nodocker.cwl")
    elif par["mode"] == "targeted":
        orig_cwl_file=os.path.join(meta["resources_dir"], "rhapsody_targeted_1.10.1_nodocker.cwl")

    # Inject computational requirements into pipeline
    if meta["memory_mb"] or meta["cpus"]:
        cwl_file = os.path.join(par["output"], "pipeline.cwl")

        # Read in the file
        with open(orig_cwl_file, 'r') as file :
            cwl_data = file.read()

        # Inject computational requirements into pipeline
        if meta["memory_mb"]:
            memory = int(meta["memory_mb"]) - 2000 # keep 2gb for OS
            cwl_data = re.sub('"ramMin": [^\n]*,\n', f'"ramMin": {memory},\n', cwl_data)
        if meta["cpus"]:
            cwl_data = re.sub('"coresMin": [^\n]*,\n', f'"coresMin": {meta["cpus"]},\n', cwl_data)

        # Write the file out again
        with open(cwl_file, 'w') as file:
            file.write(cwl_data)
    else:
        cwl_file = orig_cwl_file

    return cwl_file

def process_fasta(feature_type: str, path: str) -> pd.DataFrame:
    with open(path) as f:
        df = pd.DataFrame(data={
        'feature_type': feature_type,
        'feature_id': [line[1:].strip() for line in f if line[0] == ">"],
        'reference_file': os.path.basename(path),
        })
        return df

def process_gtf(feature_type: str, path: str) -> pd.DataFrame:
    with open(path) as f:
        data = []
        for line in f:
            if not line.startswith("#"):
                attr = dict(item.strip().split(' ') for item in line.split('\t')[8].strip('\n').split(';') if item)
                row = {
                    'feature_types': feature_type,
                    'feature_ids': attr["gene_name"].strip("\""),
                    'reference_file': os.path.basename(path),
                }
            data.append(row)
    df = pd.DataFrame(data)
    df = df.drop_duplicates()
    return df

def extract_feature_types(par: dict[str, Any]):
    feature_types = []

    if par["mode"] == "targeted":
        for file in par["reference"]:
            logger.info(f"Processing reference fasta {file}")
            feature_types.append(process_fasta("Gene Expression", file))

    if par["mode"] == "wta":
        file = par["transcriptome_annotation"]
        logger.info(f"Processing reference gtf {file}")
        feature_types.append(process_gtf("Gene Expression", file))

    if par["abseq_reference"]:
        for file in par["abseq_reference"]:
            logger.info(f"Processing abseq fasta {file}")
            feature_types.append(process_fasta("Antibody Capture", file))

    if par["supplemental_reference"]:
        for file in par["supplemental_reference"]:
            logger.info(f"Processing supp fasta {file}")
            feature_types.append(process_fasta("Other", file))
    
    return pd.concat(feature_types)

def main(par: dict[str, Any], meta: dict[str, Any]):
    # Preprocess params
    par = process_params(par)

    # Create output dir if not exists
    if not os.path.exists(par["output"]):
        os.makedirs(par["output"])

    ## Process parameters
    proc_pars = ["--no-container", "--outdir", par["output"]]

    if par["parallel"]:
        proc_pars.append("--parallel")

    if par["timestamps"]:
        proc_pars.append("--timestamps")

    with tempfile.TemporaryDirectory(prefix="cwl-bd_rhapsody_wta-", dir=meta["temp_dir"]) as temp_dir:
        # extract transcriptome gtf if need be
        if par["transcriptome_annotation"] and is_gz_file(par["transcriptome_annotation"]):
            with open(os.path.join(temp_dir, "transcriptome.gtf"), 'wb') as genes_uncompressed:
                with gzip.open(par["transcriptome_annotation"], 'rb') as genes_compressed:
                    shutil.copyfileobj(genes_compressed, genes_uncompressed)
                    par["transcriptome_annotation"] = genes_uncompressed.name

        # Create params file
        config_file = os.path.join(par["output"], "config.yml")
        config_content = generate_config(par)
        with open(config_file, "w") as f:
            f.write(config_content)

        # Create cwl file (if need be)
        cwl_file = generate_cwl_file(par, meta)

        ## Run pipeline
        if not par["dryrun"]:
            cmd = ["cwl-runner"] + proc_pars + [cwl_file, os.path.basename(config_file)]

        env = dict(os.environ)
        env["TMPDIR"] = temp_dir

        logger.info("> " + ' '.join(cmd))
        _ = subprocess.check_call(
            cmd,
            cwd=os.path.dirname(config_file),
            env=env
        )

        # extracting feature ids from references
        # extract info from reference files (while they still exist)
        feature_df = extract_feature_types(par)
        feature_types_file = os.path.join(par["output"], "feature_types.tsv")
        feature_df.to_csv(feature_types_file, sep="\t", index=False)

    # look for counts file
    if not par["sample_prefix"]:
        par["sample_prefix"] = "sample"
    counts_filename = par["sample_prefix"] + "_RSEC_MolsPerCell.csv"
    
    if par["sample_tags_version"]:
        counts_filename = "Combined_" + counts_filename
    counts_file = os.path.join(par["output"], counts_filename)
    
    if not os.path.exists(counts_file):
        raise ValueError(f"Could not find output counts file '{counts_filename}'")

    # look for metrics file
    metrics_filename = par["sample_prefix"] + "_Metrics_Summary.csv"
    metrics_file = os.path.join(par["output"], metrics_filename)
    if not os.path.exists(metrics_file):
        raise ValueError(f"Could not find output metrics file '{metrics_filename}'")

if __name__ == "__main__":
    main(par, meta)