from pathlib import Path
import tarfile
from textwrap import dedent
import re
import pytest
import sys

## VIASH START
meta = {
    'executable': './target/docker/mapping/cellranger_multi/cellranger_multi',
    'resources_dir': 'resources_test/',
    'cpus': 15,
    'memory_gb': 20,
    'config': 'src/mapping/cellranger_multi/config.vsh.yaml'
}
## VIASH END

input1_R1 = meta["resources_dir"] + "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R1_001.fastq.gz"
input1_R2 = meta["resources_dir"] + "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R2_001.fastq.gz"
input2_R1 = meta["resources_dir"] + "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R1_001.fastq.gz"
input2_R2 = meta["resources_dir"] + "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R2_001.fastq.gz"
input3_R1 = meta["resources_dir"] + "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R1_001.fastq.gz"
input3_R2 = meta["resources_dir"] + "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R2_001.fastq.gz"
gex_reference = meta["resources_dir"] + "reference_gencodev41_chr1/reference_cellranger.tar.gz"
feature_reference = meta["resources_dir"] + "10x_5k_anticmv/raw/feature_reference.csv"
vdj_reference = meta["resources_dir"] + "10x_5k_anticmv/raw/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz"


def test_cellranger_multi(run_component):
    args = [
            "--output", "output1/",
            "--input", input1_R1,
            "--input", input1_R2,
            "--input", input2_R1,
            "--input", input2_R2,
            "--input", input3_R1,
            "--input", input3_R2,
            "--gex_reference", gex_reference,
            "--vdj_reference", vdj_reference,
            "--feature_reference", feature_reference,
            "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset;5k_human_antiCMV_T_TBNK_connect_VDJ_subset",
            "--library_type", "Gene Expression;Antibody Capture;VDJ"]
    run_component(args)

    # check for raw data
    assert Path("output1/multi/count/raw_feature_bc_matrix.h5").is_file()

    # check for metrics summary
    assert Path("output1/per_sample_outs/run/metrics_summary.csv").is_file()

    # check for filtered gex+ab data
    assert Path("output1/per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file()

    # check for vdj data
    assert Path("output1/per_sample_outs/run/vdj_t/filtered_contig_annotations.csv").is_file()

def test_cellranger_multi_decompressed_reference(run_component, tmp_path):
    extracted_tar = tmp_path / "reference_extracted"
    extracted_tar.mkdir()
    with tarfile.open(gex_reference) as open_tarfile:
        open_tarfile.extractall(extracted_tar)
        run_component([
                "--output", "output2/",
                "--input", input1_R1,
                "--input", input1_R2,
                "--input", input2_R1,
                "--input", input2_R2,
                "--input", input3_R1,
                "--input", input3_R2,
                "--gex_reference", extracted_tar,
                "--vdj_reference", vdj_reference,
                "--feature_reference", feature_reference,
                "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset;5k_human_antiCMV_T_TBNK_connect_VDJ_subset",
                "--library_type", "Gene Expression;Antibody Capture;VDJ",
                "--dryrun"])

def test_cellranger_multi_directory_input(run_component):
    args=[
        "--output", "output5/",
        "--input", meta["resources_dir"] + "10x_5k_anticmv/raw/",
        "--gex_reference", gex_reference,
        "--vdj_reference", vdj_reference,
        "--feature_reference", feature_reference,
        "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset;5k_human_antiCMV_T_TBNK_connect_VDJ_subset",
        "--library_type", "Gene Expression;Antibody Capture;VDJ",
        "--gex_secondary_analysis", "true",
        "--gex_generate_bam", "false",
        "--gex_include_introns", "false",
        "--dryrun"]
    run_component(args)

def test_vdj_inner_enrichment_primers(run_component, tmp_path):
    
    enrichment_primers_file = tmp_path / "enrichment_primers.txt"
    with enrichment_primers_file.open('w') as primers_file_open:
        primers_file_open.write("AGTCTCTCAGCTGGTACACG\nTCTGATGGCTCAAACACAGC")
    args=[
        "--output", "output5/",
        "--input", meta["resources_dir"] + "10x_5k_anticmv/raw/",
        "--gex_reference", gex_reference,
        "--vdj_reference", vdj_reference,
        "--feature_reference", feature_reference,
        "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset;5k_human_antiCMV_T_TBNK_connect_VDJ_subset",
        "--library_type", "Gene Expression;Antibody Capture;VDJ",
        "--gex_secondary_analysis", "true",
        "--gex_generate_bam", "false",
        "--gex_include_introns", "false",
        "--vdj_inner_enrichment_primers", str(enrichment_primers_file),
        "--dryrun"]
    run_component(args)
    assert Path("output5/config.csv").is_file()
    with Path("output5/config.csv").open('r') as config_file:
        config_contents = config_file.read()
    expected_csv_content = r"\[vdj\]\nreference,.*\ninner-enrichment-primers,.*\n"
    assert re.search(expected_csv_content, config_contents)

def test_cellranger_multi_applies_gex_options(run_component):
    args=[
            "--output", "output3/",
            "--input", input1_R1,
            "--input", input1_R2,
            "--input", input2_R1,
            "--input", input2_R2,
            "--input", input3_R1,
            "--input", input3_R2,
            "--gex_reference", gex_reference,
            "--vdj_reference", vdj_reference,
            "--feature_reference", feature_reference,
            "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset;5k_human_antiCMV_T_TBNK_connect_VDJ_subset",
            "--library_type", "Gene Expression;Antibody Capture;VDJ",
            "--gex_secondary_analysis", "true",
            "--gex_generate_bam", "false",
            "--gex_include_introns", "false",
            "--dryrun"]
    run_component(args)

    assert Path("output3/config.csv").is_file()
    with Path("output3/config.csv").open('r') as config_file:
        config_contents = config_file.read()
    expected_csv_content = dedent(
        """\
        chemistry,auto
        no-secondary,False
        no-bam,True
        include-introns,False
        """)
    print (expected_csv_content, flush=True)
    assert expected_csv_content in config_contents

def test_cellranger_multi_no_vdj_reference(run_component):
    # GH291
    args=[
        "--output", "output4/",
        "--input", input1_R1,
        "--input", input1_R2,
        "--input", input2_R1,
        "--input", input2_R2,
        "--input", input3_R1,
        "--input", input3_R2,
        "--gex_reference", gex_reference,
        "--feature_reference", feature_reference,
        "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset",
        "--library_type", "Gene Expression;Antibody Capture",
        "--dryrun"]
    run_component(args)
    assert Path("output4/config.csv").is_file()

def test_cellranger_multi_crispt_data(run_component):
    args = [
        "--input", meta["resources_dir"] + "10x_5k_lung_crispr/raw/SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_subset_S5_L001_R1_001.fastq.gz",
        "--input", meta["resources_dir"] + "10x_5k_lung_crispr/raw/SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_subset_S5_L001_R2_001.fastq.gz",
        "--input", meta["resources_dir"] + "10x_5k_lung_crispr/raw/SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_subset_S4_L001_R1_001.fastq.gz",
        "--input", meta["resources_dir"] + "10x_5k_lung_crispr/raw/SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_subset_S4_L001_R2_001.fastq.gz",
        "--library_id", "SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_subset;SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_subset",
        "--library_type", "Gene Expression;CRISPR Guide Capture",
        "--gex_reference", gex_reference,
        "--feature_reference", meta["resources_dir"] + "10x_5k_lung_crispr/raw/SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_count_feature_reference_corrected.csv",
        "--output", "output5/"
    ]
    run_component(args)
    # check for raw data
    assert Path("output5/multi/count/raw_feature_bc_matrix.h5").is_file()
    # check for metrics summary
    assert Path("output5/per_sample_outs/run/metrics_summary.csv").is_file()
    # check for filtered gex+ab data
    assert Path("output5/per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file()
    # check for crispr data
    assert Path("output5/per_sample_outs/run/count/crispr_analysis/").is_dir()

def test_cellranger_multi_helper_input(run_component):
    args = [
            "--output", "output6/",
            "--gex_input", input1_R1,
            "--gex_input", input1_R2,
            "--abc_input", input2_R1,
            "--abc_input", input2_R2,
            "--vdj_input", input3_R1,
            "--vdj_input", input3_R2,
            "--gex_reference", gex_reference,
            "--vdj_reference", vdj_reference,
            "--feature_reference", feature_reference]
    run_component(args)

    # check for raw data
    assert Path("output6/multi/count/raw_feature_bc_matrix.h5").is_file()

    # check for metrics summary
    assert Path("output6/per_sample_outs/run/metrics_summary.csv").is_file()

    # check for filtered gex+ab data
    assert Path("output6/per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file()

    # check for vdj data
    assert Path("output6/per_sample_outs/run/vdj_t/filtered_contig_annotations.csv").is_file()

def test_cellranger_multi_combined_helper_and_global_input(run_component):
    args = [
            "--output", "output7/",
            "--input", input1_R1,
            "--input", input1_R2,
            "--abc_input", input2_R1,
            "--abc_input", input2_R2,
            "--vdj_input", input3_R1,
            "--vdj_input", input3_R2,
            "--gex_reference", gex_reference,
            "--vdj_reference", vdj_reference,
            "--feature_reference", feature_reference,
            "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset",
            "--library_type", "Gene Expression"]
    run_component(args)

    # check for raw data
    assert Path("output7/multi/count/raw_feature_bc_matrix.h5").is_file()

    # check for metrics summary
    assert Path("output7/per_sample_outs/run/metrics_summary.csv").is_file()

    # check for filtered gex+ab data
    assert Path("output7/per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file()

    # check for vdj data
    assert Path("output7/per_sample_outs/run/vdj_t/filtered_contig_annotations.csv").is_file()

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
