from unittest import TestCase, main
import subprocess
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory
from textwrap import dedent

## VIASH START
meta = {
    'executable': './target/docker/mapping/cellranger_multi/cellranger_multi',
    'resources_dir': 'resources_test/',
    'cpus': 15,
    'memory_gb': 20
}
## VIASH END

input1_R1 = meta["resources_dir"] + "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R1_001.fastq.gz"
input1_R2 = meta["resources_dir"] + "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R2_001.fastq.gz"
input2_R1 = meta["resources_dir"] + "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R1_001.fastq.gz"
input2_R2 = meta["resources_dir"] + "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R2_001.fastq.gz"
input3_R1 =  meta["resources_dir"] + "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R1_001.fastq.gz"
input3_R2 =  meta["resources_dir"] + "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R2_001.fastq.gz"
gex_reference =  meta["resources_dir"] +  "reference_gencodev41_chr1/reference_cellranger.tar.gz"
feature_reference = meta["resources_dir"] + "10x_5k_anticmv/raw/feature_reference.csv"
vdj_reference = meta["resources_dir"]+  "10x_5k_anticmv/raw/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz"


class TestCellrangerMulti(TestCase):
    def _run_and_check_output(self, args_as_list, expected_raise=False):
        try:
            subprocess.check_output([meta['executable']] + args_as_list, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            if not expected_raise:
                print(e.stdout.decode("utf-8"))
            raise e

    def test_cellranger_multi(self):
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
        if meta['cpus']:
            args.extend(["---cpus", str(meta['cpus'])])
        if meta['memory_gb']:
            args.extend(["---memory", f"{meta['memory_gb']}GB"])
        self._run_and_check_output(args)

        # check for raw data
        self.assertTrue(Path("output1/multi/count/raw_feature_bc_matrix.h5").is_file())

        # check for metrics summary
        self.assertTrue(Path("output1/per_sample_outs/run/metrics_summary.csv").is_file())

        # check for filtered gex+ab data
        self.assertTrue(Path("output1/per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file())

        # check for vdj data
        self.assertTrue(Path("output1/per_sample_outs/run/vdj_t/filtered_contig_annotations.csv").is_file())

    def test_cellranger_multi_decompressed_reference(self):
        with tarfile.open(gex_reference) as open_tarfile:
            with TemporaryDirectory() as tmpdirname:
                open_tarfile.extractall(tmpdirname)
                self._run_and_check_output([
                        "--output", "output2/",
                        "--input", input1_R1,
                        "--input", input1_R2,
                        "--input", input2_R1,
                        "--input", input2_R2,
                        "--input", input3_R1,
                        "--input", input3_R2,
                        "--gex_reference", tmpdirname,
                        "--vdj_reference", vdj_reference,
                        "--feature_reference", feature_reference,
                        "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset;5k_human_antiCMV_T_TBNK_connect_VDJ_subset",
                        "--library_type", "Gene Expression;Antibody Capture;VDJ",
                        "--dryrun"])
    
    def test_cellranger_multi_directory_input(self):
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
        if meta['cpus']:
            args.extend(["---cpus", str(meta['cpus'])])
        if meta['memory_gb']:
            args.extend(["---memory", f"{meta['memory_gb']}GB"])
        self._run_and_check_output(args)

    def test_cellranger_multi_applies_gex_options(self):
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
        if meta['cpus']:
            args.extend(["---cpus", str(meta['cpus'])])
        if meta['memory_gb']:
            args.extend(["---memory", f"{meta['memory_gb']}GB"])
        self._run_and_check_output(args)

        self.assertTrue(Path("output3/config.csv").is_file())
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

    def test_cellranger_multi_no_vdj_reference(self):
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
        if meta['cpus']:
            args.extend(["---cpus", str(meta['cpus'])])
        if meta['memory_gb']:
            args.extend(["---memory", f"{meta['memory_gb']}GB"])
        self._run_and_check_output(args)
        self.assertTrue(Path("output4/config.csv").is_file())

if __name__ == "__main__":
    main()