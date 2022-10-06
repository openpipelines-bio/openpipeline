from unittest import TestCase, main
import subprocess
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory

## VIASH START
meta = {
    'executable': './target/docker/mapping/cellranger_multi/cellranger_multi',
    'resources_dir': 'resources_test/'
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
            raise e

    def test_cellranger_multi(self):
        self._run_and_check_output([
               "--output", "output/",
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
               "--library_type", "Gene Expression;Antibody Capture;VDJ"])
        # check for raw data
        self.assertTrue(Path("output/multi/count/raw_feature_bc_matrix.h5").is_file())

        # check for metrics summary
        self.assertTrue(Path("output/per_sample_outs/run/metrics_summary.csv").is_file())
        
        # check for filtered gex+ab data
        self.assertTrue(Path("output/per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file())

        # check for vdj data
        self.assertTrue(Path("output/per_sample_outs/run/vdj_t/filtered_contig_annotations.csv").is_file())
    
    def test_cellranger_multi_decompressed_reference(self):
        with tarfile.open(gex_reference) as open_tarfile:
            with TemporaryDirectory() as tmpdirname:
                open_tarfile.extractall(tmpdirname)
                self._run_and_check_output([
                        "--output", "test/",
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

if __name__ == "__main__":
    main()