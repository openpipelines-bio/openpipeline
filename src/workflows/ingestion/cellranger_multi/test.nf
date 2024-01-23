nextflow.enable.dsl=2

include { cellranger_multi } from params.rootDir + "/target/nextflow/workflows/ingestion/cellranger_multi/main.nf"

workflow test_wf {
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList([
      [
        id: "foo",
        input:[resources_test.resolve("10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R1_001.fastq.gz"),
               resources_test.resolve("10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R2_001.fastq.gz"),
               resources_test.resolve("10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R1_001.fastq.gz"),
               resources_test.resolve("10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R2_001.fastq.gz"),
               resources_test.resolve("10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R1_001.fastq.gz"),
               resources_test.resolve("10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R2_001.fastq.gz")],
        gex_reference: resources_test.resolve("reference_gencodev41_chr1/reference_cellranger.tar.gz"),
        vdj_reference: resources_test.resolve("10x_5k_anticmv/raw/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz"),
        feature_reference: resources_test.resolve("10x_5k_anticmv/raw/feature_reference.csv"),
        library_id: ["5k_human_antiCMV_T_TBNK_connect_GEX_1_subset", "5k_human_antiCMV_T_TBNK_connect_AB_subset", "5k_human_antiCMV_T_TBNK_connect_VDJ_subset"],
        library_type: ["Gene Expression", "Antibody Capture", "VDJ"]
      ]
    ])
    | map{ state -> [state.id, state] }
    | cellranger_multi
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, out]"
      assert output[1] instanceof Map : "Output should be a Map."
      // todo: check whether output dir contains fastq files
      "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
}