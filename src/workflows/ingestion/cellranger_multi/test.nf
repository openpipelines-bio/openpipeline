nextflow.enable.dsl=2

include { cellranger_multi } from params.rootDir + "/target/nextflow/workflows/ingestion/cellranger_multi/main.nf"
include { cellranger_multi_test } from params.rootDir + "/target/_test/nextflow/test_workflows/ingestion/cellranger_multi_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "foo",
        input:[
          resources_test.resolve("10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R1_001.fastq.gz"),
          resources_test.resolve("10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R2_001.fastq.gz"),
          resources_test.resolve("10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R1_001.fastq.gz"),
          resources_test.resolve("10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R2_001.fastq.gz"),
          resources_test.resolve("10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R1_001.fastq.gz"),
          resources_test.resolve("10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R2_001.fastq.gz")
        ],
        gex_reference: resources_test.resolve("reference_gencodev41_chr1/reference_cellranger.tar.gz"),
        vdj_reference: resources_test.resolve("10x_5k_anticmv/raw/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz"),
        feature_reference: resources_test.resolve("10x_5k_anticmv/raw/feature_reference.csv"),
        library_id: [
          "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset",
          "5k_human_antiCMV_T_TBNK_connect_AB_subset",
          "5k_human_antiCMV_T_TBNK_connect_VDJ_subset"
        ],
        library_type: [
          "Gene Expression",
          "Antibody Capture",
          "VDJ"
        ]
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

    | cellranger_multi_test.run(
      fromState: ["input": "output_h5mu"]
    )

    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
}

workflow test_wf2 {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "foo",
        input: resources_test.resolve("10x_5k_fixed/raw/"),
        gex_reference: resources_test.resolve("reference_gencodev41_chr1/reference_cellranger.tar.gz"),
        feature_reference: resources_test.resolve("10x_5k_anticmv/raw/feature_reference.csv"),
        library_id: ["4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_subset"],
        library_type: ["Gene Expression"],
        library_lanes: "any",
        probe_barcode_ids: ["BC001|BC002", "BC003", "BC004"],
        gex_generate_bam: false,
        sample_force_cells: [5000, -1, -1],
        probe_set: resources_test.resolve("10x_5k_fixed/raw/Chromium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A_corrected.csv"),
        sample_ids: ["Liver_BC1andOvarian_BC2", "Colorectal_BC3", "Pancreas_BC4"]
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
      assert output_list.size() == 3 : "output channel should contain three events"
      assert (output_list.collect{it[0]} as Set) == (["Liver_BC1andOvarian_BC2", "Pancreas_BC4", "Colorectal_BC3"] as Set) : "Output ID should be same as input ID"
    }
}
