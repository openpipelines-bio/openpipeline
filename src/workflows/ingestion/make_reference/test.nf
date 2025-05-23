nextflow.enable.dsl=2

include { make_reference } from params.rootDir + "/target/_test/nextflow/workflows/ingestion/make_reference/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "gencode_v41_ercc",
        genome_fasta: resources_test.resolve("reference_gencodev41_chr1/reference.fa.gz"),
        transcriptome_gtf: resources_test.resolve("reference_gencodev41_chr1/reference.gtf.gz"),
        ercc: resources_test.resolve("reference_gencodev41_chr1/ERCC92.zip"),
        subset_regex: "(ERCC-00002|chr1)",
        target: ["cellranger", "bd_rhapsody", "star"]
      ]        
    ])
    | map{ state -> [state.id, state] }
    | make_reference
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].size() == 5 : "output data should contain 5 elements"
      // todo: check output data tuple
      "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "There should be one output"
    }
}
