nextflow.enable.dsl=2

include { bd_rhapsody } from params.rootDir + "/target/nextflow/workflows/ingestion/bd_rhapsody/main.nf"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList(
    [
      [
        id: "foo",
        mode: "wta",
        input: file("$resources_test/bdrhap_5kjrt/raw/12WTA*.fastq.gz"),
        reference: resources_test.resolve("reference_gencodev41_chr1/reference_bd_rhapsody.tar.gz"),
        transcriptome_annotation: resources_test.resolve("reference_gencodev41_chr1/reference.gtf.gz"),
        putative_cell_call: "mRNA",
        exact_cell_count: 4900
      ] 
    ])
    | map{ state -> [state.id, state] }
    | bd_rhapsody 
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"

      def id = output[0]
      def data = output[1]

      assert id == "foo" : "Output ID should be same as input ID"
      assert "output_raw" in data : "Output should contain output_raw"
      assert "output_h5mu" in data : "Output should contain output_h5mu"
      assert data.output_h5mu.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toList()
    | view { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
    }
    // | check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}