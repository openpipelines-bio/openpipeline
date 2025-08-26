nextflow.enable.dsl=2

include { pseudobulk_deseq2 } from params.rootDir + "/target/nextflow/workflows/differential_expression/pseudobulk_deseq2/main.nf"
include { pseudobulk_deseq2_test } from params.rootDir + "/target/_test/nextflow/test_workflows/differential_expression/pseudobulk_deseq2_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file(params.resources_test)

  output_ch = Channel.fromList(
    [
      [
        id: "simple_execution_test",
        input: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        obs_cell_group: "cell_type",
        obs_groups: ["treatment", "disease"],
        min_obs_per_sample: 5,
        design_formula: "~ treatment",
        contrast_column: "treatment",
        contrast_values: ["ctrl", "stim"],
        output: "simple_execution_test_output"
      ]
    ])
    | map{ state -> [state.id, state] }
    | pseudobulk_deseq2
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id.endsWith("_test") : "Output ID should be same as input ID"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'. found: ${state}"
      assert state.output.isDirectory() : "'output' should be a directory."
      def deseqFiles = state.output.listFiles().findAll { file ->
        file.name.matches(/deseq2_analysis_.*\.csv/)
      }
      assert deseqFiles.size() == 4, "Output directory should contain exactly 4 deseq2_analysis_*.csv files. Found: ${deseqFiles.size()} files: ${deseqFiles.collect { it.name }}"
    
    "Output: $output"
    }

    | flatMap { output->
        def id = output[0]
        def state = output[1]
        
        // Find all CSV files and create separate entries
        def deseqFiles = state.output.listFiles().findAll { file ->
          file.name.matches(/deseq2_analysis_.*\.csv/)
        }
        
        deseqFiles.collect { csvFile ->
          def cellType = csvFile.name.replaceAll(/deseq2_analysis_(.*)\.csv/, '$1')
          def newId = "${id}_${cellType}"
          [newId, ["input": csvFile]]
        }
    }

    | view { output -> "After flatMap: $output" }
    | pseudobulk_deseq2_test.run(
        fromState: [
          "input": "input"
        ]
    )

  output_ch
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 4 : "output channel should contain 4 events"

      def sortedIds = output_list.collect{it[0]}.sort()
      def expectedIds = [
        "simple_execution_test_classical_monocyte",
        "simple_execution_test_erythrocyte",
        "simple_execution_test_neutrophil",
        "simple_execution_test_plasma_cell"
      ]

      assert sortedIds == expectedIds : "IDs don't match. Expected: ${expectedIds}, found: ${sortedIds}"
    }
}