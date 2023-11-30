nextflow.enable.dsl=2

workflowDir = params.rootDir + "/src/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { bd_rhapsody } from targetDir + "/mapping/bd_rhapsody/main.nf"
include { from_bdrhap_to_h5mu } from targetDir + "/convert/from_bdrhap_to_h5mu/main.nf"

include { readConfig; channelFromParams; preprocessInputs; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

config = readConfig("$workflowDir/ingestion/bd_rhapsody/config.vsh.yaml")

workflow bd_rhapsody_entrypoint {
  helpMessage(config)

  channelFromParams(params, config)
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | preprocessInputs("config": config)

    // run bd rhapsody
    | bd_rhapsody.run(
      auto: [ publish: true ],
      fromState: { id, state ->
        // pass all arguments except:
        //  - remove output_h5mu and output_compression
        //  - rename output_raw to output
        def data_ = state.clone()
        data_.remove("output_h5mu")
        data_.remove("output_raw")
        data_.remove("output_compression")
        data_ + [ output: state.output_raw ]
      },
      toState: { id, data, state ->
        state + [ output_raw: data.output ]
      }
    )
    | view {"After bd_rhapsody: $it"}

    // convert to h5mu
    | from_bdrhap_to_h5mu.run(
      fromState: { id, state ->
        [
          id: id,
          input: state.output_raw,
          output: state.output_h5mu,
          output_compression: "gzip"
        ]
      },
      toState: { id, data, state ->
        [
          output_raw: state.output_h5mu,
          output_h5mu: data.output
        ]
      },
      auto: [publish: true]
    )

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    mode: "wta",
    input: params.resources_test + "/bdrhap_5kjrt/raw/12WTA*.fastq.gz",
    reference: params.resources_test + "/reference_gencodev41_chr1/reference_bd_rhapsody.tar.gz",
    transcriptome_annotation: params.resources_test + "/reference_gencodev41_chr1/reference.gtf.gz",
    putative_cell_call: "mRNA",
    exact_cell_count: 4900
  ]

  output_ch =
    channelFromParams(testParams, config)
      | view { "Input: $it" }
      | run_wf
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