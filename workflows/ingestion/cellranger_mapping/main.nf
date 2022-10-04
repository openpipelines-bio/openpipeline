nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"
include { cellranger_count } from targetDir + "/mapping/cellranger_count/main.nf"
include { cellranger_count_split } from targetDir + "/mapping/cellranger_count_split/main.nf"
include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataFlowHelper.nf"

config = readConfig("$workflowDir/ingestion/cellranger_mapping/config.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }
}

workflow run_wf {
  take:
  input_ch

  main:

  output_ch = input_ch
  
    // split params for downstream components
    | setWorkflowArguments(
      cellranger_count: [
        "input": "input",
        "expect_cells": "expect_cells",
        "chemistry": "chemistry",
        "secondary_analysis": "secondary_analysis",
        "generate_bam": "generate_bam",
        "include_introns": "include_introns"
      ],
      use_raw_or_filtered: [
        "which_10xh5": "which_10xh5"
      ],
      from_10xh5_to_h5mu: [ 
        "output": "output_h5mu",
        "obsm_metrics": "obsm_metrics",
        "min_genes": "min_genes",
        "min_counts": "min_counts"
      ]
    )

    | getWorkflowArguments(key: "cellranger_count")
    // | view { "cellranger_count: $it" }
    | cellranger_count.run(auto: [ publish: true ])

    // split output dir into map
    // | view { "cellranger_count_split: $it" }
    | cellranger_count_split

    // convert to h5mu
    | pmap { id, data, split_args -> 

      // let user toggle between filtered_h5 or raw_h5
      input = data[split_args.use_raw_or_filtered.which_10xh5]
      
      // combine new data for from_10xh5_to_h5mu
      new_data = 
        [ 
          sample_id: id, 
          input: input, 
          input_metrics_summary: data.metrics_summary
        ] +
        split_args.from_10xh5_to_h5mu

      // store output to third field to return as output
      // and remove the split_args in the meantime
      [ id, new_data, data ]
    }
    // | view { "from_10xh5_to_h5mu: $it" }
    | from_10xh5_to_h5mu.run(auto: [ publish: true ])

    // return output map
    | pmap { id, h5mu, data ->
      [ id, data + [h5mu: h5mu] ]
    }

  emit:
  output_ch
}

/*
# data at various stages of the 'test_wf':
Input: [
  foo, 
  [id:foo, input:[resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq], reference:resources_test/cellranger_tiny_fastq/cellranger_tiny_ref, output_raw:$id.$key.output_raw, output_h5mu:$id.$key.output_h5mu.h5mu, chemistry:auto, secondary_analysis:false, generate_bam:true, include_introns:true, which_10xh5:raw_h5, id_to_obs_names:false, obs_sample_id:sample_id, obsm_metrics:metrics_summary]]cellranger_count: [foo, [id:foo, reference:resources_test/cellranger_tiny_fastq/cellranger_tiny_ref, output_raw:$id.$key.output_raw, input:[resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq], chemistry:auto, generate_bam:true, include_introns:true], [which:[which_10xh5:raw_h5], from_10xh5_to_h5mu:[output:$id.$key.output_h5mu.h5mu, obs_sample_id:sample_id, obsm_metrics:metrics_summary]]
]
cellranger_count_split: [
  foo, 
  ./foo.cellranger_count.output, 
  [which:[which_10xh5:raw_h5], from_10xh5_to_h5mu:[output:$id.$key.output_h5mu.h5mu, obs_sample_id:sample_id, obsm_metrics:metrics_summary]]
]
from_10xh5_to_h5mu: [
  foo, 
  [sample_id:foo, input:./foo.cellranger_count_split.raw_h5.h5, input_metrics_summary:./foo.cellranger_count_split.metrics_summary.csv, output:$id.$key.output_h5mu.h5mu, obs_sample_id:sample_id, obsm_metrics:metrics_summary], 
  [filtered_h5:./foo.cellranger_count_split.filtered_h5.h5, metrics_summary:./foo.cellranger_count_split.metrics_summary.csv, molecule_info:./foo.cellranger_count_split.molecule_info.h5, bam:./foo.cellranger_count_split.bam.bam, bai:./foo.cellranger_count_split.bai.bai, raw_h5:./foo.cellranger_count_split.raw_h5.h5]
]
Output: [
  foo, 
  [filtered_h5:./foo.cellranger_count_split.filtered_h5.h5, metrics_summary:./foo.cellranger_count_split.metrics_summary.csv, molecule_info:./foo.cellranger_count_split.molecule_info.h5, bam:./foo.cellranger_count_split.bam.bam, bai:./foo.cellranger_count_split.bai.bai, raw_h5:./foo.cellranger_count_split.raw_h5.h5, h5mu:./foo.from_10xh5_to_h5mu.output_h5mu.h5mu]
]
*/

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    input: params.resources_test + "/cellranger_tiny_fastq/cellranger_tiny_fastq",
    reference: params.resources_test + "/cellranger_tiny_fastq/cellranger_tiny_ref"
  ]

  output_ch =
    viashChannel(testParams, config)
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, out]"
      assert output[1] instanceof Map : "Output should be a Map."
      // todo: check whether output dir contains fastq files
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}