nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { make_reference } from targetDir + "/reference/make_reference/main.nf"
include { build_bdrhap_reference } from targetDir + "/reference/build_bdrhap_reference/main.nf"
include { star_build_reference } from targetDir + "/mapping/star_build_reference/main.nf"
include { build_cellranger_reference } from targetDir + "/reference/build_cellranger_reference/main.nf"

include { readConfig; channelFromParams; preprocessInputs; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/ingestion/make_reference/config.vsh.yaml")

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  
  ref_ch = input_ch
    // split params for downstream components
    | preprocessInputs("config": config)
    | setWorkflowArguments(
      make_reference: [
        "genome_fasta": "genome_fasta", 
        "transcriptome_gtf": "transcriptome_gtf",
        "ercc": "ercc",
        "output_fasta": "output_fasta",
        "output_gtf": "output_gtf",
        "subset_regex": "subset_regex"
      ],
      cellranger: [
        "output": "output_cellranger",
        "target": "target"
      ],
      bd_rhapsody: [
        "output": "output_bd_rhapsody",
        "target": "target"
      ],
      star: [
        "output": "output_star",
        "target": "target"
      ]
    )

    // generate reference
    | getWorkflowArguments(key: "make_reference")
    | make_reference.run(auto: [ publish: true ])


  // generate cellranger index (if so desired)
  cellranger_ch = ref_ch
    | getWorkflowArguments(key: "cellranger")
    | filter{ "cellranger" in it[1].target }
    | build_cellranger_reference.run(
      renameKeys: [ genome_fasta: "output_fasta", transcriptome_gtf: "output_gtf" ], 
      auto: [ publish: true ]
    )
    | map{ tup -> tup.take(2) }

  // generate bd_rhapsody index (if so desired)
  bd_rhapsody = ref_ch
    | getWorkflowArguments(key: "bd_rhapsody")
    | filter{ "bd_rhapsody" in it[1].target }
    | build_bdrhap_reference.run(
      renameKeys: [ genome_fasta: "output_fasta", transcriptome_gtf: "output_gtf" ], 
      auto: [ publish: true ]
    )
    | map{ tup -> tup.take(2) }

  // generate star index (if so desired)
  star = ref_ch
    | getWorkflowArguments(key: "star")
    | filter{ "star" in it[1].target }
    | star_build_reference.run(
      renameKeys: [ genome_fasta: "output_fasta", transcriptome_gtf: "output_gtf" ], 
      auto: [ publish: true ]
    )
    | map{ tup -> tup.take(2) }
  
  // merge everything together
  passthr_ch = input_ch
    | map{ tup -> [ tup[0] ] + tup.drop(2) }

  output_ch = ref_ch
    | map{ tup -> tup.take(2) }
    | join(cellranger_ch, remainder: true)
    | join(bd_rhapsody, remainder: true)
    | join(star, remainder: true)
    | join(passthr_ch)
    | map{ tup -> 
      id = tup[0]
      data = tup[1] + [ output_cellranger: tup[2], output_bd_rhapsody: tup[3], output_star: tup[4] ]
      data = data.findAll{it.value != null} // remove empty fields
      psthr = tup.drop(4)
      [ id, data ] + psthr
    }

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  params.param_list = [
    [
      id: "gencode_v41_ercc",
      genome_fasta: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz",
      transcriptome_gtf: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz",
      ercc: "https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip",
      subset_regex: "(ERCC-00002|chr20)",
      target: ["cellranger", "bd_rhapsody", "star"]
    ]
  ]

  output_ch =
    channelFromParams(params, config)
    | view{ "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 3 : "outputs should contain two elements; [id, file, passthrough]"
      assert output[1].size() == 5 : "output data should contain 5 elements"
      // todo: check output data tuple
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "There should be one output"
    }
}
