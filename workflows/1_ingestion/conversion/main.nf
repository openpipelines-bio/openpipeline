nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf" 
include { from_10xmtx_to_h5mu } from targetDir + "/convert/from_10xmtx_to_h5mu/main.nf" 

include { publish } from targetDir + "/transfer/publish/main.nf"               
include { overrideOptionValue } from workflowDir + "/utils/utils.nf"                       


workflow {
    main:
    
    if (!params.containsKey("input_type") || params.input_type == "") {
        exit 1, "ERROR: Please provide a --input_type parameter for the conversion."
    }
    if (!params.containsKey("input") || params.input == "") {
        exit 1, "ERROR: Please provide a --input parameter pointing to the count matrices to be converted"
    }
    if (!params.containsKey("output") || params.output == "") {
        exit 1, "ERROR: Please provide a --output parameter."
    }

    input = Channel.fromPath(params.input)
        | map { input -> [ input.baseName, input ]}

    switch(params.input_type) { 
        case "10xh5":
            middle = input | from_10xh5_to_h5mu
            break

        case "10xmtx":
            middle = input | from_10xmtx_to_h5mu
            break

        default:
            exit 1, "ERROR: Unrecognised --input_type."
    }

    output = middle
        | map { overrideOptionValue(it, "publish", "output", "${params.output}/${it[0]}.h5mu") }
        | publish.run(
            map: { [it[0], [input: it[1], output: "${it[0]}.h5mu"]] },
            auto: [ publish: true ]
        )

}
