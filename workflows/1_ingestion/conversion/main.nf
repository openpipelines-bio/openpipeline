nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include  { convert_10x_to_h5mu }  from  targetDir + "/convert/convert_10x_to_h5mu/main.nf"  params(params)
include  { convert_10x_mtx_to_h5mu }  from  targetDir + "/convert/convert_10x_mtx_to_h5mu/main.nf"  params(params)

include  { publish }                 from  targetDir + "/transfer/publish/main.nf"                params(params)
include  { overrideOptionValue }     from  workflowDir + "/utils/utils.nf"                        params(params)


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

    switch(params.input_type) { 
        case "10x_h5":
            Channel.fromPath(params.input)
                | map { input -> [ input.name, input, params ]}
                | convert_10x_to_h5mu
                | map { overrideOptionValue(it, "publish", "output", "${params.output}/${it[0]}.h5mu") }
                | publish
            break

        case "10x_mtx":
            Channel.fromPath(params.input)
                | map { input -> [ input.name, input, params ]}
                | convert_10x_mtx_to_h5mu
                | map { overrideOptionValue(it, "publish", "output", "${params.output}/${it[0]}.h5mu") }
                | publish
            break


        default:
            exit 1, "WARNING: There was no input_type recognize. Please use the --input_type parameter to set the input's input format."
    }

}
