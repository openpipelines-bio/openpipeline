nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include  { convert_10x_h5_to_h5ad }  from  targetDir + "/convert/convert_10x_h5_to_h5ad/main.nf"  params(params)
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
    if (!params.containsKey("layer") || params.layer == "") {
        print("Setting the layer to default: rna")
        params.layer = "rna"
    }
    if (!params.containsKey("output") || params.output == "") {
        exit 1, "ERROR: Please provide a --output parameter."
    }
    
    print("Converting" + params.input)
    Channel.fromPath(params.input)
        | map { input -> [ input.name, input, params ]}
        | view {it[0]}
        | convert_10x_h5_to_h5ad
        | map { overrideOptionValue(it, "publish", "output", "${params.output}/${it[0]}.h5ad") }
        | publish
}
