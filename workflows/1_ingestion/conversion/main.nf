nextflow.enable.dsl=2

workflowDir = "${params.rootDir}/workflows"
targetDir = "${params.rootDir}/target/nextflow"

include  { convert_10x_h5_to_h5ad }  from  targetDir + "/convert/convert_10x_h5_to_h5ad/main.nf"  params(params)
include  { publish }                 from  targetDir + "/transfer/publish/main.nf"                params(params)
// include  { overrideOptionValue }     from  workflowDir + "/utils/utils.nf"                        params(params)


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

// A functional approach to 'updating' a value for an option
// in the params Map.
def overrideOptionValue(triplet, _key, _option, _value) {
    mapCopy = triplet[2].toConfigObject().toMap() // As mentioned on https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/config/CascadingConfig.groovy

    return [
        triplet[0],
        triplet[1],
        triplet[2].collectEntries{ function, v1 ->
        (function == _key)
            ? [ (function) : v1.collectEntries{ k2, v2 ->
                (k2 == "arguments")
                    ? [ (k2) : v2.collectEntries{ k3, v3 ->
                        (k3 == _option)
                            ? [ (k3) : v3 + [ "value" : _value ] ]
                            : [ (k3) : v3 ]
                    } ]
                    : [ (k2) : v2 ]
            } ]
            : [ (function), v1 ]
        }
    ]
}


