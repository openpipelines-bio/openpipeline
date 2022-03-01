nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf" params(params)
include { filter_with_scrublet } from targetDir + "/filter/filter_with_scrublet/main.nf" params(params)
include { lognorm } from targetDir + '/normalize/lognorm/main.nf' params(params)
include { hvg_scanpy } from targetDir + '/hvg/hvg_scanpy/main.nf' params(params)
include { pca } from targetDir + '/dimred/pca/main.nf' params(params)
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf' params(params)
include { umap } from targetDir + '/dimred/umap/main.nf' params(params)
include { leiden } from targetDir + '/cluster/leiden/main.nf' params(params)

include { publish } from targetDir + "/transfer/publish/main.nf" params(params)
include { overrideOptionValue } from workflowDir + "/utils/utils.nf" params(params)


workflow {
    main:
    
    if (!params.containsKey("input") || params.input == "") {
        exit 1, "ERROR: Please provide a --input parameter pointing to the count matrices to be converted"
    }
    if (!params.containsKey("output") || params.output == "") {
        exit 1, "ERROR: Please provide a --output parameter."
    }
            
    Channel.fromPath(params.input)
        | map { input -> [ input.name, input, params ]}
        | map { overrideOptionValue(it, "filter_with_counts", "modality", "rna") }
        | view { "before filter: ${it[0]} - ${it[1]}" }
        | filter_with_counts
        | view { "after filter: ${it[0]} - ${it[1]}" }
        // | scrublet
        // | lognorm           
        // | hvg_scanpy        
        // | pca               
        // | find_neighbors    
        // | leiden            
        // | umap              
        | map { overrideOptionValue(it, "publish", "output", "${params.output}/${it[0]}.h5mu") } 
        | publish
 
}
