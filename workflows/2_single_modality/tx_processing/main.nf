nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include  { filter_on_umi_genes_and_mito }  from  targetDir + "/filter/filter_on_umi_genes_and_mito/main.nf"    params(params)
include  { find_doublets_using_scrublet }  from  targetDir + "/doublets/find_doublets_using_scrublet/main.nf"  params(params)
include  { lognorm }       from  targetDir + '/normalize/lognorm/main.nf'             params(params)
include  { hvg_scanpy }    from  targetDir + '/hvg/hvg_scanpy/main.nf'                params(params)
include  { pca }           from  targetDir + '/dimred/pca/main.nf'                    params(params)
include  { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf'      params(params)
include  { umap }          from  targetDir + '/dimred/umap/main.nf'                   params(params)
include  { leiden }        from  targetDir + '/cluster/leiden/main.nf'                params(params)

include  { publish }                 from  targetDir + "/transfer/publish/main.nf"    params(params)
include  { overrideOptionValue }     from  workflowDir + "/utils/utils.nf"            params(params)


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
        | filter_on_umi_genes_and_mito
        | find_doublets_using_scrublet
        | lognorm           
        | hvg_scanpy        
        | pca               
        | find_neighbors    
        | leiden            
        | umap              
        | map { overrideOptionValue(it, "publish", "output", "${params.output}/${it[0]}.h5ad") } 
        | publish
 
}
