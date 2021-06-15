nextflow.enable.dsl=2

workflowDir = "${params.rootDir}/workflows"
targetDir = "${params.rootDir}/target/nextflow"

include  { bd_rhapsody_wta }         from  targetDir + "/mapping/bd_rhapsody_wta/main.nf"         params(params)
include  { bd_rhapsody_extracth5ad } from  targetDir + "/mapping/bd_rhapsody_extracth5ad/main.nf" params(params)
include  { publish }                 from  targetDir + "/transfer/publish/main.nf"                params(params)
include  { overrideOptionValue }     from  workflowDir + "/utils/utils.nf"                        params(params)

/* BD Rhapsody WTA - common workflow
 * 
 * consumed params:
 *   reference_genome:              a path to STAR index as a tar.gz file
 *   transcriptome_annotation:      a path to GTF annotation file
 *   output                         a publish dir for the output h5ad files
 * input format:                [ id, path_string, params ]
 *   value id:                      a sample id for one or more fastq files
 *   value path_string:             one or more fastq paths, separated with semicolons, paths may be globs
 *   value params:                  the params object, which may already have sample specific overrides
 * output format:               [ id, h5ad, params ]
 *   value id:                      same as input
 *   value h5ad:                    h5ad object of mapped fastq reads
 *   value params:                  same as input params, but the publish__output parameter has been overridden
 * publishes:
 *   the output h5ad files
 */
workflow main_wf {
    take: input_

    main:
    if (!params.containsKey("reference_genome") || params.reference_genome == "") {
        exit 1, "ERROR: Please provide a --reference_genome parameter pointing to the STAR index for tar.gz format.\nSee BD Genomics WTA Rhapsody analysis pipeline documentation for instructions to obtain pre-built STAR index file."
    }
    if (!params.containsKey("transcriptome_annotation") || params.transcriptome_annotation == "") {
        exit 1, "ERROR: Please provide a --transcriptome_annotation parameter pointing to the GTF annotation file."
    }
    if (!params.containsKey("output") || params.output == "") {
        exit 1, "ERROR: Please provide a --output parameter."
    }
    
    output_ = input_ \
        | map { it -> 
            [   it[0],
                [ "input" : it[1].split(";").collect { path -> file(path) }.flatten(),
                  "reference_genome": file(params.reference_genome),
                  "transcriptome_annotation": file(params.transcriptome_annotation)
                ],
                it[2]
            ]
        } \
        | map { overrideOptionValue(it, "publish", "output", "${params.output}/${it[0]}/${it[0]}.h5ad") } \
        | view{ it[1] } \
        | bd_rhapsody_wta \
        | bd_rhapsody_extracth5ad \
        | publish

    emit: output_
}

/* BD Rhapsody WTA - single-sample workflow
 * 
 * consumed params:
 *   id:                            a sample id for one or more fastq files
 *   input:                         one or more fastq paths, separated with semicolons, paths may be globs
 *   reference_genome:              a path to STAR index as a tar.gz file
 *   transcriptome_annotation:      a path to GTF annotation file
 *   output                         a publish dir for the output h5ad files
 * output format:               [ id, h5ad, params ]
 *   value id:                      a sample id for one or more fastq files
 *   value h5ad:                    h5ad object of mapped fastq reads
 *   value params:                  the params object, which may already have sample specific overrides
 * publishes:
 *   the output h5ad files
 */
workflow single_wf {
    main:
    if (!params.containsKey("id") || params.id == "") {
        exit 1, "ERROR: Please provide an --id parameter pointing to your fastq files. For example: --id 'mysample'."
    }
    if (!params.containsKey("input") || params.input == "") {
        exit 1, "ERROR: Please provide a --input parameter pointing to your fastq files. For example: --input 'data/*.fastq.gz'."
    }
    
    output_ = Channel.from(params.id) \
        | map { [ it, params.input, params ] } \
        | main_wf

    emit: output_
}

/* BD Rhapsody WTA - multi-sample workflow
 * 
 * consumed params:
 *   tsv:                           a tsv for processing multiple input files. tsv must contain two columns, 'id' and 'input'. 
 *                                  'id' is a sample id for one or more fastq files. 
 *                                  'input' is one or more fastq paths, separated with semicolons, paths may be globs
 *   reference_genome:              a path to STAR index as a tar.gz file
 *   transcriptome_annotation:      a path to GTF annotation file
 *   output                         a publish dir for the output h5ad files
 * output format:               [ id, h5ad, params ]
 *   value id:                      a sample id for one or more fastq files
 *   value h5ad:                    h5ad object of mapped fastq reads
 *   value params:                  the params object, which may already have sample specific overrides
 * publishes:
 *   the output h5ad files
 */
workflow multi_wf {
    main:
    if (!params.containsKey("tsv") || params.tsv == "") {
        exit 1, "ERROR: Please provide an --tsv parameter. The tsv must include two columns 'id' and 'input'."
    }
    
    print("Starting multisample workflow")
    
    output_ = Channel.fromPath(file(params.tsv)) \
        | splitCsv(header: true, sep: "\t") \
        | map { tsv -> [ tsv.id, tsv.input, params ] } \
        | main_wf

    emit: output_
}
