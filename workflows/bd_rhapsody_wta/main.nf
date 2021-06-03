nextflow.enable.dsl=2

rootDir = "$projectDir/../.."
workflowDir = "$rootDir/workflows"
targetDir = "$rootDir/target/nextflow"

include  { bd_rhapsody_wta }         from  targetDir + "/mapping/bd_rhapsody_wta/main.nf"         params(params)
include  { bd_rhapsody_extracth5ad } from  targetDir + "/mapping/bd_rhapsody_extracth5ad/main.nf" params(params)
include  { publish }                 from  targetDir + "/transfer/publish/main.nf"                params(params)
include  { overrideOptionValue }     from  workflowDir + "/utils/utils.nf"                        params(params)

def splitFile(str) {
    str.split(";").collect { file(it) }
}

workflow main_wf {
    // input format: [ id, [input, reference_genome, transcriptome_annotation], params ]
    take: input_

    main:
    
    def inputList = params.input.split(";").collect { file(it) }
    
    output_ = input_ \
        | bd_rhapsody_wta \
        | bd_rhapsody_extracth5ad

    // ouput format: [ id, output, params ]
    emit: output_
}

workflow cli_wf {
    main:
    if (!params.containsKey("id") || params.id == "") {
        exit 1, "ERROR: Please provide an --id parameter pointing to your fastq files. For example: --id 'mysample'."
    }
    if (!params.containsKey("input") || params.input == "") {
        exit 1, "ERROR: Please provide a --input parameter pointing to your fastq files. For example: --input 'data/*.fastq.gz'."
    }
    if (!params.containsKey("reference_genome") || params.reference_genome == "") {
        exit 1, "ERROR: Please provide a --reference_genome parameter pointing to the STAR index for tar.gz format.\nSee BD Genomics WTA Rhapsody analysis pipeline documentation for instructions to obtain pre-built STAR index file."
    }
    if (!params.containsKey("transcriptome_annotation") || params.transcriptome_annotation == "") {
        exit 1, "ERROR: Please provide a --transcriptome_annotation parameter pointing to the GTF annotation file."
    }
    if (!params.containsKey("output") || params.output == "") {
        exit 1, "ERROR: Please provide a --output parameter."
    }
    
    Channel.from("") \
        | map { foo -> 
            [   params.id,
                [ "input" : splitFile(params.input),
                  "reference_genome": file(params.reference_genome),
                  "transcriptome_annotation": file(params.transcriptome_annotation)
                ],
                params
            ]
        } \
        | map { overrideOptionValue(it, "publish", "output", "${params.output}/${it[0]}/${it[0]}.h5ad") } \
        | main_wf \
        | publish
}

workflow tsv_wf {
    main:
    if (!params.containsKey("tsv") || params.tsv == "") {
        exit 1, "ERROR: Please provide an --tsv parameter. The tsv must include two columns 'id' and 'input'."
    }
    if (!params.containsKey("reference_genome") || params.reference_genome == "") {
        exit 1, "ERROR: Please provide a --reference_genome parameter pointing to the STAR index for tar.gz format.\nSee BD Genomics WTA Rhapsody analysis pipeline documentation for instructions to obtain pre-built STAR index file."
    }
    if (!params.containsKey("transcriptome_annotation") || params.transcriptome_annotation == "") {
        exit 1, "ERROR: Please provide a --transcriptome_annotation parameter pointing to the GTF annotation file."
    }
    if (!params.containsKey("output") || params.output == "") {
        exit 1, "ERROR: Please provide a --output parameter."
    }
    
    Channel.fromPath(file(params.tsv)) \
        | splitCsv(header: true, sep: "\t") \
        | map { row -> 
            [   row.id,
                [ "input" : splitFile(row.input),
                  "reference_genome": file(params.reference_genome),
                  "transcriptome_annotation": file(params.transcriptome_annotation)
                ],
                params
            ]
        } \
        | map { overrideOptionValue(it, "publish", "output", "${params.output}/${it[0]}/${it[0]}.h5ad") } \
        | main_wf \
        | publish
}
