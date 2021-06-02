nextflow.enable.dsl=2

rootDir = "$projectDir/../.."
targetDir = "$rootDir/target/nextflow"

include  { bd_rhapsody_wta }         from  targetDir + "/mapping/bd_rhapsody_wta/main.nf"         params(params)
include  { bd_rhapsody_extracth5ad } from  targetDir + "/mapping/bd_rhapsody_extracth5ad/main.nf" params(params)
include  { publish }                 from  targetDir + "/transfer/publish/main.nf"                params(params)

workflow bd_rhapsody_wta_wf {

    main:
    
    if (!params.containsKey("id") || params.id == "") {
        exit 1, "ERROR: Please provide an --id parameter pointing to your fastq files. For example: --id 'mysample'."
    }
    if (!params.containsKey("input") || params.input == "") {
        exit 1, "ERROR: Please provide a --input1 parameter pointing to your fastq files. For example: --input 'data/*.fastq.gz'."
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
    
    def inputList = params.input.split(";").collect { file(it) }
    
    output_ = Channel.from("") \
        | map { foo -> 
            [   params.id,
                [ "input" : inputList,
                  "reference_genome": file(params.reference_genome),
                  "transcriptome_annotation": file(params.transcriptome_annotation)
                ],
                params
            ]
        } \
        | bd_rhapsody_wta \
        | bd_rhapsody_extracth5ad
        | publish

    emit:
    output_
}
