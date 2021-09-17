nextflow.preview.dsl=2

targetDir = "${params.rootDir}/target/nextflow"

include  { cellranger_mkfastq }  from  "${targetDir}/bcl_demux/cellranger_mkfastq/main.nf"  params(params)

workflow {
    main:

    if (!params.containsKey("input") || params.input == "") {
        exit 1, "ERROR: Please provide a --input parameter pointing to bcl dir or fastq file"
    }
    if (!params.containsKey("samplesheet") || params.samplesheet == "") {
        exit 1, "ERROR: Please provide a --samplesheet parameter"
    }
    if (!params.containsKey("publishDir") || params.publishDir == "") {
        exit 1, "ERROR: Please provide a --publishDir parameter"
    }

    def bcl_= Channel.fromPath(params.input)

    def samplesheet_ = Channel.fromPath(params.samplesheet)

    joined_ = bcl_ \
        | combine(samplesheet_) \
        | map{ it -> [ "input" , [ "input": it[0], "samplesheet": it[1] ], params ] } \
        | cellranger_mkfastq 

    emit:
    joined_
}
