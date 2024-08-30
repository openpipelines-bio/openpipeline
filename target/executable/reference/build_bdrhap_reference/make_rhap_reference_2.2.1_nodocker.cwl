requirements:
  InlineJavascriptRequirement: {}
class: CommandLineTool
label: Reference Files Generator for BD Rhapsodyâ„¢ Sequencing Analysis Pipeline
cwlVersion: v1.2
doc: >- 
    The Reference Files Generator creates an archive containing Genome Index and Transcriptome annotation files needed for the BD Rhapsodyâ„¢ Sequencing Analysis Pipeline. The app takes as input one or more FASTA and GTF files and produces a compressed archive in the form of a tar.gz file. The archive contains:\n  - STAR index\n  - Filtered GTF file


baseCommand: run_reference_generator.sh 
inputs: 
    Genome_fasta:
        type: File[]
        label: Reference Genome
        doc: |-
            Reference genome file in FASTA format. The BD Rhapsodyâ„¢ Sequencing Analysis Pipeline uses GRCh38 for Human and GRCm39 for Mouse.
        inputBinding:
            prefix: --reference-genome
            shellQuote: false
    Gtf:
        type: File[]
        label: Transcript Annotations
        doc: |-
            Transcript annotation files in GTF format. The BD Rhapsodyâ„¢ Sequencing Analysis Pipeline uses Gencode v42 for Human and M31 for Mouse.
        inputBinding:
            prefix: --gtf
            shellQuote: false
    Extra_sequences:
        type: File[]?
        label: Extra Sequences
        doc: |-
            Additional sequences in FASTA format to use when building the STAR index. (E.g. phiX genome)
        inputBinding:
            prefix: --extra-sequences
            shellQuote: false
    Mitochondrial_Contigs:
        type: string[]?
        default: ["chrM", "chrMT", "M", "MT"]
        label: Mitochondrial Contig Names
        doc: |-
            Names of the Mitochondrial contigs in the provided Reference Genome. Fragments originating from contigs other than these are identified as 'nuclear fragments' in the ATACseq analysis pipeline.
        inputBinding:
            prefix: --mitochondrial-contigs
            shellQuote: false
    Filtering_off:
        type: boolean?
        label: Turn off filtering
        doc: |-
            By default the input Transcript Annotation files are filtered based on the gene_type/gene_biotype attribute. Only features having the following attribute values are are kept:
            - protein_coding
            - lncRNA (lincRNA and antisense for Gencode < v31/M22/Ensembl97)
            - IG_LV_gene
            - IG_V_gene
            - IG_V_pseudogene
            - IG_D_gene
            - IG_J_gene
            - IG_J_pseudogene
            - IG_C_gene
            - IG_C_pseudogene
            - TR_V_gene
            - TR_V_pseudogene
            - TR_D_gene
            - TR_J_gene
            - TR_J_pseudogene
            - TR_C_gene
            If you have already pre-filtered the input Annotation files and/or wish to turn-off the filtering, please set this option to True.
        inputBinding: 
            prefix: --filtering-off
            shellQuote: false
    WTA_Only:
        type: boolean?
        label: WTA only index
        doc: Build a WTA only index, otherwise builds a WTA + ATAC index.
        inputBinding:
            prefix: --wta-only-index
            shellQuote: false
    Archive_prefix:
        type: string?
        label: Archive Prefix
        doc: |-
            A prefix for naming the compressed archive file containing the Reference genome index and annotation files. The default value is constructed based on the input Reference files.
        inputBinding:
            prefix: --archive-prefix
            shellQuote: false
    Extra_STAR_params:
        type: string?
        label: Extra STAR Params
        doc: |-
            Additional parameters to pass to STAR when building the genome index. Specify exactly like how you would on the command line.
            Example:
              --limitGenomeGenerateRAM 48000 --genomeSAindexNbases 11
        inputBinding:
            prefix: --extra-star-params 
            shellQuote: true
  
    Maximum_threads:
        type: int?
        label: Maximum Number of Threads
        doc: |-
            The maximum number of threads to use in the pipeline. By default, all available cores are used.
        inputBinding:
            prefix: --maximum-threads
            shellQuote: false

outputs:

    Archive:
        type: File
        doc: |- 
            A Compressed archive containing the Reference Genome Index and annotation GTF files. This archive is meant to be used as an input in the BD Rhapsodyâ„¢ Sequencing Analysis Pipeline.
        id: Reference_Archive
        label: Reference Files Archive
        outputBinding:
            glob: '*.tar.gz'

