requirements:
  InlineJavascriptRequirement: {}
hints:
  - class: 'DockerRequirement'
    dockerPull: bdgenomics/rhapsody:4.0
class: CommandLineTool
label: Reference Files Generator for BD Rhapsody™ Sequencing Analysis Pipeline
cwlVersion: v1.2
doc: >-
    The Reference Files Generator creates an archive containing Genome Index and Transcriptome annotation files needed for the BD Rhapsody™ Sequencing Analysis Pipeline. The app takes as input one or more FASTA and GTF files and produces a compressed archive in the form of a tar.gz file. The archive contains:\n  - STAR index\n  - Filtered GTF file


baseCommand: run_reference_generator.sh
inputs:
    Genome_fasta:
        type: File[]
        label: Reference Genome
        doc: |-
            Reference genome file in FASTA format. The pre-built BD reference archive files use GRCh38 for Human and GRCm39 for Mouse.
        inputBinding:
            prefix: --reference-genome
            shellQuote: false
    Gtf:
        type: File[]
        label: Transcript Annotations
        doc: |-
            Transcript annotation files in GTF format. The pre-built BD reference archive files use Gencode v49 for Human and M38 for Mouse.
        inputBinding:
            prefix: --gtf
            shellQuote: false
    Extra_sequences:
        type: File[]?
        label: Extra Sequences
        doc: |-
            Additional sequences in FASTA format to use when building the STAR and bwa-mem2 indices. (E.g. phiX genome)
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
    Transcription_Factor_Motif_PFM:
        type: File[]?
        label: Transcription Factor Motif PFM
        doc: |-
            Text file of Transcription Factor Motif position frequency matrices in JASPAR format.
            The pre-built BD reference archive files use the JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt file for both Human and Mouse.
            You can browse the list of all files here : https://jaspar.elixir.no/download/data/2024/CORE/
        inputBinding:
            prefix: --transcription-factor-motif-pfm
            shellQuote: false
    Disable_Biotype_Filtering:
        type: boolean?
        label: Turn off biotype filtering
        doc: |-
            By default, the input Gtf files are filtered based on the gene_type, gene_biotype, or transcript_type attributes.
            Only genes with a gene_type or gene_biotype matching one of the accepted values are retained, and each gene must have at least one transcript with a valid transcript_type.
            An exception is made to the transcript requirement if the gene does not overlap with any other gene that has a valid biotype.

            Accepted biotypes include:

            - protein_coding
            - protein_coding_LOF
            - lncRNA (lincRNA and antisense for Gencode < v31/M22/Ensembl97)
            - IG_LV_gene
            - IG_V_gene
            - IG_V_pseudogene
            - IG_D_gene
            - IG_D_pseudogene
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
            prefix: --disable-biotype-filtering
            shellQuote: false
    Disable_Readthrough_Filtering:
        type: boolean?
        label: Turn off filtering of readthrough transcripts
        doc: |-
            By default genes with only readthrough transcripts are removed. Any readthrough_transcript feature is also removed if its parent gene overlaps with another gene that meets the biotype requirement. Please set this option to True to disable this behaviour.
        inputBinding:
            prefix: --disable-readthrough-filtering
            shellQuote: false
    Filter_PARs:
        type: boolean?
        label: Turn on filtering of features in Pseudo-Autosomal Regions (PARs)
        doc: |-
            This applies to only a Human build 38 reference. If enabled, features in the 2 PARs on the Y chromosome are removed.
        inputBinding:
            prefix: --filter-PARs
            shellQuote: false
    Combine_Same_Gene_Names:
        type: boolean?
        label: Combine genes with the same gene name, even when gene ids are different
        doc: |-
            Combine genes with the same gene name, even when gene ids are different. If not enabled (default), chromosome, start, and strand information will be appended to the duplicated gene names.
        inputBinding:
            prefix: --combine-same-gene-names
            shellQuote: false
    WTA_Only:
        type: boolean?
        label: WTA only index
        doc: Build a WTA only index, otherwise builds a WTA + ATAC index.
        inputBinding:
            prefix: --wta-only-index
            shellQuote: false
    Disable_FASTA_Validation:
        type: boolean?
        label: Disable FASTA validation
        doc: |-
            If True, skip validation that FASTA files contain chromosomes referenced in GTF files. Make this option True only in rare cases where you deliberately want to bypass the chromosome compatibility check. Default is False. 
        inputBinding:
            prefix: --disable-fasta-validation
            shellQuote: false
    Repair_GTF:
        type: boolean?
        label: Repair invalid GTF files
        doc: |-
            Attempt to repair invalid GTF files using the built-in repair functionality. If True, the validator will try to fix common GTF format issues and use the repaired files for downstream processing. Default is False. 
        inputBinding:
            prefix: --repair-gtf
            shellQuote: false
    GTF_Validation_Complete_Mode:
        type: boolean?
        label: Generate a full report about GTF
        doc: |-
            Generates a full report about errors and warnings in GTF, instead of failing when the first fatal error is encountered. Default is False. 
        inputBinding:
            prefix: --gtf-validation-complete-mode
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
              --limitGenomeGenerateRAM 48000000000 --genomeSAindexNbases 11
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
            A Compressed archive containing the Reference Genome Index and annotation GTF files. This archive is meant to be used as an input in the BD Rhapsody™ Sequencing Analysis Pipeline.
        id: Reference_Archive
        label: Reference Files Archive
        outputBinding:
            glob: '*.tar.gz'

