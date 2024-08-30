workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map{ id, state ->
      // remove all fields starting with 'output_'
      def newState = state.findAll{k, v -> !k.startsWith("output_")}
      // make sure target is present
      newState.target = newState.target ?: []
      [id, newState]
    }
    | make_reference_component.run(
      fromState: [
        "input": "input",
        "genome_fasta": "genome_fasta", 
        "transcriptome_gtf": "transcriptome_gtf",
        "ercc": "ercc",
        "output_fasta": "output_fasta",
        "output_gtf": "output_gtf",
        "subset_regex": "subset_regex"
      ],
      toState: [
        "output_fasta": "output_fasta",
        "output_gtf": "output_gtf"
      ]
    )
    | build_cellranger_arc_reference.run(
      runIf: { id, state ->
        state.target.contains("cellranger_arc")
      },
      fromState: [
        "genome_fasta": "output_fasta",
        "annotation_gtf": "output_gtf",
        "output": "output_cellranger_arc",
        "motifs_file": "motifs_file",
        "non_nuclear_contigs": "non_nuclear_contigs",
      ],
      toState: [
        "output_cellranger_arc": "output"
      ],
    )
    | build_cellranger_reference.run(
      runIf: { id, state ->
        state.target.contains("cellranger")
      },
      fromState: [
        genome_fasta: "output_fasta",
        transcriptome_gtf: "output_gtf"
      ],
      toState: [
        output_cellranger: "output"
      ]
    )
    | build_star_reference.run(
      runIf: { id, state ->
        state.target.contains("star")
      },
      fromState: [
        genome_fasta: "output_fasta",
        transcriptome_gtf: "output_gtf",
        genomeSAindexNbases: "star_genome_sa_index_nbases"
      ],
      toState: [
        output_star: "output"
      ] 
    )
    | build_bdrhap_reference.run(
      runIf: { id, state ->
        state.target.contains("bd_rhapsody")
      },
      fromState: [
        genome_fasta: "output_fasta",
        gtf: "output_gtf",
        mitochondrial_contigs: "bdrhap_mitochondrial_contigs",
        filtering_off: "bdrhap_filtering_off",
        wta_only_index: "bdrhap_wta_only_index",
        rna_only_index: "bdrhap_rna_only_index"
      ],
      toState: [
        output_bd_rhapsody: "reference_archive"
      ]
    )
    | setState([
      "output_fasta",
      "output_gtf",
      "output_cellranger",
      "output_star",
      "output_bd_rhapsody",
      "output_cellranger_arc",
    ])
  emit:
  output_ch
}