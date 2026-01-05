process cellrangerMultiStub {
  input:
    tuple val(id), val(unused)

  output:
    tuple val(id), path("raw_dir"), path("*.h5mu"), path("samples.csv")

  script:
    """
    echo "This process is not meant to be run without -stub being defined."
    exit 1
    """

  stub:
    """
    mkdir raw_dir
    touch ${id}_sample_1.h5mu
    touch ${id}_sample_2.h5mu
    touch ${id}_sample_3.h5mu
    echo -e "sample_name,file\n${id}_sample_1,${id}_sample_1.h5mu\n${id}_sample_2,${id}_sample_2.h5mu\n${id}_sample_3,${id}_sample_3.h5mu" > samples.csv
    """
}

workflow run_wf {
  take:
  input_ch

  main:
  cellranger_ch = input_ch
    | cellranger_multi_component.run(
      fromState: [
        "input": "input",
        "gex_input": "gex_input",
        "abc_input": "abc_input",
        "cgc_input": "cgc_input",
        "mux_input": "mux_input",
        "vdj_input": "vdj_input",
        "vdj_t_input": "vdj_t_input",
        "vdj_t_gd_input": "vdj_t_gd_input",
        "vdj_b_input": "vdj_b_input",
        "vdj_denovo": "vdj_denovo",
        "agc_input": "agc_input",
        "library_id": "library_id",
        "library_type": "library_type",
        "library_subsample": "library_subsample",
        "library_lanes": "library_lanes",
        "library_chemistry": "library_chemistry",
        "sample_ids": "sample_ids",
        "sample_description": "sample_description",
        "sample_expect_cells": "sample_expect_cells",
        "sample_force_cells": "sample_force_cells",
        "feature_reference": "feature_reference",
        "feature_r1_length": "feature_r1_length",
        "feature_r2_length": "feature_r2_length",
        "gex_reference": "gex_reference",
        "gex_secondary_analysis": "gex_secondary_analysis",
        "gex_generate_bam": "gex_generate_bam",
        "gex_expect_cells": "gex_expect_cells",
        "gex_force_cells": "gex_force_cells",
        "gex_include_introns": "gex_include_introns",
        "gex_r1_length": "gex_r1_length",
        "gex_r2_length": "gex_r2_length",
        "gex_chemistry": "gex_chemistry",
        "vdj_reference": "vdj_reference",
        "vdj_inner_enrichment_primers": "vdj_inner_enrichment_primers",
        "vdj_r1_length": "vdj_r1_length",
        "vdj_r2_length": "vdj_r2_length",
        "cell_multiplex_oligo_ids": "cell_multiplex_oligo_ids",
        "min_assignment_confidence": "min_assignment_confidence",
        "cmo_set": "cmo_set",
        "barcode_sample_assignment": "barcode_sample_assignment",
        "ocm_barcode_ids": "ocm_barcode_ids",
        "min_crispr_umi": "min_crispr_umi",
        "emptydrops_minimum_umis": "emptydrops_minimum_umis",
        "hashtag_ids": "hashtag_ids",
        "probe_set": "probe_set",
        "filter_probes": "filter_probes",
        "probe_barcode_ids": "probe_barcode_ids",
        "control_id": "control_id",
        "mhc_allele": "mhc_allele",
        "check_library_compatibility": "check_library_compatibility",
        "output": "output_raw",
      ],
      toState: [
        "output_raw": "output", 
        "input": "output"
      ]
    )

  h5mu_ch = cellranger_ch
    | from_cellranger_multi_to_h5mu.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "uns_metrics": state.uns_metrics,
          "output_compression": "gzip"
        ]
      },
      toState: {id, output, state -> 
        [
          "sample_csv": output.sample_csv,
          "output_h5mu": output.output,
          "output_raw": state.output_raw
        ]
      },
      filter: {!workflow.stubRun},
    )
 
  h5mu_stub_ch = cellranger_ch
    | filter{workflow.stubRun}
    | map {id, state -> [id, state.input]}
    | cellrangerMultiStub
    | map {id, raw_dir, output, output_files ->
      [id, ["output_raw": raw_dir, "output_h5mu": output, "sample_csv": output_files]]
    }

  output_ch = h5mu_stub_ch.concat(h5mu_ch)
    | flatMap {id, state ->
      def h5mu_list = state.output_h5mu
      def samples = readCsv(state.sample_csv.toUriString())
      println "Samples: $samples" 
      def result = h5mu_list.collect{ h5mu_file ->
        println "H5mu: ${h5mu_file}, getName: ${h5mu_file.getName()}"
        def corresponding_csv_entry = samples.find{h5mu_file.getName() == it.file}
        print "CSV entry: $corresponding_csv_entry"
        // The cellranger component used to only be able to output a single h5mu file
        // In cases where cell multiplexing is not used (1 output sample), it uses 'run' for the sample ID as a dummy.
        // This sample ID 'run' was never used for the ID of the channel events.
        // So here we overwrite this 'run' id with the name of the input event.
        def new_id = h5mu_list.size() == 1 ? id : corresponding_csv_entry.sample_name
        return [ new_id, ["output_h5mu": h5mu_file, "output_raw": state.output_raw, "_meta": ["join_id": id]]]
      }
      return result
    }

  emit:
  output_ch
}
