#!/bin/bash

set -eo pipefail

# add additional params
extra_params=( )

if [ ! -z "$par_bam" ]; then
  extra_params+=( "--bam $par_bam" )
fi
  
if [ ! -z "$par_bam_list" ]; then
  extra_params+=( "--bam-list $par_bam_list" )
fi

if [ "$par_stdin" = true ]; then
  extra_params+=( "--stdin" )
fi

if [ ! -z "$par_targets" ] ; then
  extra_params+=( "--targets $par_targets" )
fi

if [ ! -z "$par_region" ]; then
  extra_params+=( "--region $par_region" )
fi

if [ ! -z "$par_max_coverage" ]; then
  extra_params+=( "--max-coverage $par_max_coverage" )
fi

if [ ! -z "$par_samples" ]; then
  extra_params+=( "--samples $par_samples" )
fi

if [ ! -z "$par_populations" ]; then
  extra_params+=( "--poplutions $par_populations" )
fi

if [ ! -z "$par_cnv_map" ]; then
  extra_params+=( "--cnv-map $par_cnv_map" )
fi

if [ "$par_gvcf" = true ]; then
  extra_params+=( "--gvcf" )
fi

if [ ! -z "$par_gvcf_chunk" ]; then
  extra_params+=( "--gvcf-chunk $par_gvcf_chunk" )
fi

if [ ! -z "$par_variant_input" ]; then
  extra_params+=( "--variant-input $par_variant_input" )
fi

if [ "$par_only_use_input_alleles" = true ]; then
  extra_params+=( "--only-use-input-alleles")
fi
  
if [ ! -z "$par_haplotype_basis_alleles" ]; then
  extra_params+=( "--haplotype-basis-alleles $par_haplotype_basis_alleles" )
fi

if [ "$par_report_all_haplotype_alleles" = true ]; then
  extra_params+=( "--report-all-haplotype-alleles" )
fi

if [ "$par_report_monomorphic" = true ]; then
  extra_params+=( "--report-monomorphic" )
fi

if [ "$par_strict_vcf" = true ]; then
  extra_params+=( "--strict-vcf" )
fi

if [ "$par_pooled_discrete" = true ]; then
  extra_params+=( "--pooled-discrete" )
fi

if [ "$par_pooled_continuous" = true ]; then
  extra_params+=( "--pooled-continuous" )
fi

if [ "$par_use_reference_allele" = true ]; then
  extra_params+=( "--use-reference-allele" )
fi

if [ "$par_no_snps" = true ]; then
  extra_params+=( "--no-snps" )
fi

if [ "$par_no_snps" = true ]; then
  extra_params+=( "--no-snps" )
fi

if [ "$par_no_indels" = true ]; then
  extra_params+=( "--no-indels" )
fi

if [ "$par_no_mnps" = true ]; then
  extra_params+=( "--no-mnps" )
fi

if [ "$par_no_complex" = true ]; then
  extra_params+=( "--no-complex" )
fi

if [ "$par_no_partial_observations" = true ]; then
  extra_params+=( "--no-partial-observations" )
fi

if [ "$par_dont_left_align_indels" = true ]; then
  extra_params+=( "--dont-left-align-indels" )
fi

if [ "$par_use_duplicate_reads" = true ]; then
  extra_params+=( "--use-duplicate-reads" )
fi

if [ "$par_standard_filters" = true ]; then
  extra_params+=( "--standard-filters" )
fi

if [ "$par_no_population_priors" = true ]; then
  extra_params+=( "--no-population-priors" )
fi

if [ "$par_hwe_priors_off" = true ]; then
  extra_params+=( "--hwe-priors-off" )
fi

if [ "$par_binomial_obs_priors_off" = true ]; then
  extra_params+=( "--binomial-obs-priors-off" )
fi

if [ "$par_allele_balance_priors_off" = true ]; then
  extra_params+=( "--allele-balance-priors-off" )
fi

if [ "$par_legacy_gls" = true ]; then
  extra_params+=( "--legacy-gls" )
fi

if [ "$par_report_genotype_likelihood_max" = true ]; then
  extra_params+=( "--report-genotype-likelihood-max" )
fi

if [ "$par_exclude_unobserved_genotypes" = true ]; then
  extra_params+=( "--exclude-unobserved-genotypes" )
fi

if [ "$par_use_mapping_quality" = true ]; then
  extra_params+=( "--use-mapping-quality" )
fi

if [ "$par_harmonic_indel_quality" = true ]; then
  extra_params+=( "--harmonic-indel-quality" )
fi

if [ "$par_genotype_qualities" = true ]; then
  extra_params+=( "--genotype-qualities" )
fi

if [ "$par_debug" = true ]; then
  extra_params+=( "--debug" )
fi

if [ "$par_dd" = true ]; then
  extra_params+=( "-dd" )
fi

if [ ! -z "$par_observation_bias" ]; then
  extra_params+=( "--observation-bias $par_observation_bias" )
fi

if [ ! -z "$par_read_mismatch_limit" ]; then
  extra_params+=( "--read-mismatch-limit $par_read_mismatch_limit" )
fi

if [ ! -z "$par_read_snp_limit" ]; then
  extra_params+=( "--read-snp-limit $par_read_snp_limit" )
fi

if [ ! -z "$par_read_indel_limit" ]; then
  extra_params+=( "--read-indel-limit $par_read_indel_limit" )
fi

if [ ! -z "$par_base_quality_cap" ]; then
  extra_params+=( "--base-quality-cap $par_base_quality_cap" )
fi

if [ ! -z "$par_prob_contamination" ]; then
  extra_params+=( "--prob-contamination $par_prob_contamination" )
fi

if [ ! -z "$par_contamination_estimates" ]; then
  extra_params+=( "--contamination-estimates $par_contamination_estimates" )
fi

if [ ! -z "$par_genotype_variant_threshold" ]; then
  extra_params+=( "--genotype-variant-threshold $par_genotype_variant_threshold" )
fi

if [ ! -d "$par_output" ]; then
  mkdir $par_output
fi
     
freebayes --fasta-reference $par_fasta_reference --pvar $par_pvar --theta $par_theta \
          --ploidy $par_ploidy --min-repeat-entropy $par_min_repeat_entropy \
          --reference-quality $par_reference_quality --use-best-n-alleles $par_use_best_n_alleles \
          --max-complex-gap $par_max_complex_gap --min-repeat-size $par_min_repeat_size \
          --min-mapping-quality $par_min_mapping_quality --min-base-quality $par_min_base_quality \
          --min-supporting-allele-qsum $par_min_supporting_allele_qsum \
          --min-supporting-mapping-qsum $par_min_supporting_mapping_qsum \
          --mismatch-base-quality-threshold $par_mismatch_base_quality_threshold \
          --read-max-mismatch-fraction $par_read_max_mismatch_fraction \
          --min-alternate-fraction $par_min_alternate_fraction \
          --min-alternate-count $par_min_alternate_count \
          --min-alternate-qsum $par_min_alternate_qsum \
          --min-alternate-total $par_min_alternate_total --min-coverage $par_min_coverage \
          --genotyping-max-iterations $par_genotyping_max_iterations \
          --genotyping-max-banddepth $par_genotyping_max_banddepth \
          --posterior-integration-limits $par_posterior_integration_limits \
          --read-dependence-factor $par_read_dependence_factor  \
          --vcf ${par_output}${par_vcf} ${extra_params[@]}
