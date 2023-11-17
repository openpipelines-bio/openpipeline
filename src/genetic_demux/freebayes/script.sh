#!/bin/bash
set -eo pipefail

# Unset boolean flags if their values are not 'true'
for flag in par_stdin par_gvcf par_only_use_input_alleles par_report_all_haplotype_alleles par_report_monomorphic par_strict_vcf \
             par_pooled_discrete par_pooled_continuous par_use_reference_allele par_throw_away_snp_obs par_throw_away_indel_obs par_throw_away_mnps_obs par_throw_away_complex_obs \
             par_no_partial_observations par_dont_left_align_indels par_use_duplicate_reads par_standard_filters par_no_population_priors \
             par_hwe_priors_off par_binomial_obs_priors_off par_allele_balance_priors_off par_legacy_gls par_report_genotype_likelihood_max \
             par_exclude_unobserved_genotypes par_use_mapping_quality par_harmonic_indel_quality par_genotype_qualities par_debug par_dd; do
  [[ "${!flag}" != "true" ]] && unset "$flag"
done

# Create output directory if it doesn't exist
if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

freebayes \
  --fasta-reference $par_fasta_reference \
  --pvar $par_pvar \
  --theta $par_theta \
  --ploidy $par_ploidy \
  --min-repeat-entropy $par_min_repeat_entropy \
  --reference-quality $par_reference_quality \
  --use-best-n-alleles $par_use_best_n_alleles \
  --max-complex-gap $par_max_complex_gap \
  --min-repeat-size $par_min_repeat_size \
  --min-mapping-quality $par_min_mapping_quality \
  --min-base-quality $par_min_base_quality \
  --min-supporting-allele-qsum $par_min_supporting_allele_qsum \
  --min-supporting-mapping-qsum $par_min_supporting_mapping_qsum \
  --mismatch-base-quality-threshold $par_mismatch_base_quality_threshold \
  --read-max-mismatch-fraction $par_read_max_mismatch_fraction \
  --min-alternate-fraction $par_min_alternate_fraction \
  --min-alternate-count $par_min_alternate_count \
  --min-alternate-qsum $par_min_alternate_qsum \
  --min-alternate-total $par_min_alternate_total \
  --min-coverage $par_min_coverage \
  --genotyping-max-iterations $par_genotyping_max_iterations \
  --genotyping-max-banddepth $par_genotyping_max_banddepth \
  --posterior-integration-limits $par_posterior_integration_limits \
  --read-dependence-factor $par_read_dependence_factor \
  --vcf ${par_output}/${par_vcf} \
  ${par_bam:+--bam $par_bam} \
  ${par_bam_list:+--bam-list $par_bam_list} \
  ${par_stdin:+--stdin} \
  ${par_targets:+--targets $par_targets} \
  ${par_region:+--region $par_region} \
  ${par_max_coverage:+--max-coverage $par_max_coverage} \
  ${par_samples:+--samples $par_samples} \
  ${par_populations:+--populations $par_populations} \
  ${par_cnv_map:+--cnv-map $par_cnv_map} \
  ${par_gvcf:+--gvcf} \
  ${par_gvcf_chunk:+--gvcf-chunk $par_gvcf_chunk} \
  ${par_variant_input:+--variant-input $par_variant_input} \
  ${par_only_use_input_alleles:+--only-use-input-alleles} \
  ${par_haplotype_basis_alleles:+--haplotype-basis-alleles $par_haplotype_basis_alleles} \
  ${par_report_all_haplotype_alleles:+--report-all-haplotype-alleles} \
  ${par_report_monomorphic:+--report-monomorphic} \
  ${par_strict_vcf:+--strict-vcf} \
  ${par_pooled_discrete:+--pooled-discrete} \
  ${par_pooled_continuous:+--pooled-continuous} \
  ${par_use_reference_allele:+--use-reference-allele} \
  ${par_throw_away_snp_obs:+--throw-away-snp-obs} \
  ${par_throw_away_indel_obs:+--throw-away-indel-obs} \
  ${par_throw_away_mnps_obs:+--throw-away-mnps-obs} \
  ${par_throw_away_complex_obs:+--throw-away-complex-obs} \
  ${par_no_partial_observations:+--no-partial-observations} \
  ${par_dont_left_align_indels:+--dont-left-align-indels} \
  ${par_use_duplicate_reads:+--use-duplicate-reads} \
  ${par_standard_filters:+--standard-filters} \
  ${par_no_population_priors:+--no-population-priors} \
  ${par_hwe_priors_off:+--hwe-priors-off} \
  ${par_binomial_obs_priors_off:+--binomial-obs-priors-off} \
  ${par_allele_balance_priors_off:+--allele-balance-priors-off} \
  ${par_legacy_gls:+--legacy-gls} \
  ${par_report_genotype_likelihood_max:+--report-genotype-likelihood-max} \
  ${par_exclude_unobserved_genotypes:+--exclude-unobserved-genotypes} \
  ${par_use_mapping_quality:+--use-mapping-quality} \
  ${par_harmonic_indel_quality:+--harmonic-indel-quality} \
  ${par_genotype_qualities:+--genotype-qualities} \
  ${par_debug:+--debug} \
  ${par_dd:+-dd} \
  ${par_observation_bias:+--observation-bias $par_observation_bias} \
  ${par_read_mismatch_limit:+--read-mismatch-limit $par_read_mismatch_limit} \
  ${par_read_snp_limit:+--read-snp-limit $par_read_snp_limit} \
  ${par_read_indel_limit:+--read-indel-limit $par_read_indel_limit} \
  ${par_base_quality_cap:+--base-quality-cap $par_base_quality_cap} \
  ${par_prob_contamination:+--prob-contamination $par_prob_contamination} \
  ${par_contamination_estimates:+--contamination-estimates $par_contamination_estimates} \
  ${par_genotype_variant_threshold:+--genotype-variant-threshold $par_genotype_variant_threshold}
