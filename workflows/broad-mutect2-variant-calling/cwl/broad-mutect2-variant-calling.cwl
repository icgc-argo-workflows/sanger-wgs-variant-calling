#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.1
id: broad-mutect2-variant-calling

requirements:
- class: StepInputExpressionRequirement
- class: ScatterFeatureRequirement

inputs:
  jvm_mem: int?
  ref_fa:
    type: File
    secondaryFiles: ['.fai', '^.dict']
  ref_dict: File
  intervals: File?
  scatter_count: int
  split_intervals_extra_args: string?
  unfiltered_vcf_name: string
  bam_output_name: string?
  germline_resource: File?
  pon: File?
  f1r2_tar_gz: string?
  m2_extra_args: string?
  variants_for_contamination:
    type: File
    secondaryFiles: ['.tbi']
  pileup_summary_name: string
  artifact_prior_table_name: string
  merge_vcfs_name: string
  merged_stats_name: string
  merged_pileup_name: string
  segmentation_name: string
  contamination_name: string
  filtered_vcf_name: string
  bwa_mem_index_image: File
  filtered_artifacts_vcf_name: string
  filtering_stats_name: string
  seq_format: string?
  library_strategy: string
  program: string
  donor_submitter_id: string
  normal_sample_submitter_id: string
  tumour_sample_submitter_id: string
  normal_specimen_type: string
  tumour_specimen_type: string
  object_store_endpoint_url: string
  bucket_name: string
  credentials_file: File
  payload_schema_version: string
  dna_alignment_bundle_type: string
  sequencing_experiment_bundle_type: string
  ssm_call_bundle_type: string

outputs:
  filtered_vcf:
    type: File
    outputSource: filter_alignment_artifacts/filtered_vcf
  filtering_stats:
    type: File
    outputSource: filter_mutect_calls/filtering_stats
  mutect_stats:
    type: File
    outputSource: merge_stats/merged_stats
  contamination_table:
    type: File
    outputSource: calculate_contamination/contamination_table
  segmentation_table:
    type: File
    outputSource: calculate_contamination/segmentation_table
  read_orientation_model_params:
    type: File
    outputSource: learn_read_orientation/artifact_prior_table

steps:

  get_payload_aligned_normal:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/ceph-get-payload.0.1.0/tools/ceph-get-payload/ceph-get-payload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: dna_alignment_bundle_type
      seq_format: seq_format
      library_strategy: library_strategy
      program: program
      donor_submitter_id: donor_submitter_id
      sample_submitter_id: normal_sample_submitter_id
      specimen_type: normal_specimen_type
    out: [ payload ]

  get_payload_aligned_tumour:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/ceph-get-payload.0.1.0/tools/ceph-get-payload/ceph-get-payload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: dna_alignment_bundle_type
      seq_format: seq_format
      library_strategy: library_strategy
      program: program
      donor_submitter_id: donor_submitter_id
      sample_submitter_id: tumour_sample_submitter_id
      specimen_type: tumour_specimen_type
    out: [ payload ]

  get_payload_tumour_sequencing_experiment:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/ceph-get-payload.0.1.0/tools/ceph-get-payload/ceph-get-payload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: sequencing_experiment_bundle_type
      library_strategy: library_strategy
      program: program
      donor_submitter_id: donor_submitter_id
      sample_submitter_id: tumour_sample_submitter_id
      specimen_type: tumour_specimen_type
    out: [ payload ]


  download_normal:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-download.0.1.0/tools/s3-download/s3-download.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      payload_json: get_payload_aligned_normal/payload
      s3_credential_file: credentials_file
    out: [ download_file ]

  download_tumour:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-download.0.1.0/tools/s3-download/s3-download.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      payload_json: get_payload_aligned_tumour/payload
      s3_credential_file: credentials_file
    out: [ download_file ]


  split_intervals:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-split-intervals.4.1.3.0-1.0/tools/gatk-split-intervals/gatk-split-intervals.cwl
    in:
      jvm_mem: jvm_mem
      ref_fa: ref_fa
      intervals: intervals
      scatter_count: scatter_count
      split_intervals_extra_args: split_intervals_extra_args
    out: [ interval_files ]

  mutect2_calling:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-mutect2.4.1.3.0-1.3/tools/gatk-mutect2/gatk-mutect2.cwl
    scatter: intervals
    in:
      jvm_mem: jvm_mem
      ref_fa: ref_fa
      tumour_reads: download_tumour/download_file
      normal_reads: download_normal/download_file
      output_vcf: unfiltered_vcf_name
      bam_output_name: bam_output_name
      germline_resource: germline_resource
      pon: pon
      intervals: split_intervals/interval_files
      f1r2_tar_gz: f1r2_tar_gz
      m2_extra_args: m2_extra_args
    out: [ unfiltered_vcf,  mutect_stats, bam_output, f1r2_counts]

  get_normal_pileup_summaries:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-get-pileup-summaries.4.1.3.0-1.1/tools/gatk-get-pileup-summaries/gatk-get-pileup-summaries.cwl
    scatter: intervals
    in:
      jvm_mem: jvm_mem
      ref_fa: ref_fa
      seq_file: download_normal/download_file
      variants: variants_for_contamination
      intervals: split_intervals/interval_files
      output_name: pileup_summary_name
    out: [ pileups_table ]

  get_tumour_pileup_summaries:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-get-pileup-summaries.4.1.3.0-1.1/tools/gatk-get-pileup-summaries/gatk-get-pileup-summaries.cwl
    scatter: intervals
    in:
      jvm_mem: jvm_mem
      ref_fa: ref_fa
      seq_file: download_tumour/download_file
      variants: variants_for_contamination
      intervals: split_intervals/interval_files
      output_name: pileup_summary_name
    out: [ pileups_table ]

  learn_read_orientation:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-learn-read-orientation-model.4.1.3.0-1.0/tools/gatk-learn-read-orientation-model/gatk-learn-read-orientation-model.cwl
    in:
      jvm_mem: jvm_mem
      input_f1r2_tar_gz: mutect2_calling/f1r2_counts
      output_name: artifact_prior_table_name
    out: [ artifact_prior_table ]

  merge_vcfs:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-merge-vcfs.4.1.3.0-1.1/tools/gatk-merge-vcfs/gatk-merge-vcfs.cwl
    in:
      jvm_mem: jvm_mem
      input_vcf: mutect2_calling/unfiltered_vcf
      output_name: merge_vcfs_name
    out: [ output_vcf ]

  merge_stats:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-merge-mutect-stats.4.1.3.0-1.0/tools/gatk-merge-mutect-stats/gatk-merge-mutect-stats.cwl
    in:
      jvm_mem: jvm_mem
      input_stats: mutect2_calling/mutect_stats
      output_name: merged_stats_name
    out: [ merged_stats ]

  merge_normal_pileups:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-gather-pileup-summaries.4.1.3.0-1.0/tools/gatk-gather-pileup-summaries/gatk-gather-pileup-summaries.cwl
    in:
      jvm_mem: jvm_mem
      ref_dict: ref_dict
      input_pileup: get_normal_pileup_summaries/pileups_table
      output_name: merged_pileup_name
    out: [ merged_pileup ]

  merge_tumour_pileups:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-gather-pileup-summaries.4.1.3.0-1.0/tools/gatk-gather-pileup-summaries/gatk-gather-pileup-summaries.cwl
    in:
      jvm_mem: jvm_mem
      ref_dict: ref_dict
      input_pileup: get_tumour_pileup_summaries/pileups_table
      output_name: merged_pileup_name
    out: [ merged_pileup ]

  calculate_contamination:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-calculate-contamination.4.1.3.0-1.0/tools/gatk-calculate-contamination/gatk-calculate-contamination.cwl
    in:
      jvm_mem: jvm_mem
      tumour_pileups: merge_tumour_pileups/merged_pileup
      normal_pileups: merge_normal_pileups/merged_pileup
      segmentation_output: segmentation_name
      contamination_output: contamination_name
    out: [ segmentation_table, contamination_table ]

  filter_mutect_calls:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-filter-mutect-calls.4.1.3.0-1.1/tools/gatk-filter-mutect-calls/gatk-filter-mutect-calls.cwl
    in:
      jvm_mem: jvm_mem
      unfiltered_vcf: merge_vcfs/output_vcf
      ref_fa: ref_fa
      output_vcf: filtered_vcf_name
      contamination_table: calculate_contamination/contamination_table
      segmentation_table: calculate_contamination/segmentation_table
      artifact_priors_tar_gz:
        source:
         - learn_read_orientation/artifact_prior_table
        linkMerge: merge_flattened
      mutect_stats: merge_stats/merged_stats
      filtering_stats_output: filtering_stats_name
    out: [ filtered_vcf, filtering_stats ]

  filter_alignment_artifacts:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-filter-alignment-artifacts.4.1.3.0-1.0/tools/gatk-filter-alignment-artifacts/gatk-filter-alignment-artifacts.cwl
    in:
      jvm_mem: jvm_mem
      bwa_mem_index_image: bwa_mem_index_image
      tumour_seq: download_tumour/download_file
      input_vcf: filter_mutect_calls/filtered_vcf
      output_vcf: filtered_artifacts_vcf_name
    out: [ filtered_vcf ]

  mutect2_ssm_payload_generate:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/payload-generation.0.1.3/tools/payload-generation/payload-generation.cwl
    in:
      bundle_type: ssm_call_bundle_type
      payload_schema_version: payload_schema_version
      file_to_upload: filter_alignment_artifacts/filtered_vcf
      input_metadata_aligned_seq:
        source:
          - get_payload_aligned_normal/payload
          - get_payload_aligned_tumour/payload
        linkMerge: merge_flattened
    out:
      [ payload ]

  mutect2_ssm_payload_s3_submit:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/payload-ceph-submission.0.1.4/tools/payload-ceph-submission/payload-ceph-submission.cwl
    in:
      metadata: get_payload_tumour_sequencing_experiment/payload
      payload: mutect2_ssm_payload_generate/payload
      credentials_file: credentials_file
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
    out:
      [ payload ]

  mutect2_ssm_s3_upload:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-upload.0.1.3/tools/s3-upload/s3-upload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: ssm_call_bundle_type
      upload_file: filter_alignment_artifacts/filtered_vcf
      payload_jsons:
        source:
         - mutect2_ssm_payload_s3_submit/payload
        linkMerge: merge_flattened
    out: []