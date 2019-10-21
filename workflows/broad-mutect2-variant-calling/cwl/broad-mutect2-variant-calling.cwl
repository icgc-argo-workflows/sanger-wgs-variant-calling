#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.1
id: broad-mutect2-variant-calling

requirements:
- class: StepInputExpressionRequirement
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: MultipleInputFeatureRequirement

inputs:
    bam_output_name: string?
    bucket_name: string
    bwa_mem_index_image: File
    credentials_file: File
    submitter_donor_id: string
    germline_resource: File?
    intervals: File?
    jvm_mem: int?
    library_strategy: string
    m2_extra_args: string?
    normal_submitter_sample_id: string
    object_store_endpoint_url: string
    payload_schema_version: string
    pon: File?
    program: string
    ref_dict: File
    ref_fa: 
      secondaryFiles: 
        - .fai
        - ^.dict
      type: File
    scatter_count: int
    seq_format: string?
    split_intervals_extra_args: string?
    tumour_submitter_sample_id: string
    variants_for_contamination: 
      secondaryFiles: 
        - .tbi
      type: File
    wf_version: string

outputs:
  filtered_vcf:
    type: File[]
    outputSource: mutect2_ssm_payload_gen_and_s3_submit_wf/variant_call_renamed_result
  vcf_payload:
    type: File[]
    outputSource: mutect2_ssm_payload_gen_and_s3_submit_wf/payload
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
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/ceph-get-payload.0.1.2/tools/ceph-get-payload/ceph-get-payload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: { default: 'dna_alignment' }
      seq_format: seq_format
      library_strategy: library_strategy
      program_id: program
      submitter_donor_id: submitter_donor_id
      submitter_sample_id: normal_submitter_sample_id
      tumour_normal_designation: { default: 'normal' }
    out: [ payload ]

  get_payload_aligned_tumour:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/ceph-get-payload.0.1.2/tools/ceph-get-payload/ceph-get-payload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: { default: 'dna_alignment' }
      seq_format: seq_format
      library_strategy: library_strategy
      program_id: program
      submitter_donor_id: submitter_donor_id
      submitter_sample_id: tumour_submitter_sample_id
      tumour_normal_designation: { default: 'tumour' }
    out: [ payload ]

  get_payload_tumour_sequencing_experiment:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/ceph-get-payload.0.1.2/tools/ceph-get-payload/ceph-get-payload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: { default: 'sequencing_experiment' }
      library_strategy: library_strategy
      program_id: program
      submitter_donor_id: submitter_donor_id
      submitter_sample_id: tumour_submitter_sample_id
      tumour_normal_designation: { default: 'tumour' }
    out: [ payload ]


  download_normal:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-download.0.1.2/tools/s3-download/s3-download.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      payload_json: get_payload_aligned_normal/payload
      s3_credential_file: credentials_file
    out: [ download_file ]

  download_tumour:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-download.0.1.2/tools/s3-download/s3-download.cwl
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
      output_vcf: { default: 'broad-mutect2.unfiltered.vcf.gz' }
      bam_output_name: bam_output_name
      germline_resource: germline_resource
      pon: pon
      intervals: split_intervals/interval_files
      f1r2_tar_gz: { default: 'f1r2.tar.gz' }
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
      output_name: { default: 'normal_pileup_summary.tsv' }
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
      output_name: { default: 'tumour_pileup_summary.tsv' }
    out: [ pileups_table ]

  learn_read_orientation:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-learn-read-orientation-model.4.1.3.0-1.0/tools/gatk-learn-read-orientation-model/gatk-learn-read-orientation-model.cwl
    in:
      jvm_mem: jvm_mem
      input_f1r2_tar_gz: mutect2_calling/f1r2_counts
      output_name: { default: 'artifact-priors.tar.gz'}
    out: [ artifact_prior_table ]

  merge_vcfs:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-merge-vcfs.4.1.3.0-1.1/tools/gatk-merge-vcfs/gatk-merge-vcfs.cwl
    in:
      jvm_mem: jvm_mem
      input_vcf: mutect2_calling/unfiltered_vcf
      output_name: { default: 'merged.vcf.gz' }
    out: [ output_vcf ]

  merge_stats:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-merge-mutect-stats.4.1.3.0-1.0/tools/gatk-merge-mutect-stats/gatk-merge-mutect-stats.cwl
    in:
      jvm_mem: jvm_mem
      input_stats: mutect2_calling/mutect_stats
      output_name: { default: 'merged.stats' }
    out: [ merged_stats ]

  merge_normal_pileups:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-gather-pileup-summaries.4.1.3.0-1.0/tools/gatk-gather-pileup-summaries/gatk-gather-pileup-summaries.cwl
    in:
      jvm_mem: jvm_mem
      ref_dict: ref_dict
      input_pileup: get_normal_pileup_summaries/pileups_table
      output_name: { default: 'normal_pileup_merged.tsv' }
    out: [ merged_pileup ]

  merge_tumour_pileups:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-gather-pileup-summaries.4.1.3.0-1.0/tools/gatk-gather-pileup-summaries/gatk-gather-pileup-summaries.cwl
    in:
      jvm_mem: jvm_mem
      ref_dict: ref_dict
      input_pileup: get_tumour_pileup_summaries/pileups_table
      output_name: { default: 'tumour_pileup_merged.tsv' }
    out: [ merged_pileup ]

  calculate_contamination:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-calculate-contamination.4.1.3.0-1.0/tools/gatk-calculate-contamination/gatk-calculate-contamination.cwl
    in:
      jvm_mem: jvm_mem
      tumour_pileups: merge_tumour_pileups/merged_pileup
      normal_pileups: merge_normal_pileups/merged_pileup
      segmentation_output: { default: 'segments.table' }
      contamination_output: { default: 'contamination.table' }
    out: [ segmentation_table, contamination_table ]

  filter_mutect_calls:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-filter-mutect-calls.4.1.3.0-1.1/tools/gatk-filter-mutect-calls/gatk-filter-mutect-calls.cwl
    in:
      jvm_mem: jvm_mem
      unfiltered_vcf: merge_vcfs/output_vcf
      ref_fa: ref_fa
      output_vcf: { default: 'broad-mutect2.filtered.vcf.gz' }
      contamination_table: calculate_contamination/contamination_output
      segmentation_table: calculate_contamination/segmentation_output
      artifact_priors_tar_gz:
        source:
         - learn_read_orientation/artifact_prior_table
        linkMerge: merge_flattened
      mutect_stats: merge_stats/merged_stats
      filtering_stats_output: { default: 'filtering.stats.txt' }
    out: [ filtered_vcf, filtering_stats ]

  filter_alignment_artifacts:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-filter-alignment-artifacts.4.1.3.0-1.1/tools/gatk-filter-alignment-artifacts/gatk-filter-alignment-artifacts.cwl
    in:
      jvm_mem: jvm_mem
      bwa_mem_index_image: bwa_mem_index_image
      tumour_seq: download_tumour/download_file
      input_vcf: filter_mutect_calls/filtered_vcf
      output_vcf: { default: 'broad-mutect2.snv-indel.vcf.gz' }
    out: [ filtered_vcf ]

  mutect2_ssm_payload_gen_and_s3_submit_wf:
    run: https://raw.githubusercontent.com/icgc-argo/dna-seq-processing-wfs/payload-gen-and-s3-submit-wf.0.2.0/workflows/payload-gen-and-s3-submit-wf/cwl/payload-gen-and-s3-submit-wf.cwl
    in:
      bundle_type: { default: 'somatic_variant_call' }
      payload_schema_version: payload_schema_version
      files_to_upload:
        source:
          - filter_alignment_artifacts/filtered_vcf
        linkMerge: merge_flattened
      user_submit_metadata: get_payload_tumour_sequencing_experiment/payload
      analysis_input_payload:
        source:
          - get_payload_aligned_normal/payload
          - get_payload_aligned_tumour/payload
        linkMerge: merge_flattened
      wf_short_name: { default: 'broad-mutect2'}
      wf_version: wf_version
      credentials_file: credentials_file
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
    out:
      [ payload, variant_call_renamed_result ]


  mutect2_ssm_s3_upload:
     run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-upload.0.1.5/tools/s3-upload/s3-upload.cwl
     scatter: upload_file
     in:
       endpoint_url: object_store_endpoint_url
       bucket_name: bucket_name
       s3_credential_file: credentials_file
       bundle_type: { default: 'somatic_variant_call' }
       upload_file: mutect2_ssm_payload_gen_and_s3_submit_wf/variant_call_renamed_result
       payload_jsons: mutect2_ssm_payload_gen_and_s3_submit_wf/payload
     out: []