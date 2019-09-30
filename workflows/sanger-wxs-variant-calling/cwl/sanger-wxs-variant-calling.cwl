#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.1
id: sanger-wxs-variant-calling

requirements:
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement

inputs:
  reference: File
  annot: File
  snv_indel: File
  exclude: string
  species: string?
  assembly: string?
  num_threads: int?
  ref_file:
    type: File?
    secondaryFiles: ['.fai?']
  object_store_endpoint_url: string
  bucket_name: string
  credentials_file: File
  payload_schema_version: string
  sanger_ssm_vcf_name_pattern: string
  sanger_ssm_call_bundle_type: string
  dna_alignment_bundle_type: string
  sequencing_experiment_bundle_type: string
  seq_format: string?
  library_strategy: string
  program: string
  donor_submitter_id: string
  normal_sample_submitter_id: string
  tumour_sample_submitter_id: string
  normal_specimen_type: string
  tumour_specimen_type: string

outputs:
  run_params:
    type: File
    outputSource: sanger_calling/run_params
  timings:
    type: File
    outputSource: sanger_calling/timings
  sanger_ssm:
    type: File
    outputSource: extract_sanger_ssm/output_file
  sanger_ssm_payload:
    type: File
    outputSource: sanger_ssm_payload_s3_submit/payload

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


  generate_bas_normal:
    run: https://raw.githubusercontent.com/icgc-argo/variant-calling-tools/generate-bas.0.1.0/tools/generate-bas/generate-bas.cwl
    in:
      input: download_normal/download_file
      num_threads: num_threads
      ref_file: ref_file
    out: [ bam_and_bas, bai ]

  generate_bas_tumour:
    run: https://raw.githubusercontent.com/icgc-argo/variant-calling-tools/generate-bas.0.1.0/tools/generate-bas/generate-bas.cwl
    in:
      input: download_tumour/download_file
      num_threads: num_threads
      ref_file: ref_file
    out: [ bam_and_bas, bai ]

  sanger_calling:
    run: https://raw.githubusercontent.com/icgc-argo/variant-calling-tools/sanger-wxs-variant-caller.3.1.6-1/tools/sanger-wxs-variant-caller/sanger-wxs-variant-caller.cwl
    in:
      num_threads: num_threads
      reference: reference
      annot: annot
      snv_indel: snv_indel
      tumour: generate_bas_tumour/bam_and_bas
      tumourIdx: generate_bas_tumour/bai
      normal: generate_bas_normal/bam_and_bas
      normalIdx: generate_bas_normal/bai
      exclude: exclude
      species: species
      assembly: assembly
    out:
      - run_params
      - result_archive
      - timings

  repack_sanger_results:
    run: https://raw.githubusercontent.com/icgc-argo/variant-calling-tools/repack-sanger-results.0.1.1/tools/repack-sanger-results/repack-sanger-results.cwl
    in:
      input: sanger_calling/result_archive
      library_strategy: library_strategy
    out:
      - caveman
      - pindel

  extract_sanger_ssm:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/extract-files-from-tarball.0.1.0/tools/extract-files-from-tarball/extract-files-from-tarball.cwl
    in:
      tarball: repack_sanger_results/caveman
      pattern: sanger_ssm_vcf_name_pattern
    out:
      [ output_file ]

  sanger_ssm_payload_generate:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/payload-generation.0.1.3/tools/payload-generation/payload-generation.cwl
    in:
      bundle_type: sanger_ssm_call_bundle_type
      payload_schema_version: payload_schema_version
      file_to_upload: extract_sanger_ssm/output_file
      input_metadata_aligned_seq:
        source:
          - get_payload_aligned_normal/payload
          - get_payload_aligned_tumour/payload
        linkMerge: merge_flattened
    out:
      [ payload ]

  sanger_ssm_payload_s3_submit:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/payload-ceph-submission.0.1.4/tools/payload-ceph-submission/payload-ceph-submission.cwl
    in:
      metadata: get_payload_tumour_sequencing_experiment/payload
      payload: sanger_ssm_payload_generate/payload
      credentials_file: credentials_file
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
    out:
      [ payload ]

  sanger_ssm_s3_upload:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-upload.0.1.3/tools/s3-upload/s3-upload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: sanger_ssm_call_bundle_type
      upload_file: extract_sanger_ssm/output_file
      payload_jsons:
        source:
         - sanger_ssm_payload_s3_submit/payload
        linkMerge: merge_flattened
    out: []