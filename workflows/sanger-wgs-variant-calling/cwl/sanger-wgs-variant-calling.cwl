#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.1
id: sanger-wgs-variant-calling

requirements:
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement

inputs:
  reference: File
  annot: File
  snv_indel: File
  cnv_sv: File
  qcset: File
  exclude: string
  species: string?
  assembly: string?
  skipqc: boolean?
  pindelcpu: int?
  cavereads: int?
  purity: float?
  ploidy: float?
  num_threads: int?
  ref_file:
    type: File?
    secondaryFiles: ['.fai?']
  object_store_endpoint_url: string
  bucket_name: string
  credentials_file: File
  payload_schema_version: string
  seq_format: string?
  library_strategy: string
  program_id: string
  submitter_donor_id: string
  normal_submitter_sample_id: string
  tumour_submitter_sample_id: string

outputs:
  run_params:
    type: File
    outputSource: sanger_calling/run_params
  timings:
    type: File
    outputSource: sanger_calling/timings
  global_time:
    type: File
    outputSource: sanger_calling/global_time
  sanger_snv:
    type: File
    outputSource: extract_sanger_snv/output_file
  sanger_snv_payload:
    type: File[]
    outputSource: sanger_snv_payload_gen_and_s3_submit_wf/payload
  sanger_indel:
    type: File
    outputSource: extract_sanger_indel/output_file
  sanger_indel_payload:
    type: File[]
    outputSource: sanger_indel_payload_gen_and_s3_submit_wf/payload
  sanger_cnv:
    type: File
    outputSource: extract_sanger_cnv/output_file
  sanger_cnv_payload:
    type: File[]
    outputSource: sanger_cnv_payload_gen_and_s3_submit_wf/payload
  sanger_sv:
    type: File
    outputSource: extract_sanger_sv/output_file
  sanger_sv_payload:
    type: File[]
    outputSource: sanger_sv_payload_gen_and_s3_submit_wf/payload

steps:
  get_payload_aligned_normal:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/ceph-get-payload.0.1.1/tools/ceph-get-payload/ceph-get-payload.cwl
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
      specimen_type: { default: 'normal' }
    out: [ payload ]

  get_payload_aligned_tumour:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/ceph-get-payload.0.1.1/tools/ceph-get-payload/ceph-get-payload.cwl
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
      specimen_type: { default: 'tumour' }
    out: [ payload ]

  get_payload_tumour_sequencing_experiment:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/ceph-get-payload.0.1.1/tools/ceph-get-payload/ceph-get-payload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: { default: 'sequencing_experiment' }
      library_strategy: library_strategy
      program_id: program
      submitter_donor_id: submitter_donor_id
      submitter_sample_id: tumour_submitter_sample_id
      specimen_type: { default: 'tumour' }
    out: [ payload ]


  download_normal:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-download.0.1.1/tools/s3-download/s3-download.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      payload_json: get_payload_aligned_normal/payload
      s3_credential_file: credentials_file
    out: [ download_file ]

  download_tumour:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-download.0.1.1/tools/s3-download/s3-download.cwl
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
    run: https://raw.githubusercontent.com/icgc-argo/variant-calling-tools/sanger-wgs-variant-caller.2.1.0-2/tools/sanger-wgs-variant-caller/sanger-wgs-variant-caller.cwl
    in:
      num_threads: num_threads
      reference: reference
      annot: annot
      snv_indel: snv_indel
      cnv_sv: cnv_sv
      qcset: qcset
      tumour: generate_bas_tumour/bam_and_bas
      tumourIdx: generate_bas_tumour/bai
      normal: generate_bas_normal/bam_and_bas
      normalIdx: generate_bas_normal/bai
      exclude: exclude
      species: species
      assembly: assembly
      skipqc: skipqc
      pindelcpu: pindelcpu
      cavereads: cavereads
      purity: purity
      ploidy: ploidy
    out:
      - run_params
      - result_archive
      - timings
      - global_time

  repack_sanger_results:
    run: https://raw.githubusercontent.com/icgc-argo/variant-calling-tools/repack-sanger-results.0.1.0/tools/repack-sanger-results/repack-sanger-results.cwl
    in:
      input: sanger_calling/result_archive
    out:
      - normal_contamination
      - tumour_contamination
      - ascat
      - brass
      - caveman
      - genotyped
      - pindel

  extract_sanger_snv:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/extract-files-from-tarball.0.1.0/tools/extract-files-from-tarball/extract-files-from-tarball.cwl
    in:
      tarball: repack_sanger_results/caveman
      pattern: { default: 'flagged.muts.vcf.gz'}
    out:
      [ output_file ]

  sanger_snv_payload_gen_and_s3_submit_wf:
    run: https://raw.githubusercontent.com/icgc-argo/dna-seq-processing-wfs/payload-gen-and-s3-submit-wf.0.2.0/workflows/payload-gen-and-s3-submit-wf/cwl/payload-gen-and-s3-submit-wf.cwl
    in:
      bundle_type: { default: 'somatic_variant_call' }
      payload_schema_version: payload_schema_version
      files_to_upload:
        source:
          - extract_sanger_snv/output_file
        linkMerge: merge_flattened
      user_submit_metadata: get_payload_tumour_sequencing_experiment/payload
      analysis_input_payload:
        source:
          - get_payload_aligned_normal/payload
          - get_payload_aligned_tumour/payload
        linkMerge: merge_flattened
      wf_short_name: { default: 'sanger-wgs'}
      wf_version: wf_version
      data_type: { default: 'snv'}
      credentials_file: credentials_file
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
    out:
      [ payload ]

  sanger_snv_s3_upload:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-upload.0.1.4/tools/s3-upload/s3-upload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: { default: 'somatic_variant_call' }
      upload_file: extract_sanger_snv/output_file
      payload_jsons: sanger_snv_payload_gen_and_s3_submit_wf/payload
    out: []


  extract_sanger_indel:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/extract-files-from-tarball.0.1.0/tools/extract-files-from-tarball/extract-files-from-tarball.cwl
    in:
      tarball: repack_sanger_results/pindel
      pattern: { default: 'flagged.vcf.gz'}
    out:
      [ output_file ]


  sanger_indel_payload_gen_and_s3_submit_wf:
    run: https://raw.githubusercontent.com/icgc-argo/dna-seq-processing-wfs/payload-gen-and-s3-submit-wf.0.2.0/workflows/payload-gen-and-s3-submit-wf/cwl/payload-gen-and-s3-submit-wf.cwl
    in:
      bundle_type: { default: 'somatic_variant_call' }
      payload_schema_version: payload_schema_version
      files_to_upload:
        source:
          - extract_sanger_indel/output_file
        linkMerge: merge_flattened
      user_submit_metadata: get_payload_tumour_sequencing_experiment/payload
      analysis_input_payload:
        source:
          - get_payload_aligned_normal/payload
          - get_payload_aligned_tumour/payload
        linkMerge: merge_flattened
      wf_short_name: { default: 'sanger-wgs'}
      wf_version: wf_version
      data_type: { default: 'indel'}
      credentials_file: credentials_file
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
    out:
      [ payload ]


  sanger_indel_s3_upload:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-upload.0.1.4/tools/s3-upload/s3-upload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: { default: 'somatic_variant_call' }
      upload_file: extract_sanger_indel/output_file
      payload_jsons: sanger_indel_payload_gen_and_s3_submit_wf/payload
    out: []

  extract_sanger_cnv:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/extract-files-from-tarball.0.1.0/tools/extract-files-from-tarball/extract-files-from-tarball.cwl
    in:
      tarball: repack_sanger_results/ascat
      pattern: { default: 'copynumber.caveman.vcf.gz'}
    out:
      [ output_file ]


  sanger_cnv_payload_gen_and_s3_submit_wf:
    run: https://raw.githubusercontent.com/icgc-argo/dna-seq-processing-wfs/payload-gen-and-s3-submit-wf.0.2.0/workflows/payload-gen-and-s3-submit-wf/cwl/payload-gen-and-s3-submit-wf.cwl
    in:
      bundle_type: { default: 'somatic_variant_call' }
      payload_schema_version: payload_schema_version
      files_to_upload:
        source:
          - extract_sanger_cnv/output_file
        linkMerge: merge_flattened
      user_submit_metadata: get_payload_tumour_sequencing_experiment/payload
      analysis_input_payload:
        source:
          - get_payload_aligned_normal/payload
          - get_payload_aligned_tumour/payload
        linkMerge: merge_flattened
      wf_short_name: { default: 'sanger-wgs'}
      wf_version: wf_version
      data_type: { default: 'cnv'}
      credentials_file: credentials_file
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
    out:
      [ payload ]


  sanger_cnv_s3_upload:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-upload.0.1.4/tools/s3-upload/s3-upload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: { default: 'somatic_variant_call' }
      upload_file: extract_sanger_cnv/output_file
      payload_jsons: sanger_cnv_payload_gen_and_s3_submit_wf/payload
    out: []


  extract_sanger_sv:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/extract-files-from-tarball.0.1.0/tools/extract-files-from-tarball/extract-files-from-tarball.cwl
    in:
      tarball: repack_sanger_results/brass
      pattern: { default: 'annot.vcf.gz'}
    out:
      [ output_file ]


  sanger_sv_payload_gen_and_s3_submit_wf:
    run: https://raw.githubusercontent.com/icgc-argo/dna-seq-processing-wfs/payload-gen-and-s3-submit-wf.0.2.0/workflows/payload-gen-and-s3-submit-wf/cwl/payload-gen-and-s3-submit-wf.cwl
    in:
      bundle_type: { default: 'somatic_variant_call' }
      payload_schema_version: payload_schema_version
      files_to_upload:
        source:
          - extract_sanger_sv/output_file
        linkMerge: merge_flattened
      user_submit_metadata: get_payload_tumour_sequencing_experiment/payload
      analysis_input_payload:
        source:
          - get_payload_aligned_normal/payload
          - get_payload_aligned_tumour/payload
        linkMerge: merge_flattened
      wf_short_name: { default: 'sanger-wgs'}
      wf_version: wf_version
      data_type: { default: 'sv'}
      credentials_file: credentials_file
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
    out:
      [ payload ]


  sanger_sv_s3_upload:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-upload.0.1.4/tools/s3-upload/s3-upload.cwl
    in:
      endpoint_url: object_store_endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: { default: 'somatic_variant_call' }
      upload_file: extract_sanger_sv/output_file
      payload_jsons: sanger_sv_payload_gen_and_s3_submit_wf/payload
    out: []
