class: Workflow
cwlVersion: v1.1
id: sanger-variant-calling

requirements: []

inputs:
  reference: File
  annot: File
  snv_indel: File
  cnv_sv: File
  qcset: File
  tumour: File
  tumourIdx: File
  normal: File
  normalIdx: File
  exclude: string
  species: string?
  assembly: string?
  skipqc: boolean?
  pindelcpu: int?
  cavereads: int?
  purity: float?
  ploidy: float?
  num_threads: int?
  ref_file: File?
  endpoint_url: string
  bucket_name: string
  credentials_file: File
  normal_payload: File
  tumour_payload: File
  tumour_seq_exp_payload: File
  sanger_ssm_pattern: string
  ssm_bundle_type: string


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
  sanger_ssm:
    type: File
    outputSource: extract_sanger_ssm/output_file
  sanger_ssm_payload:
    type: File
    outputSource: sanger_ssm_payload_s3_submit/payload

steps:
  download_normal:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-download.0.1.0/tools/s3-download/s3-download.cwl
    in:
      endpoint_url: endpoint_url
      bucket_name: bucket_name
      payload_json: normal_payload
      s3_credential_file: credentials_file
    out: [ download_file, download_file_index ]

  download_tumour:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/s3-download.0.1.0/tools/s3-download/s3-download.cwl
    in:
      endpoint_url: endpoint_url
      bucket_name: bucket_name
      payload_json: tumour_payload
      s3_credential_file: credentials_file
    out: [ download_file, download_file_index ]


  generate_bas_normal:
    run: https://raw.githubusercontent.com/icgc-argo/variant-calling-tools/generate-bas.0.1.0/tools/generate-bas/generate-bas.cwl
    in:
      input: download_normal/download_file
      num_threads: num_threads
      ref_file: ref_file
    out: [ bam_and_bas ]

  generate_bas_tumour:
    run: https://raw.githubusercontent.com/icgc-argo/variant-calling-tools/generate-bas.0.1.0/tools/generate-bas/generate-bas.cwl
    in:
      input: download_tumour/download_file
      num_threads: num_threads
      ref_file: ref_file
    out: [ bam_and_bas ]

  sanger_calling:
    run: https://raw.githubusercontent.com/cancerit/dockstore-cgpwgs/2.1.0/cwls/cgpwgs.cwl
    in:
      reference: reference
      annot: annot
      snv_indel: snv_indel
      cnv_sv: cnv_sv
      qcset: qcset
      tumour: generate_bas_tumour/bam_and_bas
      tumourIdx: download_tumour/download_file_index
      normal: generate_bas_normal/bam_and_bas
      normalIdx: download_normal/download_file_index
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

  extract_sanger_ssm:
    run: https://raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/extract-files-from-tarball.0.1.0/tools/extract-files-from-tarball/extract-files-from-tarball.cwl
    in:
      tarball: repack_sanger_results/caveman
      pattern: sanger_ssm_pattern
    out:
      [ output_file ]

  sanger_ssm_payload_generate:
    run: https://raw.githubusercontent.com/icgc-argo/dna-seq-processing-tools/payload-generation.0.1.3/tools/payload-generation/payload-generation.cwl
    in:
      bundle_type: ssm_bundle_type
      payload_schema_version: payload_schema_version
      file_to_upload: extract_sanger_ssm/output_file
      input_metadata_aligned_seq:
        source:
          - normal_payload
          - tumour_payload
        linkMerge: merge_flattened
    out:
      [ payload ]

  sanger_ssm_payload_s3_submit:
    run: https://raw.githubusercontent.com/icgc-argo/dna-seq-processing-tools/payload-ceph-submission.0.1.4/tools/payload-ceph-submission/payload-ceph-submission.cwl
    in:
      metadata: tumour_seq_exp_payload
      payload: sanger_ssm_payload_generate/payload
      credentials_file: credentials_file
      endpoint_url: endpoint_url
      bucket_name: bucket_name
    out:
      [ payload ]

  sanger_ssm_s3_upload:
    run: https://raw.githubusercontent.com/icgc-argo/dna-seq-processing-tools/s3-upload.0.1.3/tools/s3-upload/s3-upload.cwl
    in:
      endpoint_url: endpoint_url
      bucket_name: bucket_name
      s3_credential_file: credentials_file
      bundle_type: ssm_bundle_type
      upload_file: extract-sanger-ssm/output_file
      payload_jsons:
        source:
         - sanger_ssm_payload_s3_submit/payload
        linkMerge: merge_flattened
