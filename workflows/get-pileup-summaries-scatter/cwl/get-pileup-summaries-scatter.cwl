class: Workflow
cwlVersion: v1.1
id: get-pileup-summaries-scatter

requirements:
- class: StepInputExpressionRequirement
- class: ScatterFeatureRequirement

inputs:
  jvm_mem: int?
  interval_files: File[]
  ref_fa:
    type: File?
    secondaryFiles: ['.fai', '^.dict']
  seq_file:
    type: File
    secondaryFiles: ['.bai?', '.crai?']
  variants:
    type: File
    secondaryFiles: ['.tbi']
  output_name: string

outputs:
  pileups_table_files:
    type: File[]
    outputSource: get-pileup-summaries/pileups_table

steps:
  get-pileup-summaries:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-get-pileup-summaries.4.1.3.0-1.1/tools/gatk-get-pileup-summaries/gatk-get-pileup-summaries.cwl
    scatter: intervals
    in:
      intervals: interval_files
      jvm_mem: jvm_mem
      ref_fa: ref_fa
      seq_file: seq_file
      variants: variants
      output_name: output_name
    out: [ pileups_table ]
