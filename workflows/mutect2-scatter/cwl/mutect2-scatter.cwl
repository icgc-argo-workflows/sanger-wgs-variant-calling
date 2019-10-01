class: Workflow
cwlVersion: v1.1
id: mutect2-scatter

requirements:
- class: StepInputExpressionRequirement
- class: ScatterFeatureRequirement

inputs:
  jvm_mem: int?
  interval_files: File[]
  ref_fa:
    type: File
    secondaryFiles: ['.fai', '^.dict']
  tumour_reads:
    type: File
    secondaryFiles: ['.bai?', '.crai?']
  normal_reads:
    type: File?
    secondaryFiles: ['.bai?', '.crai?']
  output_vcf: string
  germline_resource: File?
  pon:
    type: File?
    secondaryFiles: ['.idx?', '.tbi?']
  f1r2_tar_gz: string?

outputs:
  unfiltered_vcfs:
    type: File[]
    outputSource: mutect2/unfiltered_vcf
  mutect_stats_files:
    type: File[]
    outputSource: mutect2/mutect_stats
  bam_output_files:
    type: File[]?
    outputSource: mutect2/bam_output
  f1r2_counts_files:
    type: File[]?
    outputSource: mutect2/f1r2_counts

steps:
  mutect2:
    run: https://raw.githubusercontent.com/icgc-argo/gatk-tools/gatk-mutect2.4.1.3.0-1.3/tools/gatk-mutect2/gatk-mutect2.cwl
    scatter: intervals
    in:
      intervals: interval_files
      jvm_mem: jvm_mem
      ref_fa: ref_fa
      tumour_reads: tumour_reads
      normal_reads: normal_reads
      output_vcf: output_vcf
      germline_resource: germline_resource
      pon: pon
      f1r2_tar_gz: f1r2_tar_gz
    out:
    - unfiltered_vcf
    - mutect_stats
    - bam_output
    - f1r2_counts
