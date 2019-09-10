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

outputs:
  run_params:
    type: File
    outputSource: sanger-calling/run_params
  result_archive:
    type: File
    outputSource: sanger-calling/result_archive
  timings:
    type: File
    outputSource: sanger-calling/timings
  global_time:
    type: File
    outputSource: sanger-calling/global_time

steps:
  generate-bas-normal:
    run: https://raw.githubusercontent.com/icgc-argo/variant-calling-tools/generate-bas.0.1.0/tools/generate-bas/generate-bas.cwl
    in:
      input: normal
    out: [ bam_and_bas ]
  generate-bas-tumour:
    run: https://raw.githubusercontent.com/icgc-argo/variant-calling-tools/generate-bas.0.1.0/tools/generate-bas/generate-bas.cwl
    in:
      input: tumour
    out: [ bam_and_bas ]
  sanger_calling:
    run: https://raw.githubusercontent.com/cancerit/dockstore-cgpwgs/2.1.0/cwls/cgpwgs.cwl
    in:
      reference: reference
      annot: annot
      snv_indel: snv_indel
      cnv_sv: cnv_sv
      qcset: qcset
      tumour: generate-bas-tumour/bam_and_bas
      tumourIdx: tumourIdx
      normal: generate-bas-normal/bam_and_bas
      normalIdx: normalIdx
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
