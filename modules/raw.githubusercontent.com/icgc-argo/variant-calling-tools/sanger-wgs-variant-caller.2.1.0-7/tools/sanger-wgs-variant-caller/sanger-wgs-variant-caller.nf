#!/usr/bin/env nextflow

/*
 * Copyright (c) 2020, Ontario Institute for Cancer Research (OICR).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/*
 * author Junjun Zhang <junjun.zhang@oicr.on.ca>
 *        Linda Xiang <linda.xiang@oicr.on.ca>
 */

nextflow.preview.dsl=2
version = '2.1.0-7'

params.ref_genome_tar = ""
params.vagrent_annot = ""
params.ref_snv_indel_tar = ""
params.ref_cnv_sv_tar = ""
params.qcset_tar = ""
params.tumour = ""
params.tumourIdx = ""
params.normal = ""
params.normalIdx = ""
params.exclude = "chrUn%,HLA%,%_alt,%_random,chrM,chrEBV"
params.species = "human"
params.assembly = "GRCh38"
params.skipqc = false
params.skipannot = false
params.pindelcpu = 8
params.cavereads = 350000
params.purity = ""
params.ploidy = ""
params.container_version = ""
params.cpus = 24
params.mem = 128  // GB


def getSangerWgsSecondaryFiles(main_file){  //this is kind of like CWL's secondary files
  def all_files = []
  for (ext in ['.bas']) {
    all_files.add(main_file + ext)
  }
  return all_files
}

process sangerWgsVariantCaller {
  container "quay.io/icgc-argo/sanger-wgs-variant-caller:sanger-wgs-variant-caller.${params.container_version ?: version}"

  cpus params.cpus
  memory "${params.mem} GB"

  tag "${tumour.size()}"

  input:
    path reference
    path annot
    path snv_indel
    path cnv_sv
    path qcset
    path tumour
    path tidx
    path tumour_bas
    path normal
    path nidx
    path normal_bas

  output:
    path "run.params", emit: run_params
    path "WGS_*_vs_*.result.tar.gz", emit: result_archive
    path "WGS_*_vs_*.timings.tar.gz", emit: timings
    path "WGS_*_vs_*.time", emit: global_time

  script:
    arg_skipqc = params.skipqc ? "-skipqc" : ""
    arg_skipannot = params.skipannot ? "-skipannot" : ""
    arg_pu = params.pu ? "-pu ${params.pu}" : ""
    arg_pi = params.pi ? "-pu ${params.pi}" : ""
    """
    /opt/wtsi-cgp/bin/ds-cgpwgs.pl \
      -cores ${task.cpus} \
      -reference ${reference} \
      -annot ${annot} \
      -snv_indel ${snv_indel} \
      -cnv_sv ${cnv_sv} \
      -qcset ${qcset} \
      -tumour ${tumour} \
      -tidx ${tidx} \
      -normal ${normal} \
      -nidx ${nidx} \
      -exclude ${params.exclude} \
      -species ${params.species} \
      -assembly ${params.assembly} \
      ${arg_skipqc} \
      ${arg_skipannot} \
      -pindelcpu ${params.pindelcpu} \
      -cavereads ${params.cavereads} \
      ${arg_pu} \
      ${arg_pi} \
      -outdir \$PWD
    """
}
