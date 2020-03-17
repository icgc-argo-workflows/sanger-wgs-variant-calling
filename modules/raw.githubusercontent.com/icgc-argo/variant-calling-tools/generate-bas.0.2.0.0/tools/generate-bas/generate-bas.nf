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
 * author Linda Xiang <linda.xiang@oicr.on.ca>
 */

nextflow.preview.dsl=2
version = '0.2.0.0'

params.reference = "NO_FILE"
params.seq = ""
params.seq_idx = ""
params.container_version = ""
params.cpus = 8
params.mem = 2  // GB


def getBasSecondaryFiles(main_file){  //this is kind of like CWL's secondary files
  def all_files = []
  all_files.add(main_file + '.fai')
  return all_files
}

process generateBas {
  container "quay.io/icgc-argo/generate-bas:generate-bas.${params.container_version ?: version}"

  cpus params.cpus
  memory "${params.mem} GB"

  tag "${seq.size()}"

  input:
    path seq
    path seq_idx
    path reference
    path reference_fai

  output:
    path "${seq.name}", emit: primary_file
    path "${seq.name}.bas", emit: bas_file
    path "${seq_idx.name}", emit: index_file

  script:
    arg_ref = reference.name != 'NO_FILE' ? "-r ${reference}" : ''
    """
    /opt/wtsi-cgp/bin/bam_stats \
      -i ${seq} \
      ${arg_ref} \
      --num_threads ${task.cpus} \
      -o ${seq.name}.bas
    """
}
