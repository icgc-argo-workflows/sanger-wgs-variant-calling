#!/usr/bin/env nextflow

/*
 * Copyright (c) 2019, Ontario Institute for Cancer Research (OICR).
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
 * Authors:
 *   Junjun Zhang <junjun.zhang@oicr.on.ca>
 *   Linda Xiang <linda.xiang@oicr.on.ca>
 */

nextflow.preview.dsl=2
version = '0.1.1.0'

params.normal_analysis = ""
params.tumour_analysis = ""
params.qc_files = []
params.wf_name = ""
params.wf_short_name = ""
params.wf_version = ""
params.container_version = ""
params.cpus = 1
params.mem = 1  // GB


process payloadGenSangerQc {
  container "quay.io/icgc-argo/payload-gen-sanger-qc:payload-gen-sanger-qc.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path normal_analysis
    path tumour_analysis
    path qc_files
    val wf_name
    val wf_short_name
    val wf_version

  output:
    path "*.sanger_qc.payload.json", emit: payload
    path "out/*.tgz", emit: qc_files

  script:
    """
    payload-gen-sanger-qc.py \
      -n ${normal_analysis} \
      -t ${tumour_analysis} \
      -f ${qc_files} \
      -w ${wf_name} \
      -s ${wf_short_name} \
      -r ${workflow.runName} \
      -v ${wf_version}
    """
}
