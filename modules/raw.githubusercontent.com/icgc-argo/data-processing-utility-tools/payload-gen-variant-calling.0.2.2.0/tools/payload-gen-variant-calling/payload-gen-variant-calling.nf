#!/usr/bin/env nextflow

/*
 * Copyright (c) 2019-2020, Ontario Institute for Cancer Research (OICR).
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
 * Author Junjun Zhang <junjun.zhang@oicr.on.ca>
 */

nextflow.preview.dsl = 2
version = '0.2.2.0'

params.normal_analysis = ""
params.tumour_analysis = ""
params.files_to_upload = []
params.wf_name = ""
params.wf_short_name = ""
params.wf_version = ""
params.container_version = ''
params.cpus = 1
params.mem = 1  // GB

process payloadGenVariantCalling {
  container "quay.io/icgc-argo/payload-gen-variant-calling:payload-gen-variant-calling.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path normal_analysis
    path tumour_analysis
    path files_to_upload
    val wf_name
    val wf_short_name
    val wf_version

  output:
    path "*.payload.json", emit: payload
    path "out/*{.tgz,.vcf.gz,.vcf.gz.tbi}", emit: files_to_upload

  script:
    args_tumour_analysis = !tumour_analysis.empty() ? "-t ${tumour_analysis}" : ""
    """
    payload-gen-variant-calling.py \
         -f ${files_to_upload} \
         -n ${normal_analysis} \
         -r ${workflow.runName} \
         -j ${workflow.sessionId} \
         -w ${wf_name} \
         -s ${wf_short_name} \
         -v ${wf_version} ${args_tumour_analysis}
    """
}
