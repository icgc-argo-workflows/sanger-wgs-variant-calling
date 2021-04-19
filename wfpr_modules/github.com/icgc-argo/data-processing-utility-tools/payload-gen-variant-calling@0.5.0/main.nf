#!/usr/bin/env nextflow

/*
 * Copyright (c) 2019-2021, Ontario Institute for Cancer Research (OICR).
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

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.5.0'

container = [
    'ghcr.io': 'ghcr.io/icgc-argo/data-processing-utility-tools.payload-gen-variant-calling'
]
default_container_registry = 'ghcr.io'
/********************************************************************/


// universal params go here
params.container_registry = ""
params.container_version = ""
params.container = ""

params.cpus = 1
params.mem = 1  // GB
params.publish_dir = ""  // set to empty string will disable publishDir

// tool specific parmas go here, add / change as needed
params.normal_analysis = ""
params.tumour_analysis = ""
params.files_to_upload = []
params.wf_name = ""
params.wf_short_name = ""
params.wf_version = ""


process payloadGenVariantCalling {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: "${params.publish_dir ? true : ''}"

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
    main.py \
         -f ${files_to_upload} \
         -n ${normal_analysis} \
         -r ${workflow.runName} \
         -j ${workflow.sessionId} \
         -w ${wf_name} \
         -s ${wf_short_name} \
         -v ${wf_version} ${args_tumour_analysis}
    """
}

// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  payloadGenVariantCalling(
    file(params.normal_analysis),
    file(params.tumour_analysis),
    Channel.fromPath(params.files_to_upload).collect(),
    params.wf_name,
    params.wf_short_name,
    params.wf_version
  )
}
