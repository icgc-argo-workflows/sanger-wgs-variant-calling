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
 * author Junjun Zhang <junjun.zhang@oicr.on.ca>
 *        Linda Xiang <linda.xiang@oicr.on.ca>
 */

nextflow.preview.dsl=2
version = '0.2.0.0'

params.sanger_results = ""
params.library_strategy = ""
params.container_version = ""
params.cpus = 1
params.mem = 2  // in GB


process repackSangerResults {
  container "quay.io/icgc-argo/repack-sanger-results:repack-sanger-results.${params.container_version ?: version}"
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path sanger_results
    val library_strategy

  output:
    path "*.normal.contamination.tgz", optional: true, emit: normal_contamination
    path "*.tumour.contamination.tgz", optional: true, emit: tumour_contamination
    path "*.ascat.tgz", optional: true, emit: ascat
    path "*.brass.tgz", optional: true, emit: brass
    path "*.caveman.tgz", optional: true, emit: caveman
    path "*.genotyped.tgz", optional: true, emit: genotyped
    path "*.pindel.tgz", optional: true, emit: pindel

  script:
    """
    repack-sanger-results.py \
      -i ${sanger_results} \
      -l ${library_strategy}
    """
}
