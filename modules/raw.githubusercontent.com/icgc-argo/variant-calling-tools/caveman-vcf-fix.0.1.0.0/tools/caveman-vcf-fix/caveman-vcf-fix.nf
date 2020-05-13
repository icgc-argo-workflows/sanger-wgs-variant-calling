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
version = '0.1.0.0'

params.input_tar = ""
params.container_version = ""
params.cpus = 2
params.mem = 1  // GB


process cavemanVcfFix {
  container "quay.io/icgc-argo/caveman-vcf-fix:caveman-vcf-fix.${params.container_version ?: version}"

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    path input_tar

  output:
    path "out/*.tgz", emit: fixed_tar

  script:
    """
    caveman-vcf-fix.py -i ${input_tar}
    """
}
