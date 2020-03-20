#!/bin/bash nextflow

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
 * author Junjun Zhang <junjun.zhang@oicr.on.ca>
 */

nextflow.preview.dsl=2

params.tarball = ""
params.pattern = ""

process extractFilesFromTarball {
  container 'ubuntu:18.04'

  input:
    path tarball
    val pattern

  output:
    path "*${pattern}{.bam,.cram,.vcf.gz}", emit: output_file
    path "*${pattern}{.bam.bai,.cram.crai,.vcf.gz.tbi}", emit: output_file_index
    tuple path("*${pattern}{.bam,.cram,.vcf.gz}"), path("*${pattern}{.bam.bai,.cram.crai,.vcf.gz.tbi}"), emit: extracted_files

  script:
    """
    tar -xzf ${tarball}
    """
}
