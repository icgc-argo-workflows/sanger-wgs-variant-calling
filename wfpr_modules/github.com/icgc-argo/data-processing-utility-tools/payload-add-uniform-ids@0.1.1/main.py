#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Copyright (C) 2021,  Ontario Institute for Cancer Research

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Authors:
    Junjun Zhang
"""

import os
import sys
import argparse
import csv
import json


def get_id_mapping(id_mapping_tsv):
    id_mapping = dict()
    with open(id_mapping_tsv, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            type = row['type']
            submitter_id = row['submitter_id']
            uniform_id = row['uniform_id']
            if type in id_mapping:
                sys.exit(f"Values in 'type' field duplicated. Offending value: {type}, in file: {id_mapping_tsv}")
            else:
                id_mapping[type] = dict()

            if submitter_id in id_mapping[type]:
                sys.exit(f"Values in 'submitter_id' field duplicated. Offending value: {submitter_id}, for type: {type}, in file: {id_mapping_tsv}" )
            else:
                id_mapping[type][submitter_id] = uniform_id

    if 'donor' not in id_mapping or 'specimen' not in id_mapping or 'sample' not in id_mapping:
        sys.exit(f"Provided id_mapping_tsv file '{id_mapping_tsv}' is required to have ID mappings for 'donor', 'specimen' and 'sample'")

    return id_mapping


def add_uniform_ids(payload, id_mapping):
    samples = payload.pop('samples', [])
    if not samples:
        sys.exit("Error: no 'samples' found in the input payload JSON")

    updated_samples = []
    for sample in samples:
        if id_mapping['sample'].get(sample['submitterSampleId']):
            sample['sampleId'] = id_mapping['sample'][sample['submitterSampleId']]
        else:
            sys.exit(f"Provided id_mapping_tsv misses mapping for submitter sample ID: {sample['submitterSampleId']}")

        if id_mapping['specimen'].get(sample['specimen']['submitterSpecimenId']):
            sample['specimenId'] = id_mapping['specimen'][sample['specimen']['submitterSpecimenId']]
            sample['specimen']['specimenId'] = sample["specimenId"]
        else:
            sys.exit(f"Provided id_mapping_tsv misses mapping for submitter specimen ID: {sample['specimen']['submitterSpecimenId']}")

        if id_mapping['donor'].get(sample['donor']['submitterDonorId']):
            sample['donor']['donorId'] = id_mapping['donor'][sample['donor']['submitterDonorId']]
            sample['specimen']['donorId'] = sample['donor']['donorId']
        else:
            sys.exit(f"Provided id_mapping_tsv misses mapping for submitter donor ID: {sample['donor']['submitterDonorId']}")

        updated_samples.append(sample)

    payload['samples'] = updated_samples


def main():
    """
    Add uniform IDs for donor/specimen/sample to the original input payload JSON
    """

    parser = argparse.ArgumentParser(description='Tool: payload-add-uniform-ids')
    parser.add_argument('-p', '--payload-json', type=str,
                        help='Input payload JSON', required=True)
    parser.add_argument('-i', '--id-mapping-tsv', type=str,
                        help='TSV file containing mapping between submitter IDs and uniform IDs', required=True)
    parser.add_argument('-o', '--output-dir', type=str,
                        help='Output directory', required=True)
    args = parser.parse_args()

    if not os.path.isfile(args.payload_json):
        sys.exit('Error: specified input payload JSON %s does not exist or is not accessible!' % args.payload_json)

    if not os.path.isfile(args.id_mapping_tsv):
        sys.exit('Error: specified ID mapping TSV %s does not exist or is not accessible!' % args.id_mapping_tsv)

    if not os.path.isdir(args.output_dir):
        sys.exit('Error: specified output dir %s does not exist or is not accessible!' % args.output_dir)

    with open(args.payload_json, 'r') as p:
        payload = json.loads(p.read())

    id_mapping = get_id_mapping(args.id_mapping_tsv)

    add_uniform_ids(payload, id_mapping)

    output_payload_file = "%s.uniform_id_added.json" % os.path.splitext(os.path.basename(args.payload_json))[0]

    with open(os.path.join(args.output_dir, output_payload_file), 'w') as o:
        o.write(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
