#!/usr/bin/env nextflow
nextflow.preview.dsl=2
name='sanger-wgs-variant-calling'
short_name='sanger-wgs'
version = '2.1.0-9.3.0'

/*
========================================================================================
                    ICGC ARGO Sanger WGS Variant Calling Workflow
========================================================================================
#### Homepage / Documentation
https://github.com/icgc-argo/sanger-wgs-variant-calling
#### Authors
Junjun Zhang @junjun-zhang <junjun.zhang@oicr.on.ca>
Linda Xiang @lindaxiang <linda.xiang@oicr.on.ca>
----------------------------------------------------------------------------------------

Required Parameters (no default):
--study_id                              song study ID
--tumour_aln_analysis_id                tumour WGS sequencing_alignment SONG analysis ID
--normal_aln_analysis_id                normal WGS sequencing_alignment SONG analysis ID
--ref_genome_fa                         reference genome '.fa' file, secondary file ('.fa.fai') is expected to be under the same folder
--ref_genome_tar                        tarball containing reference genome files from the same genome build
--vagrent_annot                         tarball containing VAGrENT annotation reference
--ref_cnv_sv_tar                        tarball containing CNV/SV reference
--ref_snv_indel_tar                     tarball containing SNV/Indel reference
--qcset_tar                             tarball containing QC Genotype reference
--song_url                              song server URL
--score_url                             score server URL
--api_token                             song/score API Token

General Parameters (with defaults):
--cpus                                  cpus given to all process containers (default 1)
--mem                                   memory (GB) given to all process containers (default 1)

Download Parameters (object):
--download
{
    song_container_version              song docker container version, defaults set below
    score_container_version             score docker container version, defaults set below
    song_url                            song url for download process (defaults to main song_url param)
    score_url                           score url for download process (defaults to main score_url param)
    api_token                           song/score API token for download process (defaults to main api_token param)
    song_cpus                           cpus for song container, defaults to cpus parameter
    song_mem                            memory (GB) for song container, defaults to mem parameter
    score_cpus                          cpus for score container, defaults to cpus parameter
    score_mem                           memory (GB) for score container, defaults to mem parameter
    score_transport_mem                 memory (GB) for score_transport, defaults to mem parameter
}

generateBas Parameters (object):
--generateBas
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for generateBas container, defaults to cpus parameter
    mem                                 memory (GB) for generateBas container, defaults to mem parameter
    ref_genome_fa                       reference genome '.fa' file, secondary file ('.fa.fai') is expected to be under the same folder
    tumour_normal                       Tumour or Normal
}

sangerWgsVariantCaller Parameters (object):
--sangerWgsVariantCaller
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for sangerWgsVariantCaller container, defaults to cpus parameter
    mem                                 memory (GB) for sangerWgsVariantCaller container, defaults to mem parameter
    species                             species name, default: 'human'
    assembly                            reference assembly, default: 'GRCh38'
    cavereads                           reads target per caveman split section, default: 800000
    exclude                             reference contigs to exclude, default: 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV'
    ploidy                              ploidy estimate of the genome, need to set purity if set ploidy
    purity                              purity estimate of the genome, need to set ploidy if set purity
    skipqc                              disable genotype, gender and verifyBamID, default: false
    skipannot                           disable annotation step for pindel and caveman results, default: true
    pindelcpu                           max CPUs for pindel analysis, >8 ignored, default: 6
    ref_genome_tar                      core reference tar.gz
    vagrent_annot                       VAGrENT*.tar.gz
    ref_snv_indel_tar                   SNV_INDEL*.tar.gz
    ref_cnv_sv_tar                      CNV_SV*.tar.gz
    qcset_tar                           qcGenotype*.tar.gz
}

repackSangerResults Parameters (object):
--repackSangerResults
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for repackSangerResults container, defaults to cpus parameter
    mem                                 memory (GB) for repackSangerResults container, defaults to mem parameter
    library_strategy                    library strategy of input reads, default: WGS
}

cavemanVcfFix Parameters (object):
--cavemanVcfFix
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for cavemanVcfFix container, defaults to cpus parameter
    mem                                 memory (GB) for cavemanVcfFix container, defaults to mem parameter
}

prepSangerSupplement Parameters (object):
--prepSangerSupplement
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for prepSangerSupplement container, defaults to cpus parameter
    mem                                 memory (GB) for prepSangerSupplement container, defaults to mem parameter
}

prepSangerQc Parameters (object):
--prepSangerQc
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for prepSangerQc container, defaults to cpus parameter
    mem                                 memory (GB) for prepSangerQc container, defaults to mem parameter
}


payloadGenVariantCall (object):
--payloadGenVariantCall
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for payloadGenVariantCall container, defaults to cpus parameter
    mem                                 memory (GB) for payloadGenVariantCall container, defaults to mem parameter
    wf_name                             workflow full name
    wf_short_name                       workflow short name embedded in results filename
    wf_version                          workflow version
}

Upload Parameters (object):
--upload
{
    song_container_version              song docker container version, defaults set below
    score_container_version             score docker container version, defaults set below
    song_url                            song url for upload process (defaults to main song_url param)
    score_url                           score url for upload process (defaults to main score_url param)
    api_token                           song/score API token for upload process (defaults to main api_token param)
    song_cpus                           cpus for song container, defaults to cpus parameter
    song_mem                            memory (GB) for song container, defaults to mem parameter
    score_cpus                          cpus for score container, defaults to cpus parameter
    score_mem                           memory (GB) for score container, defaults to mem parameter
    score_transport_mem                 memory (GB) for score_transport, defaults to mem parameter
    extract_cpus                        cpus for extract container, defaults to cpus parameter
    extract_mem                         memory (GB) extract score container, defaults to mem parameter
}

*/

params.study_id = ""
params.tumour_aln_analysis_id = ""
params.normal_aln_analysis_id = ""

// the following params if provided local files will be used
params.tumour_aln_metadata = "NO_FILE1"
params.tumour_aln_cram = "NO_FILE2"
params.tumour_extra_info = "NO_FILE3"
params.normal_aln_metadata = "NO_FILE4"
params.normal_aln_cram = "NO_FILE5"
params.normal_extra_info = "NO_FILE6"

// dir for outputs, must be set when running in local mode
params.publish_dir = ""
params.cleanup = true

params.api_token = ""
params.song_url = ""
params.score_url = ""

params.cpus = 1
params.mem = 1
params.max_retries = 5  // set to 0 will disable retry
params.first_retry_wait_time = 1  // in seconds

params.download = [:]
params.generateBas = [:]
params.sangerWgsVariantCaller = [:]
params.repackSangerResults = [:]
params.prepSangerSupplement = [:]
params.prepSangerQc = [:]
params.cavemanVcfFix = [:]
params.extractSangerCal = [:]
params.payloadGenVariantCall = [:]
params.upload = [:]
params.extractSangerCall = [:]

download_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'max_retries': params.max_retries,
    'first_retry_wait_time': params.first_retry_wait_time,
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    *:(params.download ?: [:])
]

generateBas_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'ref_genome_fa': '',
    *:(params.generateBas ?: [:])
]

sangerWgsVariantCaller_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'species': 'human',
    'assembly': 'GRCh38',
    'cavereads': 800000,
    'exclude': 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV',
    'skipqc': false,
    'skipannot': true,
    'pindelcpu': 6,
    'purity': null,
    'ploidy': null,
    'ref_genome_tar': '',
    'vagrent_annot': '',
    'ref_snv_indel_tar': '',
    'ref_cnv_sv_tar': '',
    'qcset_tar': '',
    *:(params.sangerWgsVariantCaller ?: [:])
]

repackSangerResults_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'library_strategy': 'WGS',
    *:(params.repackSangerResults ?: [:])
]

cavemanVcfFix_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.cavemanVcfFix ?: [:])
]

extractSangerCall_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.extractSangerCall ?: [:])
]

prepSangerSupplement_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.prepSangerSupplement ?: [:])
]

prepSangerQc_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.prepSangerQc ?: [:])
]

payloadGenVariantCall_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'publish_dir': params.publish_dir,
    *:(params.payloadGenVariantCall ?: [:])
]

upload_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'max_retries': params.max_retries,
    'first_retry_wait_time': params.first_retry_wait_time,
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    *:(params.upload ?: [:])
]


// Include all modules and pass params
include { SongScoreDownload as dnldT; SongScoreDownload as dnldN } from './wfpr_modules/github.com/icgc-argo/nextflow-data-processing-utility-tools/song-score-download@2.6.1/main.nf' params(download_params)
include { generateBas as basT; generateBas as basN; } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/generate-bas.0.2.1.0/tools/generate-bas/generate-bas' params(generateBas_params)
include { sangerWgsVariantCaller as sangerWgs } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/sanger-wgs-variant-caller.2.1.0-9/tools/sanger-wgs-variant-caller/sanger-wgs-variant-caller' params(sangerWgsVariantCaller_params)
include { repackSangerResults as repack } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/repack-sanger-results.0.2.0.0/tools/repack-sanger-results/repack-sanger-results' params(repackSangerResults_params)
include { cavemanVcfFix as cavemanFix } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/caveman-vcf-fix.0.1.0.0/tools/caveman-vcf-fix/caveman-vcf-fix' params(cavemanVcfFix_params)
include { prepSangerSupplement as prepSupp } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/prep-sanger-supplement.0.1.2.0/tools/prep-sanger-supplement/prep-sanger-supplement' params(prepSangerSupplement_params)
include { prepSangerQc as prepQc } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/prep-sanger-qc.0.1.2.0/tools/prep-sanger-qc/prep-sanger-qc' params(prepSangerQc_params)
include { extractFilesFromTarball as extractVarSnv; extractFilesFromTarball as extractVarIndel; extractFilesFromTarball as extractVarCnv; extractFilesFromTarball as extractVarSv } from './modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/extract-files-from-tarball.0.2.0.0/tools/extract-files-from-tarball/extract-files-from-tarball' params(extractSangerCall_params)
include { payloadGenVariantCalling as pGenVarSnv; payloadGenVariantCalling as pGenVarIndel; payloadGenVariantCalling as pGenVarCnv; payloadGenVariantCalling as pGenVarSv; payloadGenVariantCalling as pGenVarSupp; payloadGenVariantCalling as pGenQc } from "./modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/payload-gen-variant-calling.0.3.6.0/tools/payload-gen-variant-calling/payload-gen-variant-calling" params(payloadGenVariantCall_params)
include { SongScoreUpload as upSnv; SongScoreUpload as upIndel; SongScoreUpload as upCnv; SongScoreUpload as upSv; SongScoreUpload as upQc; SongScoreUpload as upSupp} from './wfpr_modules/github.com/icgc-argo/nextflow-data-processing-utility-tools/song-score-upload@2.6.0/main.nf' params(upload_params)
include { cleanupWorkdir as cleanup } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/cleanup-workdir@1.0.0/main'
include { getSecondaryFiles } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/helper-functions@1.0.0/main'
include { payloadAddUniformIds as pAddIdT; payloadAddUniformIds as pAddIdN } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/payload-add-uniform-ids@0.1.1/main'


workflow SangerWgs {
    take:
        study_id
        tumour_aln_analysis_id
        normal_aln_analysis_id
        tumour_aln_metadata
        tumour_aln_cram
        tumour_extra_info
        normal_aln_metadata
        normal_aln_cram
        normal_extra_info

    main:
        local_mode = false

        tumour_aln_seq = Channel.from()
        tumour_aln_seq_idx = Channel.from()
        normal_aln_seq = Channel.from()
        normal_aln_seq_idx = Channel.from()

        if (tumour_aln_analysis_id && normal_aln_analysis_id) {
            log.info "Run the workflow using input metadata / sequencing data from SONG/SCORE servers"

            // download tumour aligned seq and metadata from song/score (analysis type: sequencing_alignment)
            dnldT(study_id, tumour_aln_analysis_id)
            tumour_aln_seq = dnldT.out.files.flatten().first()
            tumour_aln_seq_idx = dnldT.out.files.flatten().last()
            tumour_aln_meta = dnldT.out.analysis_json

            // download normal aligned seq and metadata from song/score (analysis type: sequencing_alignment)
            dnldN(study_id, normal_aln_analysis_id)
            normal_aln_seq = dnldN.out.files.flatten().first()
            normal_aln_seq_idx = dnldN.out.files.flatten().last()
            normal_aln_meta = dnldN.out.analysis_json
        } else if (
            !tumour_aln_metadata.startsWith('NO_FILE') && \
            !tumour_aln_cram.startsWith('NO_FILE') && \
            !tumour_extra_info.startsWith('NO_FILE') && \
            !normal_aln_metadata.startsWith('NO_FILE') && \
            !normal_aln_cram.startsWith('NO_FILE') && \
            !normal_extra_info.startsWith('NO_FILE')
        ) {
            if (!params.publish_dir) {
                exit 1, "When use local inputs, params.publish_dir must be specified."
            } else {
                log.info "Use local inputs, outputs will be in: ${params.publish_dir}"
            }

            local_mode = true

            tumour_aln_seq = file(tumour_aln_cram)
            tumour_aln_seq_idx = Channel.fromPath(getSecondaryFiles(tumour_aln_cram, ['crai', 'bai']))
            pAddIdT(file(tumour_aln_metadata), file(tumour_extra_info))
            tumour_aln_meta = pAddIdT.out.payload

            normal_aln_seq = file(normal_aln_cram)
            normal_aln_seq_idx = Channel.fromPath(getSecondaryFiles(normal_aln_cram, ['crai', 'bai']))
            pAddIdN(file(normal_aln_metadata), file(normal_extra_info))
            normal_aln_meta = pAddIdN.out.payload
        } else {
            exit 1, "To download input aligned seq files from SONG/SCORE, please provide `params.tumour_aln_analysis_id` and `params.normal_aln_analysis_id`.\n" +
                "Or please provide `params.tumour_aln_metadata`, `params.tumour_aln_cram`, `params.tumour_extra_info`, `params.normal_aln_metadata`, `params.normal_aln_cram` and `params.normal_extra_info` to use local files as input."
        }

        // generate Bas for tumour
        basT(
            'tumour', tumour_aln_seq, tumour_aln_seq_idx.collect(),
            file(params.generateBas.ref_genome_fa), Channel.fromPath(getSecondaryFiles(params.generateBas.ref_genome_fa, ['fai']), checkIfExists: false).collect())

        // generate Bas for normal
        basN(
            'normal', normal_aln_seq, normal_aln_seq_idx.collect(),
            file(params.generateBas.ref_genome_fa), Channel.fromPath(getSecondaryFiles(params.generateBas.ref_genome_fa, ['fai']), checkIfExists: false).collect())


        sangerWgs(
            file(params.sangerWgsVariantCaller.ref_genome_tar),
            file(params.sangerWgsVariantCaller.vagrent_annot),
            file(params.sangerWgsVariantCaller.ref_snv_indel_tar),
            file(params.sangerWgsVariantCaller.ref_cnv_sv_tar),
            file(params.sangerWgsVariantCaller.qcset_tar),
            tumour_aln_seq,  // aln seq
            tumour_aln_seq_idx.collect(),  // idx
            basT.out.bas_file,  // bas
            normal_aln_seq,  // aln seq
            normal_aln_seq_idx.collect(),  // idx
            basN.out.bas_file  // bas
        )

        // repack results
        repack(sangerWgs.out.result_archive, 'WGS')

        // fix caveman bug
        cavemanFix(repack.out.caveman)

        // extract variant calls
        extractVarSnv(cavemanFix.out.fixed_tar, 'flagged.muts')
        extractVarIndel(repack.out.pindel, 'flagged')
        extractVarCnv(repack.out.ascat, 'copynumber.caveman')
        extractVarSv(repack.out.brass, 'annot')

        // prepare variant call supplements
        prepSupp(cavemanFix.out.fixed_tar.concat(
                    repack.out.pindel, repack.out.ascat, repack.out.brass, sangerWgs.out.timings).collect())

        pGenVarSnv(
            normal_aln_meta, tumour_aln_meta,
            extractVarSnv.out.extracted_files,
            name, short_name, version
        )
        pGenVarIndel(
            normal_aln_meta, tumour_aln_meta,
            extractVarIndel.out.extracted_files,
            name, short_name, version
        )
        pGenVarCnv(
            normal_aln_meta, tumour_aln_meta,
            extractVarCnv.out.extracted_files,
            name, short_name, version
        )
        pGenVarSv(
            normal_aln_meta, tumour_aln_meta,
            extractVarSv.out.extracted_files,
            name, short_name, version
        )
        pGenVarSupp(
            normal_aln_meta, tumour_aln_meta,
            prepSupp.out.supplement_tar.collect(),
            name, short_name, version
        )

        // prepare and upload sanger qc metrics
        prepQc(
            basN.out.bas_file_with_tn.concat(
                basT.out.bas_file_with_tn, repack.out.normal_contamination, repack.out.tumour_contamination,
                repack.out.genotyped, repack.out.ascat).collect()
        )
        pGenQc(
            normal_aln_meta, tumour_aln_meta,
            prepQc.out.qc_metrics_tar,
            name, short_name, version
        )

        // upload variant results
        if (!local_mode) {
            upSnv(study_id, pGenVarSnv.out.payload, pGenVarSnv.out.files_to_upload)
            upIndel(study_id, pGenVarIndel.out.payload, pGenVarIndel.out.files_to_upload)
            upCnv(study_id, pGenVarCnv.out.payload, pGenVarCnv.out.files_to_upload)
            upSv(study_id, pGenVarSv.out.payload, pGenVarSv.out.files_to_upload)
            upSupp(study_id, pGenVarSupp.out.payload, pGenVarSupp.out.files_to_upload)
            upQc(study_id, pGenQc.out.payload, pGenQc.out.files_to_upload)
        }

        if (params.cleanup) {
            if (local_mode) {
                cleanup(
                    basT.out.bas_file.concat(basN.out, sangerWgs.out, repack.out.pindel,
                        cavemanFix.out.fixed_tar, extractVarSnv.out.extracted_files, extractVarIndel.out.extracted_files,
                        extractVarCnv.out.extracted_files, extractVarSv.out.extracted_files,
                        prepQc.out, prepSupp.out, pGenVarSnv.out, pGenVarIndel.out,
                        pGenVarCnv.out, pGenVarSv.out, pGenQc.out, pGenVarSupp.out).collect(),
                    true
                )
            } else {
                cleanup(
                    dnldT.out.files.concat(dnldN.out, basT.out, basN.out, sangerWgs.out, repack.out.pindel,
                        cavemanFix.out.fixed_tar, extractVarSnv.out.extracted_files, extractVarIndel.out.extracted_files,
                        extractVarCnv.out.extracted_files, extractVarSv.out.extracted_files,
                        prepQc.out, prepSupp.out, pGenVarSnv.out, pGenVarIndel.out,
                        pGenVarCnv.out, pGenVarSv.out, pGenQc.out, pGenVarSupp.out).collect(),
                    upSv.out.analysis_id.concat(
                        upSnv.out.analysis_id, upIndel.out.analysis_id, upCnv.out.analysis_id,
                        upSupp.out.analysis_id, upQc.out.analysis_id).collect()
                )
            }
        }

}


workflow {
    SangerWgs(
        params.study_id,
        params.tumour_aln_analysis_id,
        params.normal_aln_analysis_id,
        params.tumour_aln_metadata,
        params.tumour_aln_cram,
        params.tumour_extra_info,
        params.normal_aln_metadata,
        params.normal_aln_cram,
        params.normal_extra_info
    )
}
