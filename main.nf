#!/usr/bin/env nextflow
nextflow.preview.dsl=2


/*
========================================================================================
                    ICGC-ARGO Sanger WGS Variant Calling Pipeline
========================================================================================
#### Homepage / Documentation
https://github.com/icgc-argo/sanger-wgs-variant-calling
#### Authors
Junjun Zhang @junjun-zhang <junjun.zhang@oicr.on.ca>
Linda Xiang @lindaxiang <linda.xiang@oicr.on.ca>
----------------------------------------------------------------------------------------

Required Parameters (no default):
--study_id                              SONG study ID
--tumour_aln_analysis_id                Tumour WGS sequencing_alignment SONG analysis ID
--normal_aln_analysis_id                Normal WGS sequencing_alignment SONG analysis ID
--ref_genome_fa                         Reference genome '.fa' file, secondary file ('.fa.fai') is expected to be under the same folder
--sanger_ref_genome_tar                 Tarball containing reference genome files from the same genome build
--sanger_vagrent_annot                  Tarball containing VAGrENT annotation reference
--sanger_ref_cnv_sv_tar                 Tarball containing CNV/SV reference
--sanger_ref_snv_indel_tar              Tarball containing SNV/Indel reference
--sanger_qcset_tar                      Tarball containing QC Genotype reference
--song_url                              SONG server URL
--score_url                             SCORE server URL
--api_token                             SONG/SCORE API Token

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
    song_cpus
    song_mem
    score_cpus
    score_mem
    score_transport_mem
}

generateBas Parameters (object):
--generateBas
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for seqDataToLaneBam container, defaults to cpus parameter
    mem                                 memory (GB) for seqDataToLaneBam container, defaults to mem parameter
}

sangerWgsVariantCaller Parameters (object):
--sangerWgsVariantCaller
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for bwaMemAligner container, defaults to cpus parameter
    mem                                 memory (GB) for bwaMemAligner container, defaults to mem parameter
    exclude                             reference contigs to exclude, default: 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV'
    ploidy                              ploidy estimate of the genome, default: 2.0
    purity                              purity estimate of the genome, default: 1.0
    seq_format                          input aligned sequence format: default: cram
}

repackSangerResults Parameters (object):
--repackSangerResults
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for bamMergeSortMarkdup container, defaults to cpus parameter
    mem                                 memory (GB) for bamMergeSortMarkdup container, defaults to mem parameter
    library_strategy                    library strategy of input reads, default: WGS
}

payloadGenVariantCall (object):
--payloadGenVariantCall
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for align container, defaults to cpus parameter
    mem                                 memory (GB) for align container, defaults to mem parameter
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
params.ref_genome_fa = ""
params.sanger_ref_genome_tar = ""
params.sanger_vagrent_annot = ""
params.sanger_ref_cnv_sv_tar = ""
params.sanger_ref_snv_indel_tar = ""
params.sanger_qcset_tar = ""
params.cleanup = true

params.cpus = 1
params.mem = 1
params.pindelcpu = 1

download_params = [
    'song_container_version': '4.0.0',
    'score_container_version': '3.0.1',
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    *:(params.download ?: [:])
]

generateBas_params = [
    *:(params.generateBas ?: [:])
]

sangerWgsVariantCaller_params = [
    'species': 'human',
    'assembly': 'GRCh38',
    'cavereads': 800000,
    'exclude': 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV',
    'ploidy': 2.0,
    'purity': 1.0,
    'seq_format': 'cram',
    *:(params.sangerWgsVariantCaller ?: [:])
]

repackSangerResults_params = [
    'library_strategy': 'WGS',
    *:(params.repackSangerResults ?: [:])
]

extractSangerCall_params = [
    *:(params.extractSangerCall ?: [:])
]

payloadGenVariantCall_params = [
    *:(params.payloadGenVariantCall ?: [:])
]

payloadGenQcMetrics_params = [
    *:(params.payloadGenVariantCall ?: [:])
]

upload_params = [
    'song_container_version': '4.0.0',
    'score_container_version': '3.0.1',
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    *:(params.uploadVariant ?: [:])
]


// Include all modules and pass params
include { songScoreDownload as dnldT; songScoreDownload as dnldN } from './song-score-utils/song-score-download' params(download_params)
include { generateBas as basT; generateBas as basN; } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/generate-bas.0.1.1/tools/generate-bas/generate-bas' params(generateBas_params)
include sangerWgsVariantCaller as sangerWgs from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/sanger-wgs-variant-caller.2.1.0-4/tools/sanger-wgs-variant-caller/sanger-wgs-variant-caller' params(sangerWgsVariantCaller_params)
include repackSangerResults as repack from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/repack-sanger-results.0.1.2/tools/repack-sanger-results/repack-sanger-results' params(repackSangerResults_params)
include { extractFilesFromTarball as extractVar; extractFilesFromTarball as extractQC } from './modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/extract-files-from-tarball.0.1.1/tools/extract-files-from-tarball/extract-files-from-tarball' params(extractSangerCall_params)
include { payloadGenSangerVariant as pGenVar } from "./modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/payload-gen-dna-seq-qc.0.2.0.0/tools/payload-gen-dna-seq-qc/payload-gen-dna-seq-qc.nf" params(payloadGenVariantCall_params)
include { payloadGenSangerQC as pGenQC } from "./modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/payload-gen-dna-seq-qc.0.2.0.0/tools/payload-gen-dna-seq-qc/payload-gen-dna-seq-qc.nf" params(payloadGenQcMetrics_params)
include { songScoreUpload as upVar; songScoreUpload as upQC} from './song-score-utils/song-score-upload' params(upload_params)


process cleanup {
    container "ubuntu:18.04"

    input:
        path files_to_delete
        val aligned_seq_analysis_id
        val qc_metrics_analysis_id

    script:
        """
        IFS=" "
        read -a files <<< "${files_to_delete}"
        for f in "\${files[@]}"
            do rm -fr \$(dirname \$(readlink -f \$f))/*  # delete all files and subdirs but not hidden ones
        done
        """
}


workflow SangerWgs {
    take:
        study_id
        tumour_aln_analysis_id
        normal_aln_analysis_id
        ref_genome_fa

    main:
        // download tumour aligned seq and metadata from song/score (analysis type: sequencing_alignment)
        dnldT(study_id, normal_aln_analysis_id)

        // download normal aligned seq and metadata from song/score (analysis type: sequencing_alignment)
        dnldN(study_id, normal_aln_analysis_id)

        // generate Bas for tumour
        basT(dnld.out.files.collect(), file(ref_genome_fa))

        // generate Bas for normal
        basN(dnld.out.files.collect(), file(ref_genome_fa))

        sangerWgs(
            params.mem,
            file(params.sanger_ref_genome_tar),
            file(params.sanger_vagrent_annot),
            file(params.sanger_ref_snv_indel_tar),
            file(params.sanger_ref_cnv_sv_tar),
            file(params.sanger_qcset_tar),
            dnldT.out.files.toSortedList().first().concat(basT.out.bas).collect(),  // aln seq + bas
            dnldT.out.files.toSortedList().last(),  // aln seq idx
            dnldN.out.files.toSortedList().first().concat(basN.out.bas).collect(),  // aln seq + bas
            dnldN.out.files.toSortedList().last(),  // aln seq idx
            params.exclude,
            params.species,
            params.assembly,
            params.skipqc,
            params.pindelcpu,
            params.cavereads,
            params.purity,
            params.ploidy
        )

        // repack results
        repack(sangerWgs.out.result_archive, 'WGS')

        // extract variant calls
        result_patterns = Channel.from(
            'flagged.muts.vcf.gz', 'flagged.vcf.gz', 'copynumber.caveman.vcf.gz', 'annot.vcf.gz')
        extractVar(repack.out.collect(), result_patterns.flatten())

        payloadGenSangerVariant(dnldT.out.song_analysis, extractVar.out)

        // it's not that easy, needs more work here
        upVar(study_id, payloadGenSangerVariant.out.payload, extractVar.out)

        // extract QC
        qc_result_patterns = Channel.from(
            '???', '???')
        extractQC(repack.out.collect(), qc_result_patterns.flatten())

        payloadGenSangerQC(dnldT.out.song_analysis, extractQC.out)

        upQC(study_id, payloadGenSangerQC.out.payload, extractQC.out)


        if (params.cleanup) {
            cleanup(
                dnld.out.files.concat(toLaneBam.out, bwaMemAligner.out, merSorMkdup.out,
                    alignedSeqQC.out, oxog.out, rgQC.out).collect(),
                upAln.out.analysis_id, upQc.out.analysis_id)
        }

    emit:
        analysis_id = upAln.out.analysis_id
        alignment_files = pGenDnaAln.out.alignment_files
}


workflow {
    SangerWgs(
        params.study_id,
        params.analysis_id,
        params.ref_genome_fa
    )
}
