
process MANTLE_STAGE_INPUTS {
    tag "${pipeline_run_id}_mantleSDK_stageInputs"

    secret 'MANTLE_USER'
    secret 'MANTLE_PASSWORD'

    container 'public.ecr.aws/c7j2m0e6/mantle-sdk:latest'

    publishDir "${outdir}/stage_inputs", mode: 'copy'

    input:
    val outdir
    val pipeline_run_id

    output:
    path('*.fastq.gz'), emit: fastq_ch

    script:
    def stage_directory = "./"

    """
    get_data.py ${pipeline_run_id} ${stage_directory} \
                --tenant ${TENANT} \
                --mantle_env ${ENV}
    """
}

process MERGE {
    tag "${pipeline_run_id}_AssemblyPipeline"

    publishDir "${params.outdir}/phage-assembly", mode: 'copy'

    container '663344187369.dkr.ecr.eu-central-1.amazonaws.com/phage_assembly_mantle:latest'

    input:
    path outdir
    val fastqfile

    output:
    path('*merged_to_assemble.fastq.gz'), emit: merged_fq

    script:

    """
    merge_fastqs.sh -i ${fastqfolder} -o ${outdir}
    """
}

process MANTLE_UPLOAD_RESULTS {
    tag "${pipeline_run_id}-mantleSDK_uploadResults"

    publishDir "${params.outdir}/mantle_upload_results", mode: 'copy'

    secret 'MANTLE_USER'
    secret 'MANTLE_PASSWORD'

    container 'public.ecr.aws/c7j2m0e6/mantle-sdk:latest'

    input:
    val pipeline_run_id
    path outdir, stageAs: 'results/*'
    val merged_fq

    output:
    tuple val(pipeline_run_id), path('*.txt'), emit: completion_timestamp

    script:

    """
    upload_data.py ${pipeline_run_id} ${outdir} \
                --tenant ${TENANT} \
                --mantle_env ${ENV}

    date > results_uploaded_mantle.txt
    """
}

workflow {
     // Get FatsQs and sample metadata using pipeline Run ID from mantle SDK
    MANTLE_STAGE_INPUTS (
        params.outdir,
        params.pipeline_run_id
    )

    MERGE (
        params.outdir,
        MANTLE_STAGE_INPUTS.out.fastq_ch
    )
    // Sync outputs back into mantle
    MANTLE_UPLOAD_RESULTS (
        params.pipeline_run_id,
        params.outdir,
        MERGE.out.merged.fq
    )
}
