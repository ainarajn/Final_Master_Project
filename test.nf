#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//Define parameters
params.genome = "hg19" //Options: hg19, hg19tohg38, hg38, hg38tohg19
params.analysis = "JPALM" //Options: PALM, JPALM
params.GWAS_list = "/homes/users/ajimenez/scratch/prueba/lista_gwas.txt"
params.GWASdir = "/homes/users/ajimenez/scratch/prueba"
params.ref_effect = true
params.population = "GCAT"
params.Relate_SNPs = "$projectDir/gcat/relate_snps_${params.population}.rds"

params.VCFdir = "/gpfs/projects/lab_dcomas/1000genomes_phase3_dcomas/vcf"
params.Hapmap3 = (params.genome in ["hg19", "hg19tohg38"]) ? "$projectDir/gcat/ldscore/hm3_SNPs_variant_id.tsv" : 
                 (params.genome in ["hg38", "hg38tohg19"]) ? "$projectDir/gcat/ldscore/hm3_SNPs_variant_id_hg38.tsv" : null
params.LDblocks = (params.genome in ["hg19", "hg19tohg38"]) ? "$projectDir/gcat/ld_blocks.bed" : 
                  (params.genome in ["hg38", "hg38tohg19"]) ? "$projectDir/gcat/ld_blocks_hg38.bed" : null

params.maxp = 5e-8
params.N_per_batch = 10
params.min_SI_SNPs = 25

log.info """\
    G W A S   ${params.analysis}  ANALYSIS -  N F   P I P E L I N E
    =================================================================
    Genome	    : ${params.genome}
    Population	: ${params.population}
    GWAS_list	: ${params.GWAS_list}
    Relate      : ${params.Relate_SNPs}

    Max. p-val	: ${params.maxp}
    """
    .stripIndent()

process Unlist_GWAS {

    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 8

    input:
    path GWAS_list

    output:
    path 'GWAS_*'

    script:
    """
    unlist_GWAS.sh ${GWAS_list}
    """

}

process GWAS_format {
    module = 'R/4.2.0-foss-2021b'
    clusterOptions = '--partition=haswell'

    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 8

    input:
    path GWAS
    path GWASdir
    path VCFdir
    val ref_effect
    path Hapmap3

    output:
    path 'format_*.txt'

    script:
    """
    gwas_format.R ${GWAS} ${GWASdir} ${VCFdir} ${ref_effect} ${Hapmap3} 
    """

}

process LD_blocks {
    module = 'R/4.2.0-foss-2021b'
    clusterOptions = '--partition=haswell'

    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 3

    //publishDir InterDir, mode: 'copy'

    input:
    path GWAS
    path GWASdir
    path LD_Blocks
    val MaxPval

    output:
    path '*-+*'

    script:
    """
    LDblock_filter.R ${GWAS} ${GWASdir} ${LD_Blocks} ${MaxPval}
    """

}

process LiftOver {

    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 5

    input:
    path GWAS
    path GWAS_dir
    val genome

    output:
    path '*_*.txt'

    when: 
    (params.genome in ["hg19tohg38", "hg38tohg19"]) &&
    GWAS.exists() &&
        GWAS.name.tokenize('-+').size() == 2 &&
        Float.parseFloat(GWAS.name.tokenize('-+')[1]) >= params.min_SI_SNPs

    script:
    """
    gwas_liftover.sh ${GWAS} ${GWAS_dir} ${genome}
    """
}

process Obtain_N_Comparisons {

    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 1

    input:
    path GWAS_list

    output:
    path 'n_pairwise_comparisons.txt'

    when:
    params.analysis == "JPALM"
 
    script:
    """
    n_comparisons.sh ${GWAS_list}
    """

}

process Obtain_GWAS_Pairs {

    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 1

    input:
    path GWAS_list

    output:
    path 'comparison_batch_*'

    when:
    params.analysis == "JPALM"
 
    script:
    """
    gwas_pairs_batch_obtention.sh ${GWAS_list}
    """

}

process Munge_Sumstats {

    clusterOptions = '--partition=haswell'

    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path GWAS
    path GWASdir

    output:
    path '*.sumstats.info'

    when:
    params.analysis == "JPALM" &&
    GWAS.exists() &&
        GWAS.name.tokenize('-+').size() == 2 &&
        Float.parseFloat(GWAS.name.tokenize('-+')[1]) >= params.min_SI_SNPs

    script:
    """
    munge_sumstats.sh ${GWAS} ${GWASdir}
    """

}

process Genetic_Correlation {

    clusterOptions = '--partition=haswell'

    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    val ready
    path GWAS_Pair
    path N_Comparisons

    output:
    path '*-+*-+*-+*'

    when:
    params.analysis == "JPALM"

    script:
    """
    genetic_correlation_batch.sh ${GWAS_Pair} ${N_Comparisons}
    """

}

process SelectSNPs {

    module = 'R/4.2.0-foss-2021b'
    clusterOptions = '--partition=haswell'

    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    path relate_SNPs
    path GWAS
    path GWAS_dir

    output:
    path '*_selected_SNPs.tsv'

    when:
    params.analysis == "PALM" &&
    (
        (params.genome in ["hg19", "hg38"] && 
        GWAS.exists() && 
        GWAS.name.tokenize('-+').size() == 2 && 
        Float.parseFloat(GWAS.name.tokenize('-+')[1]) >= params.min_SI_SNPs) 
        || 
        params.genome in ["hg19tohg38", "hg38tohg39"]
    )

    script:
    """
    snp_selection.R ${relate_SNPs} ${GWAS} ${GWAS_dir}
    """

}

process SelectSNPs_JPALM {

    module = 'R/4.2.0-foss-2021b'
    clusterOptions = '--partition=haswell'

    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path rg
    path relate_SNPs
    path GWAS_dir
    val MaxPval

    output:
    path '*-*_selected_SNPs.tsv'

    when:
    
    params.analysis == "JPALM" &&
    params.genome in ["hg19", "hg38"] && 
    rg.exists()  &&
        rg.name.tokenize('-+').size() == 4 &&
        Math.abs(Float.parseFloat(rg.name.tokenize('-+')[2])) > 0.2 &&
        Float.parseFloat(rg.name.tokenize('-+')[3]) < 0.005

    script:
    """
    snp_selection_jpalm.R ${rg} ${relate_SNPs} ${GWAS_dir} ${MaxPval}
    """

}

process SelectSNPs_JPALM_lifted {

    module = 'R/4.2.0-foss-2021b'
    clusterOptions = '--partition=haswell'

    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    val ready
    path rg
    path relate_SNPs
    path GWAS_dir
    val MaxPval
    val genome

    output:
    path '*-*_selected_SNPs.tsv'

    when:
    params.analysis == "JPALM" &&
    params.genome in ["hg19tohg38", "hg38tohg19"] && 
    rg.exists()  &&
        rg.name.tokenize('-+').size() == 4 &&
        Math.abs(Float.parseFloat(rg.name.tokenize('-+')[2])) > 0.2 &&
        Float.parseFloat(rg.name.tokenize('-+')[3]) < 0.005

    script:
    """
    snp_selection_jpalm_hg38.R ${rg} ${relate_SNPs} ${GWAS_dir} ${MaxPval} ${genome}
    """
}

process MergeSNPs {

    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    path Selected_SNPs_per_GWAS

    output:
    path 'merged_SNPs.tsv'

    script:
    """
    merge_SNPs.sh ${Selected_SNPs_per_GWAS}
    """

}

process TrimSNPs {

    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 1

    input:
    path Selected_SNPs
    val MaxPval
    val analysis
    val N_per_batch

    output:
    path 'batch_*_*.tsv'

    script:
    """
    trim_snps_batch.sh ${Selected_SNPs} ${MaxPval} ${analysis} ${N_per_batch}
    """

}

process SBL_and_SNP_Likelihood {

    clusterOptions = '--partition=haswell'

    time '1h'
    memory 8.GB

    input:
    path SNP_batch
    val population
    val genome

    output:
    path 'bp*.quad_fit.npy'

    script:
    """
    relate_and_snp_likelihood.sh ${SNP_batch} ${population} ${genome}
    """

}

process Apply_PALM {

    clusterOptions = '--partition=haswell'

    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    val ready
    path Filtered_GWAS
    val MaxPvalue
    val genome
    val Population

    output:
    path '*_marginal_PALM.txt'

    when:
    params.analysis == "PALM"

    script:
    """
    palm.sh ${Filtered_GWAS} ${MaxPvalue} ${genome} ${Population}
    """

}

process Apply_JPALM {

    clusterOptions = '--partition=haswell'

    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    val ready
    path Filtered_GWAS_Pair
    val MaxPvalue
    val genome
    val Population

    output:
    path '*_J_PALM.txt'
    path '*_significant_independent_SNPs.txt'

    when:
    params.analysis == "JPALM"

    script:
    """
    j_palm.sh ${Filtered_GWAS_Pair} ${MaxPvalue} ${genome} ${Population}
    """

}

workflow {
    Unlist_GWAS_ch = Unlist_GWAS(params.GWAS_list)
    GWAS_format_ch = GWAS_format(Unlist_GWAS_ch.flatten(), params.GWASdir, params.VCFdir, params.ref_effect, params.Hapmap3)
    LDblocks_ch = LD_blocks(GWAS_format_ch, params.GWASdir, params.LDblocks, params.maxp)
    Lifted_GWAS_ch = LiftOver(LDblocks_ch, params.GWASdir, params.genome)
    GWAS_N_ch = Obtain_N_Comparisons(LDblocks_ch.collect())
    GWAS_Pairs_ch = Obtain_GWAS_Pairs(LDblocks_ch.collect())
    Sumstats_ch = Munge_Sumstats(LDblocks_ch, params.GWASdir)
    RG_ch = Genetic_Correlation(Sumstats_ch.collect(), GWAS_Pairs_ch.flatten(), GWAS_N_ch)
    if (params.genome in ['hg19tohg38', 'hg38tohg19']) {
        if (params.analysis == "JPALM") {
            SNPs_ch = SelectSNPs_JPALM_lifted(Lifted_GWAS_ch.collect(), RG_ch.flatten(), params.Relate_SNPs, "${params.GWASdir}_${params.genome}", params.maxp, params.genome)
        } else {
            SNPs_ch = SelectSNPs(params.Relate_SNPs, Lifted_GWAS_ch, "${params.GWASdir}_${params.genome}")
        }  
    } else {
         if (params.analysis == "JPALM") {
            SNPs_ch = SelectSNPs_JPALM(RG_ch.flatten(), params.Relate_SNPs, params.GWASdir, params.maxp)
        } else {
            SNPs_ch = SelectSNPs(params.Relate_SNPs, LDblocks_ch, params.GWASdir)
        }     
    }
    Merged_SNPs_ch = MergeSNPs(SNPs_ch.collect())
    Trimmed_SNPs_ch = TrimSNPs(Merged_SNPs_ch, params.maxp, params.analysis, params.N_per_batch)
    SBL_Lik_ch = SBL_and_SNP_Likelihood(Trimmed_SNPs_ch.flatten(), params.population, params.genome)
    PALM_ch = Apply_PALM(SBL_Lik_ch.collect(), SNPs_ch, params.maxp, params.genome, params.population)
    J_PALM_ch = Apply_JPALM(SBL_Lik_ch.collect(), SNPs_ch, params.maxp, params.genome, params.population)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops... something went wrong" )
}
