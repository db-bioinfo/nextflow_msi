#!/usr/bin/env nextflow

/*
========================================================================================
    MSI Analysis Pipeline
========================================================================================
    A Nextflow pipeline for Microsatellite Instability (MSI) analysis.
    
    This pipeline processes sequencing data to analyze microsatellite instability,
    which is an important biomarker in oncology.
========================================================================================
*/

// Define parameters with defaults
params.fastqDir = "$projectDir/data"
params.outdir = "$projectDir/results"
params.refGenome = "/home/administrator/lifecode/genomes/bwa_hg38/hg38.fa"
params.cgpMsisensorpro = "/home/administrator/lifecode/genomes/msi/sensorPRO_hg38/CGP/sensorPRO_microsatellites_targeted"
params.cgp60Msisensorpro = "/home/administrator/lifecode/genomes/msi/sensorPRO_hg38/CGP_60/sensorPRO_microsatellites_targeted_CGP60"

// Clinical thresholds
params.minMappingQuality = 30
params.minCoverage = 100
params.threads = 32

// Print parameters to console
log.info """\
    MSI ANALYSIS PIPELINE
    ===================================
    fastqDir: ${params.fastqDir}
    outdir: ${params.outdir}
    refGenome: ${params.refGenome}
    minMappingQuality: ${params.minMappingQuality}
    minCoverage: ${params.minCoverage}
    threads: ${params.threads}
    """
    .stripIndent()

// Create channels for input files
Channel
    .fromFilePairs("${params.fastqDir}/*_{1,2}.fq.gz", checkIfExists: true)
    .ifEmpty { error "No FASTQ files found in ${params.fastqDir}" }
    .set { read_pairs_ch }

// Include required Python scripts in the pipeline
// These will be staged appropriately by Nextflow
params.python_dir = "$projectDir/bin"
params.lcmsian_script = "${params.python_dir}/lcmsian.py"
params.lcmsipdf_script = "${params.python_dir}/lcmsipdf.py"
params.lcmsirep_script = "${params.python_dir}/lcmsirep.py"

/*
========================================================================================
    PROCESSES
========================================================================================
*/

// Alignment with BWA
process ALIGNMENT {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/bam", mode: 'copy', 
        saveAs: { filename -> filename.endsWith('.bam') || filename.endsWith('.bai') ? filename : null }
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_aligned_rg.bam"), path("${sample_id}_aligned_rg.bam.bai"), emit: aligned_bam
    
    script:
    """
    bwa mem -R "@RG\\tID:${sample_id}\\tLB:lib_${sample_id}\\tPL:MGISEQ\\tPU:unit1\\tSM:${sample_id}" \
        -t ${params.threads} \
        ${params.refGenome} ${reads[0]} ${reads[1]} | \
        samtools view -@ ${params.threads} -bS -q ${params.minMappingQuality} | \
        samtools sort -@ ${params.threads} -o "${sample_id}_aligned_rg.bam"
    
    samtools index "${sample_id}_aligned_rg.bam"
    """
}

// Mark duplicates with GATK
process MARK_DUPLICATES {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/bam", mode: 'copy',
        saveAs: { filename -> 
            if (filename.endsWith('.bam') || filename.endsWith('.bai')) return filename
            else if (filename.endsWith('_metrics.txt')) return "qc/$filename"
            else null
        }
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}_aligned_marked.bam"), path("${sample_id}_aligned_marked.bai"), emit: marked_bam
    path "${sample_id}_dedup_metrics.txt", emit: metrics
    
    script:
    """
    gatk MarkDuplicates \
        -I ${bam} \
        -O "${sample_id}_aligned_marked.bam" \
        -M "${sample_id}_dedup_metrics.txt" \
        --CREATE_INDEX true \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --REMOVE_DUPLICATES false \
        --VALIDATION_STRINGENCY STRICT \
        --ASSUME_SORT_ORDER coordinate
    """
}

// Calculate coverage metrics
process CALCULATE_COVERAGE {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path cgp_msisensorpro
    
    output:
    tuple val(sample_id), path("${sample_id}_coverage_metrics.txt"), emit: coverage_metrics
    tuple val(sample_id), path("${sample_id}_ms_coverage.txt"), emit: ms_coverage
    tuple val(sample_id), path("${sample_id}_quality_metrics.txt"), path("ms_regions.bed"), emit: quality_metrics
    tuple val(sample_id), path(bam), path(bai), env(i_threshold), emit: bam_with_threshold
    
    script:
    """
    samtools coverage ${bam} > "${sample_id}_coverage_metrics.txt"
    
    # Create a bed file for MS regions with padding
    awk 'BEGIN {OFS="\\t"} NR>1 {
        start=\$2-20  # 20bp padding on each side
        end=\$2+20
        if(start<0) start=0
        print \$1,start,end
    }' ${cgp_msisensorpro} > ms_regions.bed
    
    # Get coverage specifically at MS sites
    samtools depth -b ms_regions.bed ${bam} > "${sample_id}_ms_coverage.txt"
    
    # Calculate MS-specific metrics
    ms_stats=\$(awk '
        BEGIN {sites=0; total_depth=0; covered_sites=0}
        {
            total_depth+=\$3
            sites++
            if(\$3 >= ${params.minCoverage}) covered_sites++
        }
        END {
            printf "%.2f\\t%.2f", total_depth/sites, (covered_sites/sites)*100
        }' "${sample_id}_ms_coverage.txt")
    
    # Read the stats
    read mean_depth pct_covered <<< "\$ms_stats"
    
    # Adjust -i threshold based on percentage of MS sites with >100x coverage
    if (( \$(echo "\$pct_covered >= 80" | bc -l) )); then
        i_threshold=0.026  # Standard threshold for high coverage
        echo "High MS site coverage detected (\$pct_covered%): using threshold -i=\$i_threshold"
    else
        i_threshold=0.022  # More stringent threshold for lower coverage
        echo "Lower MS site coverage detected (\$pct_covered%): using conservative threshold -i=\$i_threshold"
    fi
    
    # Add to reporting
    echo "MSI Analysis Quality Metrics" > "${sample_id}_quality_metrics.txt"
    echo "Mean MS site coverage: \$mean_depth" >> "${sample_id}_quality_metrics.txt"
    echo "Percent MS sites >${params.minCoverage}x: \$pct_covered%" >> "${sample_id}_quality_metrics.txt"
    echo "Used threshold (-i): \$i_threshold" >> "${sample_id}_quality_metrics.txt"
    """
}

// Run MSI sensor PRO
process RUN_MSISENSOR {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/msi", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai), val(i_threshold)
    
    output:
    tuple val(sample_id), path("${sample_id}_MSI"), path("${sample_id}_MSI_dis"), 
          path("${sample_id}_MSI_all"), path("${sample_id}_MSI_unstable"), emit: msi_results
    
    script:
    """
    msisensor-pro pro \
        -d ${params.cgpMsisensorpro} \
        -t ${bam} \
        -g ${params.refGenome} \
        -o ${sample_id}_MSI \
        -c ${params.minCoverage} \
        -b ${params.threads} \
        -i ${i_threshold}
    """
}

// Run MSI analysis
process MSI_ANALYSIS {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/results", mode: 'copy'
    
    input:
    tuple val(sample_id), path(msi), path(dis), path(all), path(unstable)
    path lcmsian_script
    
    output:
    tuple val(sample_id), path("${sample_id}_MSI_analysis.json"), emit: analysis_json
    
    script:
    """
    python ${lcmsian_script} \
        --sample ${sample_id} \
        --output ${sample_id}_MSI_analysis.json \
        --msi ${msi} \
        --dist ${dis} \
        --all ${all} \
        --unstable ${unstable}
    """
}

// Generate PDF report
process GENERATE_PDF {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/reports", mode: 'copy'
    
    input:
    tuple val(sample_id), path(analysis_json)
    path lcmsipdf_script
    
    output:
    tuple val(sample_id), path("${sample_id}_MSI_report.pdf"), emit: pdf_report
    
    script:
    """
    python ${lcmsipdf_script} \
        --input ${analysis_json} \
        --output ${sample_id}_MSI_report.pdf
    """
}

// Generate HTML report
process GENERATE_HTML {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/reports", mode: 'copy'
    
    input:
    tuple val(sample_id), path(analysis_json, stageAs: "input_analysis.json")
    path lcmsirep_script
    
    output:
    tuple val(sample_id), path("${sample_id}_MSI_report.html"), emit: html_report
    
    script:
    """
    # Copy the input file to the name expected by the script
    cp input_analysis.json ${sample_id}_MSI_analysis.json
    
    python ${lcmsirep_script} ${sample_id}
    """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {
    // Python scripts
    lcmsian_script_ch = Channel.value(file(params.lcmsian_script))
    lcmsipdf_script_ch = Channel.value(file(params.lcmsipdf_script))
    lcmsirep_script_ch = Channel.value(file(params.lcmsirep_script))
    
    // Execute pipeline processes
    ALIGNMENT(read_pairs_ch)
    MARK_DUPLICATES(ALIGNMENT.out.aligned_bam)
    CALCULATE_COVERAGE(MARK_DUPLICATES.out.marked_bam, file(params.cgpMsisensorpro))
    RUN_MSISENSOR(CALCULATE_COVERAGE.out.bam_with_threshold)
    MSI_ANALYSIS(RUN_MSISENSOR.out.msi_results, lcmsian_script_ch)
    GENERATE_PDF(MSI_ANALYSIS.out.analysis_json, lcmsipdf_script_ch)
    GENERATE_HTML(MSI_ANALYSIS.out.analysis_json, lcmsirep_script_ch)
}

/*
========================================================================================
    COMPLETION SUMMARY
========================================================================================
*/

workflow.onComplete {
    log.info """
    Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Results Dir  : ${params.outdir}
    """
}
