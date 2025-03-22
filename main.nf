#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl=2

// Print pipeline header
log.info """
=================================================
MSI Analysis Nextflow Pipeline
=================================================
Input directory   : ${params.fastq_dir}
Output directory  : ${params.outdir}
Reference genome  : ${params.ref_genome}
Min coverage      : ${params.min_coverage}x
Mapping quality   : ${params.min_mapping_quality}
Threads           : ${params.threads}
=================================================
"""

// Define input channel for paired-end reads
Channel
    .fromFilePairs("${params.fastq_dir}/*_{1,2}.fq.gz", checkIfExists: true)
    .ifEmpty { error "Cannot find any FASTQ pairs in ${params.fastq_dir}" }
    .set { read_pairs_ch }

// Process 1: BWA alignment
process bwa_align {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam")
    
    script:
    """
    bwa mem -R "@RG\\tID:${sample_id}\\tLB:lib_${sample_id}\\tPL:MGISEQ\\tPU:unit1\\tSM:${sample_id}" \
        -t ${params.threads} \
        ${params.ref_genome} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    """
}

// Process 2: Samtools filtering and sorting
process samtools_process {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(sam_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_aligned_rg.bam"), path("${sample_id}_aligned_rg.bam.bai")
    
    script:
    """
    # Filter by mapping quality and sort
    samtools view -@ ${params.threads} -bS -q ${params.min_mapping_quality} ${sam_file} | \
        samtools sort -@ ${params.threads} -o "${sample_id}_aligned_rg.bam"
    
    # Index BAM
    samtools index "${sample_id}_aligned_rg.bam"
    
    # Remove SAM file to save space
    rm ${sam_file}
    """
}

// Process 3: Mark duplicates with GATK
process mark_duplicates {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(bam_file), path(bai_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_aligned_marked.bam"), path("${sample_id}_aligned_marked.bai")
    
    script:
    """
    gatk MarkDuplicates \
        -I ${bam_file} \
        -O "${sample_id}_aligned_marked.bam" \
        -M "${sample_id}_dedup_metrics.txt" \
        --CREATE_INDEX true \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --REMOVE_DUPLICATES false \
        --VALIDATION_STRINGENCY STRICT \
        --ASSUME_SORT_ORDER coordinate
    """
}

// Process 4: Calculate coverage metrics
process calculate_coverage {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(bam_file), path(bai_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_coverage_metrics.txt")
    
    script:
    """
    samtools coverage ${bam_file} > "${sample_id}_coverage_metrics.txt"
    """
}

// Process 5: Prepare MS regions bed file
process prepare_ms_regions {
    output:
    path "ms_regions.bed"
    
    script:
    """
    awk 'BEGIN {OFS="\\t"} NR>1 {
        start=\$2-20  # 20bp padding on each side
        end=\$2+20
        if(start<0) start=0
        print \$1,start,end
    }' ${params.cgp_msisensorpro} > ms_regions.bed
    """
}

// Process 6: Calculate MS coverage and determine threshold
process calculate_ms_coverage {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(bam_file), path(bai_file)
    path ms_regions
    
    output:
    tuple val(sample_id), path("${sample_id}_ms_coverage.txt"), path("${sample_id}_quality_metrics.txt"), path("i_threshold.txt")
    
    script:
    """
    samtools depth -b ${ms_regions} ${bam_file} > "${sample_id}_ms_coverage.txt"
    
    # Calculate MS-specific metrics
    ms_stats=\$(awk '
        BEGIN {sites=0; total_depth=0; covered_sites=0}
        {
            total_depth+=\$3
            sites++
            if(\$3 >= ${params.min_coverage}) covered_sites++
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
    echo "Percent MS sites >${params.min_coverage}x: \$pct_covered%" >> "${sample_id}_quality_metrics.txt"
    echo "Used threshold (-i): \$i_threshold" >> "${sample_id}_quality_metrics.txt"
    
    echo "\$i_threshold" > i_threshold.txt
    """
}

// Process 7: Run MSIsensor-pro
process run_msisensor {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(bam_file), path(bai_file), path(ms_coverage), path(quality_metrics), path(i_threshold_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_MSI"), path("${sample_id}_MSI_dis"), path("${sample_id}_MSI_all"), path("${sample_id}_MSI_unstable")
    
    script:
    """
    i_threshold=\$(cat ${i_threshold_file})
    
    msisensor-pro pro -d ${params.cgp_msisensorpro} \
        -t ${bam_file} \
        -g ${params.ref_genome} \
        -o ${sample_id}_MSI \
        -c ${params.min_coverage} \
        -b ${params.threads} \
        -i \$i_threshold
    """
}

// Process 8: Run MSI analysis
process run_analysis {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(msi), path(msi_dis), path(msi_all), path(msi_unstable)
    
    output:
    tuple val(sample_id), path("${sample_id}_MSI_analysis.json")
    
    script:
    """
    python ${params.lcmsian_script} --sample ${sample_id} \
        --output ${sample_id}_MSI_analysis.json \
        --msi ${msi} \
        --dist ${msi_dis} \
        --all ${msi_all} \
        --unstable ${msi_unstable}
    """
}

// Process 9: Generate reports
process generate_reports {
    tag "${sample_id}"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(analysis_json)
    
    output:
    tuple val(sample_id), path("${sample_id}_MSI_report.pdf"), path("${sample_id}*.html")
    
    script:
    """
    python ${params.lcmsipdf_script} --input ${analysis_json} \
        --output ${sample_id}_MSI_report.pdf
    
    python ${params.lcmsirep_script} ${sample_id}
    """
}

// Define the workflow
workflow {
    // Alignment and BAM processing
    bwa_align(read_pairs_ch)
    samtools_process(bwa_align.out)
    mark_duplicates(samtools_process.out)
    
    // Coverage and MS region analysis
    calculate_coverage(mark_duplicates.out)
    prepare_ms_regions()
    calculate_ms_coverage(mark_duplicates.out, prepare_ms_regions.out)
    
    // MSI analysis and reporting
    ms_coverage_out = calculate_ms_coverage.out
    
    // Join the appropriate channels to create the input for msisensor
    run_msisensor(
        mark_duplicates.out.join(ms_coverage_out, by: 0)
    )
    
    run_analysis(run_msisensor.out)
    generate_reports(run_analysis.out)
}
