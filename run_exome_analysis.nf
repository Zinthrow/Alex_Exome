#!/usr/bin/env nextflow
#! params.base_name = "10847101"
params.base_name = "demo"
params.reads_dir = "/home/alarsen/Projects/Alex_Exome/data" // Path to directory containing paired-end reads
params.ref_index = "GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index"   // Path to reference genome
params.ref_fasta = "GCA_000001405.15_GRCh38_full_analysis_set.fna"
params.output = "./results"  // Output directory

process qualityControl {
    container 'annovar_bioinformatics'

    input:
    path reads_dir
    val base_name

    output:
    val "${base_name}_trimmed_R1.fastq.gz"
    val "${base_name}_trimmed_R2.fastq.gz"

    script:
    """
    mkdir fastqc_reports

    # Quality control using FastQC
    fastqc -o fastqc_reports ${reads_dir}/${base_name}_R1.fastq.gz ${reads_dir}/${base_name}_R2.fastq.gz

    # Trim the reads using Trimmomatic
    java -jar trimmomatic-0.39.jar PE -phred33 ${reads_dir}/${base_name}_R1.fastq.gz ${reads_dir}/${base_name}_R2.fastq.gz ${base_name}_trimmed_R1.fastq.gz ${base_name}_trimmed_unpaired_R1.fastq.gz ${base_name}_trimmed_R2.fastq.gz ${base_name}_trimmed_unpaired_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process alignAndCall {

    input:
    path reads_dir
    val ref_index
    val trimmed_reads1
    val trimmed_reads2
    val base_name

    output:
    val "${base_name}_aligned.bam"

    script:
    """
    # Align paired-end reads using Bowtie2
    bowtie2 -x $read_dir/$ref_index \
    -1 $reads_dir/${trimmed_reads1} \
    -2 $reads_dir${trimmed_reads2} | samtools view -Sb - > $reads_dir/${base_name}_aligned.bam
    """
}

process callVariants {
    container 'broadinstitute_gatk'

    input:
    path reads_dir
    path aligned_bam
    val ref_fasta
    val base_name

    output:
    val "${base_name}_variants.vcf"

    script:
    """
    # Call variants using GATK or other variant caller
    gatk HaplotypeCaller -R $reads_dir/$ref_fasta -I $reads_dir/${base_name}_aligned.bam -O $reads_dir/${base_name}_variants.vcf
    """
}

process annotateVariants {
    container 'annovar_bioinformatics'

    input:
    path variants_vcf
    val base_name

    output:
    val "${base_name}_annotated.vcf"

    script:
    """
    # Placeholder for ANNOVAR or another annotation tool
    # table_annovar.pl ${base_name}_variants.vcf ...  > ${base_name}_annotated.vcf
    """
}

workflow {
    reads = file(params.reads_dir)
    trimmed_reads1, trimmed_reads2, qualityControl(reads, base_name: params.base_name)
    alignAndCall(reads, trimmed_reads1: trimmed_reads1, trimmed_reads2: trimmed_reads2, base_name: params.base_name)
    callVariants(reads, ref_fasta: params.ref_fasta, alignAndCall.out.aligned_bam, base_name: params.base_name)
    annotateVariants(callVariants.out.variants_vcf, base_name: params.base_name)
}
