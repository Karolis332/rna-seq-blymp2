import os

# Samples
samples = [f.replace(".fastq.gz", "") for f in os.listdir("data") if f.endswith(".fastq.gz")]

# Paths
reference_dir = "references/hg38"
reference_fa = os.path.join(reference_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
reference_gtf = os.path.join(reference_dir, "Homo_sapiens.GRCh38.110.gtf")
hisat2_index_base = os.path.join(reference_dir, "hg38_index")

rule all:
    input:
        expand("results/fastqc/{sample}_fastqc.html", sample=samples),
        expand("results/trimmed/{sample}_trimmed.fastq.gz", sample=samples),
        expand("results/bam/{sample}.bam", sample=samples),
        expand("results/counts/{sample}.counts.txt", sample=samples),
        "results/multiqc/multiqc_report.html"

rule fastqc:
    input:
        "data/{sample}.fastq.gz"
    output:
        html = "results/fastqc/{sample}_fastqc.html",
        zip = "results/fastqc/{sample}_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    shell:
        "mkdir -p results/fastqc && fastqc {input} --outdir results/fastqc"

rule cutadapt:
    input:
        "data/{sample}.fastq.gz"
    output:
        "results/trimmed/{sample}_trimmed.fastq.gz"
    conda:
        "envs/cutadapt.yaml"
    shell:
        "mkdir -p results/trimmed && cutadapt -q 20 -o {output} {input}"

rule hisat2_align:
    input:
        fq = "results/trimmed/{sample}_trimmed.fastq.gz",
        index = hisat2_index_base + ".1.ht2"
    output:
        temp("results/bam/{sample}.sam")
    conda:
        "envs/hisat2.yaml"
    threads: 4
    shell:
        """
        mkdir -p results/bam
        hisat2 -x {hisat2_index_base} -U {input.fq} -S {output} --threads {threads}
        """

rule sam_to_bam:
    input:
        "results/bam/{sample}.sam"
    output:
        "results/bam/{sample}.bam"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools sort -o {output} {input} && samtools index {output}"

rule featurecounts:
    input:
        bam = "results/bam/{sample}.bam",
        gtf = reference_gtf
    output:
        "results/counts/{sample}.counts.txt"
    conda:
        "envs/subread.yaml"
    threads: 2
    shell:
        """
        mkdir -p results/counts
        featureCounts -a {input.gtf} -o {output} -T {threads} {input.bam}
        """

rule multiqc:
    input:
        expand("results/fastqc/{sample}_fastqc.zip", sample=samples)
    output:
        "results/multiqc/multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        "mkdir -p results/multiqc && multiqc results -o results/multiqc"
