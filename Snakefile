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
        "results/multiqc/multiqc_report.html",
        "results/counts/merged_counts.tsv"

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

rule fastp_trim:
    input:
        "data/{sample}.fastq.gz"
    output:
        trimmed = "results/trimmed/{sample}_trimmed.fastq.gz",
        json = "results/trimmed/{sample}_fastp.json",
        html = "results/trimmed/{sample}_fastp.html"
    conda:
        "envs/fastp.yaml"
    shell:
        """
        mkdir -p results/trimmed
        fastp -i {input} -o {output.trimmed} --detect_adapter_for_pe -h {output.html} -j {output.json}
        """

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

rule merge_counts:
    input:
        counts = expand("results/counts/{sample}.counts.txt", sample=samples),
        logs = expand("results/counts/{sample}.counts.txt.summary", sample=samples),
        gtf = reference_gtf
    output:
        matrix = "results/counts/merged_counts.tsv",
        gene_list = "results/counts/merged_gene_names.txt",
        summary = "results/counts/merged_counts_summary.tsv",
        gene_mapped = "results/counts/merged_counts_with_genes.tsv",
        tpm = "results/counts/tpm_normalized.tsv",
        cpm = "results/counts/cpm_normalized.tsv"
    conda:
        "envs/python.yaml"
    run:
        import pandas as pd
        import os
        import re

        # --- Merge raw counts ---
        dfs = []
        for f in input.counts:
            sample_name = os.path.basename(f).replace(".counts.txt", "")
            df = pd.read_csv(f, sep="\t", comment='#', index_col=0)
            df = df.iloc[:, [-1]]
            df.columns = [sample_name]
            dfs.append(df)

        merged = pd.concat(dfs, axis=1)
        merged.to_csv(output.matrix, sep="\t")

        with open(output.gene_list, "w") as out_f:
            for gene in merged.index:
                out_f.write(gene + "\n")

        # --- Merge summary files ---
        summary_frames = []
        for log_file in input.logs:
            sample = os.path.basename(log_file).replace(".counts.txt.summary", "")
            df = pd.read_csv(log_file, sep="\t", index_col=0, header=None)
            df.columns = [sample]
            summary_frames.append(df)

        summary_merged = pd.concat(summary_frames, axis=1)
        summary_merged.to_csv(output.summary, sep="\t")

        # --- Map Ensembl IDs to gene names using GTF ---
        gene_map = {}
        with open(input.gtf, "r") as gtf_file:
            for line in gtf_file:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if fields[2] != "gene":
                    continue
                info = fields[8]
                gene_id = re.search('gene_id "(.*?)"', info)
                gene_name = re.search('gene_name "(.*?)"', info)
                if gene_id and gene_name:
                    gene_map[gene_id.group(1)] = gene_name.group(1)

        renamed = merged.copy()
        renamed.index = [gene_map.get(gid, gid) for gid in renamed.index]
        renamed.to_csv(output.gene_mapped, sep="\t")

        # --- TPM normalization ---
        def tpm(df):
            rpk = df.div(1)  # assume length = 1 for simplicity; modify if gene lengths are known
            scaling_factor = rpk.sum(axis=0) / 1e6
            return rpk.div(scaling_factor, axis=1)

        tpm_df = tpm(merged)
        tpm_df.to_csv(output.tpm, sep="\t")

        # --- CPM normalization ---
        def cpm(df):
            scaling_factor = df.sum(axis=0) / 1e6
            return df.div(scaling_factor, axis=1)

        cpm_df = cpm(merged)
        cpm_df.to_csv(output.cpm, sep="\t")

