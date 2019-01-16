def get_samples(wildcards):
    return config["samples"][wildcards.sample]


rule minimap2_align:
    input:
        fq = get_samples,
        genome = config["genome"]
    output:
        "minimap2/alignment/{sample}.bam"
    threads:
        8
    log:
        "logs/minimap2/{sample}.log"
    shell:
        "minimap2 --MD -ax map-ont -t {threads} {input.genome} {input.fq}/*.fastq.gz | \
         samtools sort -@ {threads} -o {output} - 2> {log}"

rule ngmlr_align:
    input:
        fq = get_samples,
        genome = config["genome"]
    output:
        protected("ngmlr/alignment/{sample}.bam")
    threads:
        36
    log:
        "logs/ngmlr/{sample}.log"
    shell:
        "zcat {input.fq}/*.fastq.gz | \
         ngmlr --presets ont -t {threads} -r {input.genome} | \
         samtools sort -@ {threads} -o {output} - 2> {log}"

rule minimap2_last_like_align:
    input:
        fq = get_samples,
        genome = config["genome"]
    output:
        "minimap2_last_like/alignment/{sample}.bam"
    threads:
        8
    log:
        "logs/minimap2_last_like/{sample}.log"
    shell:
        "minimap2 --MD -ax map-ont --no-long-join -r50 -t {threads} {input.genome} {input.fq}/*.fastq.gz | \
         samtools sort -@ {threads} -o {output} - 2> {log}"


rule samtools_index:
    input:
        "{aligner}/alignment/{sample}.bam"
    output:
        "{aligner}/alignment/{sample}.bam.bai"
    threads: 4
    log:
        "logs/{aligner}/samtools_index/{sample}.log"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"

rule alignment_stats:
    input:
        bam = expand("{{aligner}}/alignment/{sample}.bam", sample=config["samples"]),
        bai = expand("{{aligner}}/alignment/{sample}.bam.bai", sample=config["samples"])
    output:
        "{aligner}/alignment_stats/{sample}.txt"
    log:
        "logs/{aligner}/alignment_stats/{sample}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/alignment_stats.py") + \
            " -o {output} {input.bam} 2> {log}"
