configfile: "config.yaml"


def get_samples(wildcards):
    return config["samples"][wildcards.sample]


rule all:
    input:
        expand("ngmlr_alignment/{sample}.bam.bai", sample=config["samples"]),
        expand("sniffles_calls/{sample}.vcf", sample=config["samples"])

rule ngmlr:
    input:
        fq = get_samples,
        genome = "/home/wdecoster/databases/Homo_sapiens/genome_hg19.fa"
    output:
        protected("ngmlr_alignment/{sample}.bam")
    threads:
        12
    log:
        "logs/ngmlr/{sample}.log"
    shell:
        "zcat {input.fq}/*.fastq.gz | \
         ngmlr -x ont -t {threads} -r {input.genome} | \
         samtools sort -@ {threads} -o {output} - 2> {log}"


rule samtools_index:
    input:
        "ngmlr_alignment/{sample}.bam"
    output:
        "ngmlr_alignment/{sample}.bam.bai"
    threads: 12
    log:
        "logs/samtools_index/{sample}.log"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"


rule sniffles:
    input:
        "ngmlr_alignment/{sample}.bam"
    output:
        protected("sniffles_calls/{sample}.vcf")
    threads: 12
    log:
        "logs/sniffles/{sample}.log"
    shell:
        "sniffles --mapped_reads {input} --vcf {output} --threads {threads} 2> {log}"
