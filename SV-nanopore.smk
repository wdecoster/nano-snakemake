configfile: "config.yaml"


def get_samples(wildcards):
    return config["samples"][wildcards.sample]


rule all:
    input:
        expand("ngmlr_alignment/{sample}.bam.bai", sample=config["samples"]),
        expand("sniffles_genotypes/{sample}.vcf", sample=config["samples"]),
        expand("mosdepth/{sample}.mosdepth.dist.txt", sample=config["samples"]),


rule ngmlr:
    input:
        fq = get_samples,
        genome = "/home/wdecoster/databases/Homo_sapiens/GRCh38_recommended/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    output:
        protected("ngmlr_alignment/{sample}.bam")
    threads:
        24
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


rule sniffles_call:
    input:
        "ngmlr_alignment/{sample}.bam"
    output:
        protected("sniffles_calls/{sample}.vcf")
    threads: 24
    log:
        "logs/sniffles_call/{sample}.log"
    shell:
        "sniffles --mapped_reads {input} --vcf {output} --threads {threads} 2> {log}"


rule sniffles_genotype:
    input:
        bam = "ngmlr_alignment/{sample}.bam",
        ivcf = "sniffles_combined/calls.vcf"
    output:
        protected("sniffles_genotypes/{sample}.vcf")
    threads: 24
    log:
        "logs/sniffles_genotype/{sample}.log"
    shell:
        "sniffles --mapped_reads {input.bam} --vcf {output} --threads {threads} --Ivcf {input.ivcf} 2> {log}"


rule survivor_calls:
    input:
        expand("sniffles_calls/{sample}.vcf", sample=config["samples"])
    output:
        protected("sniffles_combined/calls.vcf")
    params:
        distance = 1000,
        caller_support = 1,
        same_type = 1,
        same_strand = 0,
        estimate_distance = 0,
        minimum_size = 0,
    log:
        "logs/sniffles/combine_calls.log"
    shell:
        "ls {input} > sniffles_calls/samples.fofn ; \
        SURVIVOR merge sniffles_calls/samples.fofn {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output} 2> {log}"


rule survivor_genotypes:
    input:
        expand("sniffles_genotypes/{sample}.vcf", sample=config["samples"])
    output:
        protected("sniffles_combined/genotypes.vcf")
    params:
        distance = 1000,
        caller_support = 0,
        same_type = 1,
        same_strand = 0,
        estimate_distance = 0,
        minimum_size = 0,
    log:
        "logs/sniffles/combine_genotypes.log"
    shell:
        "ls {input} > sniffles_calls/samples.fofn ; \
        SURVIVOR merge sniffles_calls/samples.fofn {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output} 2> {log}"


rule mosdepth:
    input:
        bam = "ngmlr_alignment/{sample}.bam"
    threads: 4
    output:
        protected("mosdepth/{sample}.mosdepth.dist.txt")
    params:
        windowsize = 1000,
        prefix = "{sample}",
    log:
        "logs/mosdepth/mosdepth_{sample}.log"
    shell:
        "mosdepth --threads {threads} --by {params.windowsize} mosdepth/{params.prefix} {input} 2> {log}"


# annotate vcf
# add mosdepth information and plots on called sites
