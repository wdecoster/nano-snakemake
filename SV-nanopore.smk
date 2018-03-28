configfile: "config.yaml"


def get_samples(wildcards):
    return config["samples"][wildcards.sample]


rule all:
    input:
        expand("SV-plots/SV-length_genotypes_{sample}.png", sample=config["samples"]),
        expand("SV-plots/SV-length_calls_{sample}.png", sample=config["samples"]),
        expand("mosdepth/{sample}.regions.bed.gz", sample=config["samples"]),
        expand("mosdepth/{sample}.mosdepth.dist.txt", sample=config["samples"]),
        "sniffles_combined/annot_genotypes.vcf",
        "mosdepth/regions.combined.gz",


rule minimap2:
    input:
        fq = get_samples,
        genome = "/home/wdecoster/databases/Homo_sapiens/GRCh38_recommended/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    output:
        "minimap2_alignment/{sample}.bam"
    threads:
        8
    log:
        "logs/minimap2/{sample}.log"
    shell:
        "minimap2 -a -t {threads} {input.genome} {input.fq} | samtools sort -@ {threads} -o {output} - 2> {log}"

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
        "sniffles --mapped_reads {input.bam} \
                  --vcf {output} \
                  --threads {threads} \
                  --Ivcf {input.ivcf} 2> {log}"


rule survivor_calls:
    input:
        expand("sniffles_calls/{sample}.vcf", sample=config["samples"])
    output:
        protected("sniffles_combined/calls.vcf")
    params:
        distance = 1000,
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
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
        temp("sniffles_combined/genotypes.vcf")
    params:
        distance = 1000,
        caller_support = 0,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        "logs/sniffles/combine_genotypes.log"
    shell:
        "ls {input} > sniffles_calls/samples.fofn ; \
        SURVIVOR merge sniffles_calls/samples.fofn {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output} 2> {log}"


rule nanosv:
    input:
        bam = "ngmlr_alignment/{sample}.bam",
        bai = "ngmlr_alignment/{sample}.bam.bai"
    output:
        "nanosv/{sample}.vcf"
    params:
        bed = "/home/wdecoster/databases/Homo_sapiens/GRCh38_recommended/GRCh38_full_annotation.bed.gz",
        samtools = "samtools"
    log:
        "logs/nanosv/{sample}.log"
    shell:
        "NanoSV --bed {params.bed} -s {params.samtools} {input.bam} -o {output} 2> {log}"


rule survivor_nanosv:
    input:
        expand("nanosv/{sample}.vcf", sample=config["samples"])
    output:
        temp("nanosv_combined/nanosv.vcf")
    params:
        distance = 1000,
        caller_support = 0,
        same_type = 1,
        same_strand = 0,
        estimate_distance = 0,
        minimum_size = 0,
    log:
        "logs/survivor/combine_nanosv.log"
    shell:
        "ls {input} > sniffles_calls/samples.fofn ; \
        SURVIVOR merge sniffles_calls/samples.fofn {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output} 2> {log}"


rule mosdepth:
    input:
        bam = "ngmlr_alignment/{sample}.bam",
        bai = "ngmlr_alignment/{sample}.bam.bai"
    threads: 4
    output:  # change if mosdepth 0.2.2
        protected("mosdepth/{sample}.mosdepth.dist.txt"),
        protected("mosdepth/{sample}.regions.bed.gz"),
    params:
        windowsize = 500,
        prefix = "{sample}",
    log:
        "logs/mosdepth/mosdepth_{sample}.log"
    shell:
        "mosdepth --threads {threads} \
                  -n \
                  --by {params.windowsize} \
                  mosdepth/{params.prefix} {input.bam} 2> {log}"


rule mosdepth_combine:
    input:
        expand("mosdepth/{sample}.regions.bed.gz", sample=config["samples"])
    output:
        "mosdepth/regions.combined.gz"
    log:
        "logs/mosdepth/mosdepth_combine.log"
    shell:
        "python ~/projects/SV-snakemake/scripts/combine_mosdepth.py {input} -o {output} 2> {log}"


rule mosdepth_global_plot:
    input:   # change if mosdepth 0.2.2
        expand("mosdepth/{sample}.mosdepth.dist.txt", sample=config["samples"])
    output:
        protected("mosdepth_global_plot/global.html")
    log:
        "logs/mosdepth/mosdepth_global_plot.log"
    shell:
        "python scripts/plot_dist.py {input} -o {output} 2> {log}"


rule SV_length_plot_genotypes:
    input:
        "sniffles_genotypes/{sample}.vcf"
    output:
        "SV-plots/SV-length_genotypes_{sample}.png"
    log:
        "logs/svplot/svlength_{sample}.log"
    shell:
        "python ~/projects/SV-snakemake/scripts/SV-length-plot.py {input} {output} 2> {log}"


rule SV_length_plot_calls:
    input:
        "sniffles_calls/{sample}.vcf"
    output:
        "SV-plots/SV-length_calls_{sample}.png"
    log:
        "logs/svplot/svlength_{sample}.log"
    shell:
        "python ~/projects/SV-snakemake/scripts/SV-length-plot.py {input} {output} 2> {log}"


rule SV_plot_carriers:
    input:
        "sniffles_combined/genotypes.vcf"
    output:
        "SV-plots/SV-carriers.png"
    log:
        "logs/svplot/svcarriers.log"
    shell:
        "python ~/projects/SV-snakemake/scripts/SV-carriers-plot.py {input} {output} 2> {log}"


rule sort_vcf:
    input:
        "sniffles_combined/genotypes.vcf"
    output:
        temp("sniffles_combined/sorted_genotypes.vcf")
    log:
        "logs/sort_vcf/sorting_combined_genotypes.log"
    threads: 8
    shell:
        "vcf-sort --parallel {threads} {input} > {output} 2> {log}"


rule annotate_vcf:
    input:
        "sniffles_combined/sorted_genotypes.vcf"
    output:
        protected("sniffles_combined/annot_genotypes.vcf")
    log:
        "logs/annotate_vcf/annotate_genotypes.log"
    params:
        conf = "/home/wdecoster/projects/SV-snakemake/configuration/vcfanno_conf.toml"
    threads: 8
    shell:
        "vcfanno -ends -p {threads} {params.conf} {input} > {output} 2> {log}"


# add mosdepth information and plots on called sites
