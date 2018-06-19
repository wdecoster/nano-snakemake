import os
import gzip

configfile: "config.yaml"


def get_samples(wildcards):
    return config["samples"][wildcards.sample]


def get_chromosomes(genome, annotation):
    '''
    Gets the chromosome identifiers from the fasta genome and bed annotation
    and returns the intersection of both
    '''
    fai = genome + ".fai"
    if not os.path.isfile(fai):
        sys.exit("Fasta index {} not found".format(fai))
    fa_chr = set([i.split('\t')[0] for i in open(fai)])
    if not os.path.isfile(annotation):
        sys.exit("Annotation file {} not found".format(annotation))
    annot_chr = set([line.split('\t')[0] for line in gzip.open(annotation, 'rt')])
    return list(fa_chr & annot_chr)


CHROMOSOMES = get_chromosomes(config["genome"], config["annotbed"])


rule all:
    input:
        expand("SV-plots/SV-length_{caller}_genotypes_{sample}.png",
               sample=config["samples"],
               caller=["sniffles", "nanosv"]),
        expand("SV-plots/SV-{caller}_carriers.png", caller=["sniffles", "nanosv"]),
        expand("{caller}_combined/annot_genotypes.vcf", caller=["sniffles", "nanosv"]),
        "alignment_stats/ngmlr.txt",
        "all_combined/annot_genotypes.vcf",
        "high_confidence_combined/annot_genotypes.vcf",
        "mosdepth/regions.combined.gz",
        "mosdepth_global_plot/global.html",


rule sniffles:
    input:
        "sniffles_combined/annot_genotypes.vcf",
        expand("SV-plots/SV-length_sniffles_genotypes_{sample}.png", sample=config["samples"]),
        "SV-plots/SV-sniffles_carriers.png"

rule nanosv:
    input:
        "nanosv_combined/annot_genotypes.vcf",
        expand("SV-plots/SV-length_nanosv_genotypes_{sample}.png", sample=config["samples"]),
        "SV-plots/SV-nanosv_carriers.png"

rule mosdepth:
    input:
        "mosdepth/regions.combined.gz",
        "mosdepth_global_plot/global.html",

rule minimap2:
    input:
        expand("minimap2_alignment/{sample}.bam.bai", sample=config["samples"]),
        "alignment_stats/minimap2.txt",


#######################################################################################

rule minimap2_align:
    input:
        fq = get_samples,
        genome = config["genome"]
    output:
        "minimap2_alignment/{sample}.bam"
    threads:
        8
    log:
        "logs/minimap2/{sample}.log"
    shell:
        "minimap2 --MD -a -t {threads} {input.genome} {input.fq}/*.fastq.gz | \
         samtools sort -@ {threads} -o {output} - 2> {log}"

rule ngmlr:
    input:
        fq = get_samples,
        genome = config["genome"]
    output:
        protected("ngmlr_alignment/{sample}.bam")
    threads:
        36
    log:
        "logs/ngmlr/{sample}.log"
    shell:
        "zcat {input.fq}/*.fastq.gz | \
         ngmlr --presets ont -t {threads} -r {input.genome} | \
         samtools sort -@ {threads} -o {output} - 2> {log}"


rule samtools_index:
    input:
        "{aligner}_alignment/{sample}.bam"
    output:
        "{aligner}_alignment/{sample}.bam.bai"
    threads: 4
    log:
        "logs/samtools_index/{aligner}_{sample}.log"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"


rule alignment_stats:
    input:
        bam = expand("{{aligner}}_alignment/{sample}.bam", sample=config["samples"]),
        bai = expand("{{aligner}}_alignment/{sample}.bam.bai", sample=config["samples"])
    output:
        "alignment_stats/{aligner}.txt"
    log:
        "logs/alignment_stats/{aligner}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/alignment_stats.py") + \
            " -o {output} {input.bam} 2> {log}"

rule sniffles_call:
    input:
        "ngmlr_alignment/{sample}.bam"
    output:
        "sniffles_calls/{sample}.vcf"
    threads: 8
    log:
        "logs/sniffles_call/{sample}.log"
    shell:
        "sniffles --mapped_reads {input} --vcf {output} --threads {threads} 2> {log}"


rule sniffles_genotype:
    input:
        bam = "ngmlr_alignment/{sample}.bam",
        ivcf = "sniffles_combined/calls.vcf"
    output:
        "sniffles_genotypes/{sample}.vcf"
    threads: 8
    log:
        "logs/sniffles_genotype/{sample}.log"
    shell:
        "sniffles --mapped_reads {input.bam} \
                  --vcf {output} \
                  --threads {threads} \
                  --Ivcf {input.ivcf} 2> {log}"

rule samtools_split:
    input:
        bam = "ngmlr_alignment/{sample}.bam",
        bai = "ngmlr_alignment/{sample}.bam.bai",
    output:
        temp("split_ngmlr_alignment/{sample}-{chromosome}.bam")
    params:
        chrom = "{chromosome}"
    log:
        "logs/samtools_split/{sample}-{chromosome}.log"
    shell:
        "samtools view {input.bam} {params.chrom} -o {output} 2> {log}"

rule split_bed:
    input:
        bed = config["annotbed"]
    output:
        expand("split_annotation_bed/{chromosome}.bed", chromosome=CHROMOSOMES)
    log:
        "logs/split_bed/split_annotation.log"
    shell:
        '''awk '{print $0 >> $1".bed"}' {input.bed} 2> {log}'''


rule nanosv_call:
    input:
        bam = "split_ngmlr_alignment/{sample}-{chromosome}.bam",
        bai = "split_ngmlr_alignment/{sample}-{chromosome}.bam.bai",
        bed = config["annotbed"]
        # bed = "split_annotation_bed/{chromosome}.bed"
    output:
        temp("split_nanosv_genotypes/{sample}-{chromosome}.vcf")
    params:
        samtools = "samtools"
    log:
        "logs/nanosv/{sample}-{chromosome}.log"
    shell:
        "NanoSV --bed {input.bed} -s {params.samtools} {input.bam} -o {output} 2> {log}"


rule cat_vcfs:
    input:
        expand("split_nanosv_genotypes/{{sample}}-{chromosome}.vcf", chromosome=CHROMOSOMES)
    output:
        "nanosv_genotypes/{sample}.vcf"
    log:
        "logs/vcf-concat/{sample}.log"
    shell:
        "vcf-concat {input} > {output} 2> {log}"

rule survivor:
    input:
        expand("{{caller}}_{{stage}}/{sample}.vcf", sample=config["samples"])
    output:
        vcf = temp("{caller}_combined/{stage}.vcf"),
        fofn = temp("{caller}_{stage}/samples.fofn")
    params:
        distance = config["parameters"]["survivor_distance"],
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        "logs/{caller}/surivor_{stage}.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"

rule survivor_all:
    input:
        expand("{caller}_genotypes/{sample}.vcf",
               sample=config["samples"],
               caller=["sniffles", "nanosv"])
    output:
        vcf = temp("all_combined/genotypes.vcf"),
        fofn = temp("all_combined/samples.fofn")
    params:
        distance = config["parameters"]["survivor_distance"],
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        "logs/all/surivor.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"

rule survivor_pairwise:
    input:
        expand("{caller}_genotypes/{{sample}}.vcf", caller=["sniffles", "nanosv"])
    output:
        vcf = temp("high_confidence/{sample}.vcf"),
        fofn = temp("high_confidence/{sample}.fofn")
    params:
        distance = config["parameters"]["survivor_distance"],
        caller_support = 2,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        "logs/high_confidence/surivor_pairwise_{sample}.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"

rule survivor_high_confidence:
    input:
        expand("high_confidence/{sample}.vcf", sample=config["samples"])
    output:
        vcf = temp("high_confidence_combined/genotypes.vcf"),
        fofn = temp("high_confidence_combined/samples.fofn")
    params:
        distance = config["parameters"]["survivor_distance"],
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        "logs/high_confidence_combined/surivor_high_confidence.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"


rule mosdepth_get:
    input:
        bam = "ngmlr_alignment/{sample}.bam",
        bai = "ngmlr_alignment/{sample}.bam.bai"
    threads: 4
    output:
        protected("mosdepth/{sample}.mosdepth.global.dist.txt"),
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
        os.path.join(workflow.basedir, "scripts/combine_mosdepth.py") + \
            " {input} -o {output} 2> {log}"


rule mosdepth_global_plot:
    input:
        expand("mosdepth/{sample}.mosdepth.global.dist.txt", sample=config["samples"])
    output:
        "mosdepth_global_plot/global.html"
    log:
        "logs/mosdepth/mosdepth_global_plot.log"
    shell:
        os.path.join(workflow.basedir, "scripts/mosdepth_plot-dist.py") + \
            " {input} -o {output} 2> {log}"


rule SV_length_plot:
    input:
        "{caller}_{stage}/{sample}.vcf"
    output:
        "SV-plots/SV-length_{caller}_{stage}_{sample}.png"
    log:
        "logs/svplot/svlength_{caller}_{stage}_{sample}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-length-plot.py") + " {input} {output} 2> {log}"


rule SV_plot_carriers:
    input:
        "{caller}_combined/annot_genotypes.vcf"
    output:
        "SV-plots/SV-{caller}_carriers.png"
    log:
        "logs/svplot/svcarriers_{caller}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-carriers-plot.py") + \
            " {input} {output} 2> {log}"


rule sort_vcf:
    input:
        "{caller}_combined/genotypes.vcf"
    output:
        temp("{caller}_combined/sorted_genotypes.vcf")
    log:
        "logs/sort_vcf/sorting_{caller}.log"
    threads: 8
    shell:
        "vcf-sort {input} > {output} 2> {log}"


rule annotate_vcf:
    input:
        "{caller}_combined/sorted_genotypes.vcf"
    output:
        "{caller}_combined/annot_genotypes.vcf"
    log:
        "logs/annotate_vcf/annotate_{caller}.log"
    params:
        conf = "/home/wdecoster/projects/SV-snakemake/configuration/vcfanno_conf.toml"
    threads: 8
    shell:
        "vcfanno -ends -p {threads} {params.conf} {input} > {output} 2> {log}"
