import os
import gzip
import sys

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
if not CHROMOSOMES:
    sys.exit("\n\nUnexpectedly no chromosomes were found for SV calling.\n"
             "Is your annotation matching to your fasta file?\n\n")


rule minimap2:
    input:
        expand("minimap2/SV-plots/SV-length_{caller}_genotypes_{sample}.png",
               sample=config["samples"],
               caller=["sniffles", "nanosv"]),
        expand("minimap2/SV-plots/SV-{caller}_carriers.png", caller=["sniffles", "nanosv"]),
        expand("minimap2/{caller}_combined/annot_genotypes.vcf", caller=["sniffles", "nanosv"]),
        expand("minimap2/alignment_stats/{sample}.txt", sample=config["samples"]),
        "minimap2/all_combined/annot_genotypes.vcf",
        "minimap2/high_confidence_combined/annot_genotypes.vcf",
        "minimap2/high_sensitivity_combined/annot_genotypes.vcf",
        "minimap2/mosdepth/regions.combined.gz",
        "minimap2/mosdepth_global_plot/global.html",


rule minimap2_last_like:
    input:
        expand("minimap2_last_like/SV-plots/SV-length_{caller}_genotypes_{sample}.png",
               sample=config["samples"],
               caller=["sniffles", "nanosv"]),
        expand(
            "minimap2_last_like/SV-plots/SV-{caller}_carriers.png", caller=["sniffles", "nanosv"]),
        expand("minimap2_last_like/{caller}_combined/annot_genotypes.vcf",
               caller=["sniffles", "nanosv"]),
        expand("minimap2_last_like/alignment_stats/{sample}.txt", sample=config["samples"]),
        "minimap2_last_like/all_combined/annot_genotypes.vcf",
        "minimap2_last_like/high_confidence_combined/annot_genotypes.vcf",
        "minimap2_last_like/high_sensitivity_combined/annot_genotypes.vcf",
        "minimap2_last_like/mosdepth/regions.combined.gz",
        "minimap2_last_like/mosdepth_global_plot/global.html",

rule ngmlr:
    input:
        expand("ngmlr/SV-plots/SV-length_{caller}_genotypes_{sample}.png",
               sample=config["samples"],
               caller=["sniffles", "nanosv"]),
        expand("ngmlr/SV-plots/SV-{caller}_carriers.png", caller=["sniffles", "nanosv"]),
        expand("ngmlr/{caller}_combined/annot_genotypes.vcf", caller=["sniffles", "nanosv"]),
        expand("ngmlr/alignment_stats/{sample}.txt", sample=config["samples"]),
        "ngmlr/all_combined/annot_genotypes.vcf",
        "ngmlr/high_confidence_combined/annot_genotypes.vcf",
        "ngmlr/high_sensitivity_combined/annot_genotypes.vcf",
        "ngmlr/mosdepth/regions.combined.gz",
        "ngmlr/mosdepth_global_plot/global.html",


#######################################################################################

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
        "minimap2 --MD -a -t {threads} {input.genome} {input.fq}/*.fastq.gz | \
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
        "minimap2 --MD -a  --no-long-join -r50 -t {threads} {input.genome} {input.fq}/*.fastq.gz | \
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

rule sniffles_call:
    input:
        "{aligner}/alignment/{sample}.bam"
    output:
        "{aligner}/sniffles_calls/{sample}.vcf"
    threads: 8
    log:
        "logs/{aligner}/sniffles_call/{sample}.log"
    shell:
        "sniffles --mapped_reads {input} --vcf {output} --threads {threads} 2> {log}"

rule sniffles_genotype:
    input:
        bam = "{aligner}/alignment/{sample}.bam",
        ivcf = "{aligner}/sniffles_combined/calls.vcf"
    output:
        "{aligner}/sniffles_genotypes/{sample}.vcf"
    threads: 8
    log:
        "logs/{aligner}/sniffles_genotype/{sample}.log"
    shell:
        "sniffles --mapped_reads {input.bam} \
                  --vcf {output} \
                  --threads {threads} \
                  --Ivcf {input.ivcf} 2> {log}"

rule samtools_split:
    input:
        bam = "{aligner}/alignment/{sample}.bam",
        bai = "{aligner}/alignment/{sample}.bam.bai",
    output:
        temp("{aligner}/alignment/{sample}-{chromosome}.bam")
    params:
        chrom = "{chromosome}"
    log:
        "logs/{aligner}/samtools_split/{sample}-{chromosome}.log"
    shell:
        "samtools view {input.bam} {params.chrom} -o {output} 2> {log}"


rule nanosv_call:
    '''

    call variants using NanoSV on separate chromosomes
    the shell command will first check if there are reads in this chromosome
    and if not, will just touch the output and leave it empty
    without raising an error
    '''
    input:
        bam = "{aligner}/alignment/{sample}-{chromosome}.bam",
        bai = "{aligner}/alignment/{sample}-{chromosome}.bam.bai",
        bed = config["annotbed"]
    output:
        temp("{aligner}/split_nanosv_genotypes/{sample}-{chromosome}.vcf")
    params:
        samtools = "samtools"
    threads:
        2
    log:
        "logs/{aligner}/nanosv/{sample}-{chromosome}.log"
    shell:
        """
        reads=$(samtools idxstats {input.bam} | awk 'BEGIN {{FS = "\\t"}} ; {{sum+=$3}} END {{print sum}}')
        if [ "$reads" -eq "0" ]; then
            echo "##fileformat=VCFv4.1" > {output} && \
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> {output}
            echo "NanoSV: No reads in {input.bam}" >> exceptions.txt 2> {log}
        else
            NanoSV --bed {input.bed} \
                    --threads {threads} \
                    --sambamba {params.samtools} {input.bam} \
                    -o {output} 2> {log}
        fi
        """

rule bcftools_reheader:
    input:
        "{aligner}/split_nanosv_genotypes/{sample}-{chromosome}.vcf",
    output:
        vcf = temp("{aligner}/split_nanosv_genotypes_renamed/{sample}-{chromosome}.vcf"),
        sample = temp("{aligner}/split_nanosv_genotypes_renamed/sample_{sample}-{chromosome}.txt")
    params:
        sample = "{sample}"
    log:
        "logs/{aligner}/bcftools_reheader/{sample}-{chromosome}.log"
    shell:
        """
        echo {params.sample} > {output.sample} &&
        bcftools reheader -s {output.sample} {input} > {output.vcf} 2> {log}
        """

rule cat_vcfs:
    input:
        expand("{{aligner}}/split_nanosv_genotypes_renamed/{{sample}}-{chromosome}.vcf",
               chromosome=CHROMOSOMES)
    output:
        "{aligner}/nanosv_genotypes/{sample}.vcf"
    log:
        "logs/{aligner}/bcftools-concat/{sample}.log"
    shell:
        "bcftools concat {input} | bcftools sort - -o {output} 2> {log}"

rule survivor:
    input:
        expand("{{aligner}}/{{caller}}_{{stage}}/{sample}.vcf", sample=config["samples"])
    output:
        vcf = temp("{aligner}/{caller}_combined/{stage}.vcf"),
        fofn = temp("{aligner}/{caller}_{stage}/samples.fofn")
    params:
        distance = config["parameters"]["survivor_distance"],
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        "logs/{aligner}/{caller}/surivor_{stage}.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"

rule survivor_all:
    input:
        expand("{{aligner}}/{caller}_genotypes/{sample}.vcf",
               sample=config["samples"],
               caller=["sniffles", "nanosv"])
    output:
        vcf = temp("{aligner}/all_combined/genotypes.vcf"),
        fofn = temp("{aligner}/all_combined/samples.fofn")
    params:
        distance = config["parameters"]["survivor_distance"],
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        "logs/{aligner}/all/surivor.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"

rule survivor_pairwise_high_confidence:
    input:
        expand("{{aligner}}/{caller}_genotypes/{{sample}}.vcf", caller=["sniffles", "nanosv"])
    output:
        vcf = temp("{aligner}/high_confidence/{sample}.vcf"),
        fofn = temp("{aligner}/high_confidence/{sample}.fofn")
    params:
        distance = config["parameters"]["survivor_high_confidence_distance"],
        caller_support = 2,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        "logs/{aligner}/high_confidence/surivor_pairwise_{sample}.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"

rule survivor_combine_high_confidence:
    input:
        expand("{{aligner}}/high_confidence/{sample}.vcf", sample=config["samples"])
    output:
        vcf = temp("{aligner}/high_confidence_combined/genotypes.vcf"),
        fofn = temp("{aligner}/high_confidence_combined/samples.fofn")
    params:
        distance = config["parameters"]["survivor_high_confidence_distance"],
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        "logs/{aligner}/high_confidence_combined/surivor_high_confidence.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"

rule survivor_pairwise_high_sensitivity:
    input:
        expand("{{aligner}}/{caller}_genotypes/{{sample}}.vcf", caller=["sniffles", "nanosv"])
    output:
        vcf = temp("{aligner}/high_sensitivity/{sample}.vcf"),
        vcf_unmerged = temp("{aligner}/high_sensitivity/{sample}_unmerged.vcf"),
        fofn = temp("{aligner}/high_sensitivity/{sample}.fofn")
    params:
        distance = config["parameters"]["survivor_high_sensitivity_distance"],
        caller_support = 1,
        same_type = -1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        "logs/{aligner}/high_sensitivity/surivor_pairwise_{sample}.log"
    shell:
        "vcf-concat {input} | vcf-sort > {output.vcf_unmerged} ; \
        ls {output.vcf_unmerged} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"

rule survivor_combine_high_sensitivity:
    input:
        expand("{{aligner}}/high_sensitivity/{sample}.vcf", sample=config["samples"])
    output:
        vcf = temp("{aligner}/high_sensitivity_combined/genotypes.vcf"),
        fofn = temp("{aligner}/high_sensitivity_combined/samples.fofn")
    params:
        distance = config["parameters"]["survivor_high_sensitivity_distance"],
        caller_support = 1,
        same_type = -1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        "logs/{aligner}/high_sensitivity_combined/surivor_high_sensitivity.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"


rule mosdepth_get:
    input:
        bam = "{aligner}/alignment/{sample}.bam",
        bai = "{aligner}/alignment/{sample}.bam.bai"
    threads: 4
    output:
        protected("{aligner}/mosdepth/{sample}.mosdepth.global.dist.txt"),
        protected("{aligner}/mosdepth/{sample}.regions.bed.gz"),
    params:
        windowsize = 500,
        prefix = "{sample}",
        aligner = "{aligner}"
    log:
        "logs/{aligner}/mosdepth/mosdepth_{sample}.log"
    shell:
        "mosdepth --threads {threads} \
                  -n \
                  --by {params.windowsize} \
                  {params.aligner}/mosdepth/{params.prefix} {input.bam} 2> {log}"


rule mosdepth_combine:
    input:
        expand("{{aligner}}/mosdepth/{sample}.regions.bed.gz", sample=config["samples"])
    output:
        "{aligner}/mosdepth/regions.combined.gz"
    log:
        "logs/{aligner}/mosdepth/mosdepth_combine.log"
    shell:
        os.path.join(workflow.basedir, "scripts/combine_mosdepth.py") + \
            " {input} -o {output} 2> {log}"


rule mosdepth_global_plot:
    input:
        expand("{{aligner}}/mosdepth/{sample}.mosdepth.global.dist.txt", sample=config["samples"])
    output:
        "{aligner}/mosdepth_global_plot/global.html"
    log:
        "logs/{aligner}/mosdepth/mosdepth_global_plot.log"
    shell:
        os.path.join(workflow.basedir, "scripts/mosdepth_plot-dist.py") + \
            " {input} -o {output} 2> {log}"


rule SV_length_plot:
    input:
        "{aligner}/{caller}_{stage}/{sample}.vcf"
    output:
        plot = "{aligner}/SV-plots/SV-length_{caller}_{stage}_{sample}.png",
        counts = "{aligner}/SV-plots/SV-nucleotides_affected_{caller}_{stage}_{sample}.txt",
    log:
        "logs/{aligner}/svplot/svlength_{caller}_{stage}_{sample}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-length-plot.py") + \
            " {input} --output {output.plot} --counts {output.counts} 2> {log}"


rule SV_plot_carriers:
    input:
        "{aligner}/{caller}_combined/annot_genotypes.vcf"
    output:
        "{aligner}/SV-plots/SV-{caller}_carriers.png"
    log:
        "logs/{aligner}/svplot/svcarriers_{caller}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-carriers-plot.py") + \
            " {input} {output} 2> {log}"


rule sort_vcf:
    input:
        "{aligner}/{caller}_combined/genotypes.vcf"
    output:
        temp("{aligner}/{caller}_combined/sorted_genotypes.vcf")
    log:
        "logs/{aligner}/bcftools_sort/sorting_{caller}.log"
    threads: 8
    shell:
        "bcftools sort {input} -o {output} 2> {log}"


rule annotate_vcf:
    input:
        "{aligner}/{caller}_combined/sorted_genotypes.vcf"
    output:
        "{aligner}/{caller}_combined/annot_genotypes.vcf"
    log:
        "logs/{aligner}/annotate_vcf/annotate_{caller}.log"
    params:
        conf = config["vcfanno_conf"],
    threads: 8
    shell:
        "vcfanno -ends -p {threads} {params.conf} {input} > {output} 2> {log}"
