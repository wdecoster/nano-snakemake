import os
import gzip
import sys

include: "rules/mosdepth.smk"
include: "rules/plots.smk"
include: "rules/align.smk"
include: "rules/survivor.smk"
include: "rules/callers.smk"
include: "rules/vcf.smk"

configfile: "config.yaml"

# Target rules #

rule fast:
    input:
        expand("minimap2/SV-plots/SV-length_{caller}_genotypes_{sample}.png",
               sample=config["samples"],
               caller=["sniffles"]),
        expand("minimap2/SV-plots/SV-{caller}_carriers.png",
               caller=["sniffles"]),
        expand("minimap2/{caller}_combined/annot_genotypes.vcf",
               caller=["sniffles", "svim"]),
        expand("minimap2/alignment_stats/{sample}.txt",
               sample=config["samples"]),
        expand("minimap2/npinv/{sample}.vcf",
               sample=config["samples"]),

rule precise:
    input:
        expand("ngmlr/SV-plots/SV-length_{caller}_genotypes_{sample}.png",
               sample=config["samples"],
               caller=["sniffles"]),
        expand("ngmlr/SV-plots/SV-{caller}_carriers.png",
               caller=["sniffles"]),
        expand("ngmlr/{caller}_combined/annot_genotypes.vcf",
               caller=["sniffles", "svim"]),
        expand("ngmlr/alignment_stats/{sample}.txt",
               sample=config["samples"]),
        expand("ngmlr/npinv/{sample}.vcf",
               sample=config["samples"]),
        "ngmlr/mosdepth/regions.combined.gz",
        "ngmlr/mosdepth_global_plot/global.html",


rule minimap2:
    input:
        expand("minimap2/SV-plots/SV-length_{caller}_genotypes_{sample}.png",
               sample=config["samples"],
               caller=["sniffles", "nanosv"]),
        expand("minimap2/SV-plots/SV-{caller}_carriers.png",
               caller=["sniffles", "nanosv"]),
        expand("minimap2/{caller}_combined/annot_genotypes.vcf",
               caller=["sniffles", "nanosv", "svim"]),
        expand("minimap2/alignment_stats/{sample}.txt",
               sample=config["samples"]),
        "minimap2/all_combined/annot_genotypes.vcf",
        "minimap2/mosdepth/regions.combined.gz",
        "minimap2/mosdepth_global_plot/global.html",
        expand("minimap2/npinv/{sample}.vcf",
               sample=config["samples"]),

rule minimap2_pbsv:
    input:
        expand("minimap2_pbsv/SV-plots/SV-length_{caller}_genotypes_{sample}.png",
               sample=config["samples"],
               caller=["sniffles", "nanosv", "pbsv"]),
        expand("minimap2_pbsv/SV-plots/SV-{caller}_carriers.png",
               caller=["sniffles", "nanosv", "pbsv"]),
        expand("minimap2_pbsv/{caller}_combined/annot_genotypes.vcf",
               caller=["sniffles", "nanosv", "svim", "pbsv"]),
        expand("minimap2_pbsv/alignment_stats/{sample}.txt",
               sample=config["samples"]),
        "minimap2_pbsv/all_combined/annot_genotypes.vcf",
        "minimap2_pbsv/mosdepth/regions.combined.gz",
        "minimap2_pbsv/mosdepth_global_plot/global.html",
        expand("minimap2_pbsv/npinv/{sample}.vcf",
               sample=config["samples"]),

rule ngmlr:
    input:
        expand("ngmlr/SV-plots/SV-length_{caller}_genotypes_{sample}.png",
               sample=config["samples"],
               caller=["sniffles", "nanosv"]),
        expand("ngmlr/SV-plots/SV-{caller}_carriers.png",
               caller=["sniffles", "nanosv"]),
        expand("ngmlr/{caller}_combined/annot_genotypes.vcf",
               caller=["sniffles", "nanosv", "svim"]),
        expand("ngmlr/alignment_stats/{sample}.txt",
               sample=config["samples"]),
        expand("ngmlr/npinv/{sample}.vcf",
               sample=config["samples"]),
        "ngmlr/all_combined/annot_genotypes.vcf",
        "ngmlr/mosdepth/regions.combined.gz",
        "ngmlr/mosdepth_global_plot/global.html",

rule last_prepare:
    input:
        "last/index/last-train.params"

rule last:
    input:
        "last/tandem_genotypes_reformatted/combined.txt"
