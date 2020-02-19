import os
import gzip
import sys

configfile: "config.yaml"

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

include: "rules/common.smk"
include: "rules/mosdepth.smk"
include: "rules/plots.smk"
include: "rules/align.smk"
include: "rules/survivor.smk"
include: "rules/callers.smk"
include: "rules/vcf.smk"

# Target rules #

rule fast:
    input:
        expand(f"{OUTDIR}/minimap2/SV-plots/SV-length_{{caller}}_genotypes_{{sample}}.png",
               sample=config["samples"],
               caller=["sniffles"]),
        expand(f"{OUTDIR}/minimap2/SV-plots/SV-{{caller}}_carriers.png",
               caller=["sniffles"]),
        expand(f"{OUTDIR}/minimap2/{{caller}}_combined/annot_genotypes.vcf",
               caller=["sniffles", "svim"]),
        expand(f"{OUTDIR}/minimap2/alignment_stats/{{sample}}.txt",
               sample=config["samples"]),
        expand(f"{OUTDIR}/minimap2/npinv/{{sample}}.vcf",
               sample=config["samples"]),

rule precise:
    input:
        expand(f"{OUTDIR}/ngmlr/SV-plots/SV-length_{{caller}}_genotypes_{{sample}}.png",
               sample=config["samples"],
               caller=["sniffles"]),
        expand(f"{OUTDIR}/ngmlr/SV-plots/SV-{{caller}}_carriers.png",
               caller=["sniffles"]),
        expand(f"{OUTDIR}/ngmlr/{{caller}}_combined/annot_genotypes.vcf",
               caller=["sniffles", "svim"]),
        expand(f"{OUTDIR}/ngmlr/alignment_stats/{{sample}}.txt",
               sample=config["samples"]),
        expand(f"{OUTDIR}/ngmlr/npinv/{{sample}}.vcf",
               sample=config["samples"]),
        f"{OUTDIR}/ngmlr/mosdepth/regions.combined.gz",
        f"{OUTDIR}/ngmlr/mosdepth_global_plot/global.html",


rule minimap2:
    input:
        expand(f"{OUTDIR}/minimap2/SV-plots/SV-length_{{caller}}_genotypes_{{sample}}.png",
               sample=config["samples"],
               caller=["sniffles", "nanosv"]),
        expand(f"{OUTDIR}/minimap2/SV-plots/SV-{{caller}}_carriers.png",
               caller=["sniffles", "nanosv"]),
        expand(f"{OUTDIR}/minimap2/{{caller}}_combined/annot_genotypes.vcf",
               caller=["sniffles", "nanosv", "svim"]),
        expand(f"{OUTDIR}/minimap2/alignment_stats/{{sample}}.txt",
               sample=config["samples"]),
        f"{OUTDIR}/minimap2/all_combined/annot_genotypes.vcf",
        f"{OUTDIR}/minimap2/mosdepth/regions.combined.gz",
        f"{OUTDIR}/minimap2/mosdepth_global_plot/global.html",
        expand(f"{OUTDIR}/minimap2/npinv/{{sample}}.vcf",
               sample=config["samples"]),

rule minimap2_pbsv:
    input:
        expand(f"{OUTDIR}/minimap2_pbsv/SV-plots/SV-length_{{caller}}_genotypes_{{sample}}.png",
               sample=config["samples"],
               caller=["sniffles", "nanosv", "pbsv"]),
        expand(f"{OUTDIR}/minimap2_pbsv/SV-plots/SV-{{caller}}_carriers.png",
               caller=["sniffles", "nanosv", "pbsv"]),
        expand(f"{OUTDIR}/minimap2_pbsv/{{caller}}_combined/annot_genotypes.vcf",
               caller=["sniffles", "nanosv", "svim", "pbsv"]),
        expand(f"{OUTDIR}/minimap2_pbsv/alignment_stats/{{sample}}.txt",
               sample=config["samples"]),
        f"{OUTDIR}/minimap2_pbsv/all_combined/annot_genotypes.vcf",
        f"{OUTDIR}/minimap2_pbsv/mosdepth/regions.combined.gz",
        f"{OUTDIR}/minimap2_pbsv/mosdepth_global_plot/global.html",
        expand(f"{OUTDIR}/minimap2_pbsv/npinv/{{sample}}.vcf",
               sample=config["samples"]),

rule ngmlr:
    input:
        expand(f"{OUTDIR}/ngmlr/SV-plots/SV-length_{{caller}}_genotypes_{{sample}}.png",
               sample=config["samples"],
               caller=["sniffles", "nanosv"]),
        expand(f"{OUTDIR}/ngmlr/SV-plots/SV-{{caller}}_carriers.png",
               caller=["sniffles", "nanosv"]),
        expand(f"{OUTDIR}/ngmlr/{{caller}}_combined/annot_genotypes.vcf",
               caller=["sniffles", "nanosv", "svim"]),
        expand(f"{OUTDIR}/ngmlr/alignment_stats/{{sample}}.txt",
               sample=config["samples"]),
        expand(f"{OUTDIR}/ngmlr/npinv/{{sample}}.vcf",
               sample=config["samples"]),
        f"{OUTDIR}/ngmlr/all_combined/annot_genotypes.vcf",
        f"{OUTDIR}/ngmlr/mosdepth/regions.combined.gz",
        f"{OUTDIR}/ngmlr/mosdepth_global_plot/global.html",

rule last_prepare:
    input:
        f"{OUTDIR}/last/index/last-train.params"

rule last:
    input:
        f"{OUTDIR}/last/tandem_genotypes_reformatted/combined.txt"
