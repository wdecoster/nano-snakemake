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

##### Target rules #####

rule minimap2:
    input:
        expand("minimap2/SV-plots/SV-length_{caller}_genotypes_{sample}.png",
               sample=config["samples"],
               caller=["sniffles", "nanosv"]),
        expand("minimap2/SV-plots/SV-{caller}_carriers.png",
               caller=["sniffles", "nanosv"]),
        expand("minimap2/{caller}_combined/annot_genotypes.vcf",
               caller=["sniffles", "nanosv"]),
        expand("minimap2/alignment_stats/{sample}.txt",
               sample=config["samples"]),
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
        expand("minimap2_last_like/SV-plots/SV-{caller}_carriers.png",
               caller=["sniffles", "nanosv"]),
        expand("minimap2_last_like/{caller}_combined/annot_genotypes.vcf",
               caller=["sniffles", "nanosv"]),
        expand("minimap2_last_like/alignment_stats/{sample}.txt",
               sample=config["samples"]),
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
        expand("ngmlr/SV-plots/SV-{caller}_carriers.png",
               caller=["sniffles", "nanosv"]),
        expand("ngmlr/{caller}_combined/annot_genotypes.vcf",
               caller=["sniffles", "nanosv"]),
        expand("ngmlr/alignment_stats/{sample}.txt",
               sample=config["samples"]),
        "ngmlr/all_combined/annot_genotypes.vcf",
        "ngmlr/high_confidence_combined/annot_genotypes.vcf",
        "ngmlr/high_sensitivity_combined/annot_genotypes.vcf",
        "ngmlr/mosdepth/regions.combined.gz",
        "ngmlr/mosdepth_global_plot/global.html",
