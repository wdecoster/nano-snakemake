rule survivor:
    input:
        [f"{OUTDIR}/{{aligner}}/{{caller}}_{{stage}}/{sample}.vcf" for sample in config["samples"]]
    output:
        vcf = temp(f"{OUTDIR}/{{aligner}}/{{caller}}_combined/{{stage}}.vcf"),
        fofn = temp(f"{OUTDIR}/{{aligner}}/{{caller}}_{{stage}}/samples.fofn")
    threads: get_resource("survivor", "threads")
    resources:
        mem=get_resource("survivor", "mem"),
        walltime=get_resource("survivor", "walltime")
    params:
        distance = config["parameters"]["survivor_distance"],
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    conda: "../envs/survivor.yaml"
    log:
        f"{LOGDIR}/{{aligner}}/{{caller}}/survivor_{{stage}}.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"

rule survivor_all:
    input:
        [f"{OUTDIR}/{{aligner}}/{caller}_genotypes/{sample}.vcf"
               for sample in config["samples"]
               for caller in ["sniffles", "svim", "nanosv"]]
    output:
        vcf = temp(f"{OUTDIR}/{{aligner}}/all_combined/genotypes.vcf"),
        fofn = temp(f"{OUTDIR}/{{aligner}}/all_combined/samples.fofn")
    threads: get_resource("survivor_all", "threads")
    resources:
        mem=get_resource("survivor_all", "mem"),
        walltime=get_resource("survivor_all", "walltime")
    params:
        distance = config["parameters"]["survivor_distance"],
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    conda: "../envs/survivor.yaml"
    log:
        f"{LOGDIR}/{{aligner}}/all/survivor.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"


rule bgzip_and_tabix:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_genotypes/{{sample}}.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/{{caller}}_genotypes/{{sample}}.vcf.gz"
    threads: get_resource("bgzip_and_tabix", "threads")
    resources:
        mem=get_resource("bgzip_and_tabix", "mem"),
        walltime=get_resource("bgzip_and_tabix", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/bgzip-tabix/{{caller}}-{{sample}}.log"
    shell:
        """
        vcf-sort {input} | bgzip > {output} 2> {log} &&
        tabix -p vcf {output} 2>> {log}
        """
