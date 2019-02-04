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
        "logs/{aligner}/{caller}/survivor_{stage}.log"
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
        "logs/{aligner}/all/survivor.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"

rule bgzip_and_tabix:
    input:
        "{aligner}/{caller}_genotypes/{sample}.vcf"
    output:
        "{aligner}/{caller}_genotypes/{sample}.vcf.gz"
    log:
        "logs/{aligner}/bgzip-tabix/{caller}-{sample}.log"
    shell:
        """
        vcf-sort {input} | bgzip > {output} 2> {log} &&
        tabix -p vcf {output} 2>> {log}
        """
