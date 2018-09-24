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
        "logs/{aligner}/high_confidence/survivor_pairwise_{sample}.log"
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
        "logs/{aligner}/high_confidence_combined/survivor_high_confidence.log"
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
        sort -k1,1d -k2,2n {input} | bgzip > {output} 2> {log} &&
        tabix -p vcf {output} 2>> {log}
        """

rule survivor_pairwise_high_sensitivity:
    input:
        expand("{{aligner}}/{caller}_genotypes/{{sample}}.vcf.gz", caller=["sniffles", "nanosv"])
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
        "logs/{aligner}/high_sensitivity/survivor_pairwise_{sample}.log"
    shell:
        """
        bcftools concat -a {input} | bcftools sort - -o {output.vcf_unmerged} 2> {log}; \
        ls {output.vcf_unmerged} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2>> {log}
        """

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
        "logs/{aligner}/high_sensitivity_combined/survivor_high_sensitivity.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"
