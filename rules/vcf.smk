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
        bcftools reheader -s {output.sample} {input} -o {output.vcf} 2> {log}
        """

rule bcftools_reheader_sniffles:
    """Rule to be deleted as soon as ngmlr uses read groups correctly"""
    input:
        "{aligner}/sniffles_genotypes_temp/{sample}.vcf"
    output:
        vcf = temp("{aligner}/sniffles_genotypes/{sample}.vcf"),
        sample = temp("{aligner}/sniffles_genotypes/sample_{sample}.txt")
    params:
        sample = "{sample}"
    log:
        "logs/{aligner}/bcftools_reheader/{sample}.log"
    shell:
        """
        echo {params.sample} > {output.sample} &&
        bcftools reheader -s {output.sample} {input} -o {output.vcf} 2> {log}
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
