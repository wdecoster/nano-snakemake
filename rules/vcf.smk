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
