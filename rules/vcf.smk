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
        f"{OUTDIR}/{{aligner}}/split_nanosv_genotypes/{{sample}}-{{chromosome}}.vcf",
    output:
        vcf = temp(f"{OUTDIR}/{{aligner}}/split_nanosv_genotypes_renamed/{{sample}}-{{chromosome}}.vcf"),
        sample = temp(f"{OUTDIR}/{{aligner}}/split_nanosv_genotypes_renamed/sample_{{sample}}-{{chromosome}}.txt")
    threads: get_resource("bcftools_reheader", "threads")
    resources:
        mem=get_resource("bcftools_reheader", "mem"),
        walltime=get_resource("bcftools_reheader", "walltime")
    log:
        "logs/{{aligner}}/bcftools_reheader/{{sample}}-{{chromosome}}.log"
    conda: "../envs/bcftools.yaml"
    shell:
        """
        echo {wildcards.sample} > {output.sample} &&
        bcftools reheader -s {output.sample} {input} -o {output.vcf} 2> {log}
        """

rule bcftools_reheader_sniffles:
    """Rule to be deleted as soon as ngmlr uses read groups correctly"""
    input:
        f"{OUTDIR}/{{aligner}}/sniffles_genotypes_temp/{{sample}}.vcf"
    output:
        vcf = f"{OUTDIR}/{{aligner}}/sniffles_genotypes/{{sample}}.vcf",
        sample = temp(f"{OUTDIR}/{{aligner}}/sniffles_genotypes/sample_{{sample}}.txt")
    threads: get_resource("bcftools_reheader_sniffles", "threads")
    resources:
        mem=get_resource("bcftools_reheader_sniffles", "mem"),
        walltime=get_resource("bcftools_reheader_sniffles", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/bcftools_reheader/{{sample}}.log"
    conda: "../envs/bcftools.yaml"
    shell:
        """
        echo {wildcards.sample} > {output.sample} &&
        bcftools reheader -s {output.sample} {input} -o {output.vcf} 2> {log}
        """

rule cat_vcfs:
    input:
        [f"{OUTDIR}/{{aligner}}/split_nanosv_genotypes_renamed/{{sample}}-{chromosome}.vcf" for chromosome in CHROMOSOMES]
    output:
        f"{OUTDIR}/{{aligner}}/nanosv_genotypes/{{sample}}.vcf"
    threads: get_resource("cat_vcfs", "threads")
    resources:
        mem=get_resource("cat_vcfs", "mem"),
        walltime=get_resource("cat_vcfs", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/bcftools-concat/{{sample}}.log"
    conda: "../envs/bcftools.yaml"
    shell:
        "bcftools concat {input} | bcftools sort - -o {output} 2> {log}"

rule sort_vcf:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_combined/genotypes.vcf"
    output:
        temp(f"{OUTDIR}/{{aligner}}/{{caller}}_combined/sorted_genotypes.vcf")
    threads: get_resource("sort_vcf", "threads")
    resources:
        mem=get_resource("sort_vcf", "mem"),
        walltime=get_resource("sort_vcf", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/bcftools_sort/sorting_{{caller}}.log"
    conda: "../envs/bcftools.yaml"
    shell:
        "bcftools sort {input} > {output} 2> {log}"


rule annotate_vcf:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_combined/sorted_genotypes.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/{{caller}}_combined/annot_genotypes.vcf"
    log:
        f"{LOGDIR}/{{aligner}}/annotate_vcf/annotate_{{caller}}.log"
    params:
        conf = config["vcfanno_conf"],
    threads: get_resource("annotate_vcf", "threads")
    resources:
        mem=get_resource("annotate_vcf", "mem"),
        walltime=get_resource("annotate_vcf", "walltime")
    conda: "../envs/vcfanno.yaml"
    shell:
        "vcfanno -ends -p {threads} {params.conf} {input} > {output} 2> {log}"
