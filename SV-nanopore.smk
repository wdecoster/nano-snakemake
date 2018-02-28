configfile: "config.yaml"


def get_samples():
    return config["fqfolder"]


rule all:
    input:
        "ngmlr_alignment/" + config["samplename"] + ".bam.bai",
        protected("sniffles_calls/" + config["samplename"] + ".vcf")

rule ngmlr:
    input:
        fq = get_samples(),
        genome = "/home/wdecoster/databases/Homo_sapiens/genome_hg19.fa"
    output:  protected("ngmlr_alignment/" + config["samplename"] + ".bam")
    threads: 48
    log:
        "logs/ngmlr/" + config["samplename"] + ".log"
    shell:
        "zcat {input.fq}/*.fastq.gz | \
         ngmlr -x ont -t {threads} -r {input.genome} | \
         samtools sort -@ {threads} -o {output} - 2> {log}"


rule samtools_index:
    input:
        "ngmlr_alignment/" + config["samplename"] + ".bam"
    output:
        "ngmlr_alignment/" + config["samplename"] + ".bam.bai"
    threads: 24
    log:
        "logs/samtools_index/" + config["samplename"] + ".log"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"


rule sniffles:
    input:
        "ngmlr_alignment/" + config["samplename"] + ".bam"
    output:
        protected("sniffles_calls/" + config["samplename"] + ".bedpe")
    threads: 48
    log:
        "logs/sniffles/" + config["samplename"] + ".log"
    shell:
        "sniffles --mapped_reads {input} --bedpe {output} --threads {threads} 2> {log}"
