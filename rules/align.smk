def get_samples(wildcards):
    return config["samples"][wildcards.sample]

def get_all_samples(wildcards):
    return config["samples"].values()

rule minimap2_align:
    input:
        fq = get_samples,
        genome = config["genome"]
    output:
        f"{OUTDIR}/minimap2/alignment/{{sample}}.bam"
    threads: get_resource("minimap2_align", "threads")
    resources:
        mem=get_resource("minimap2_align", "mem"),
        walltime=get_resource("minimap2_align", "walltime")
    log:
        f"{LOGDIR}/minimap2/{{sample}}.log"
    conda: "../envs/minimap.yaml"
    shell:
        """
        minimap2 --MD -ax map-ont -t {threads} \
         -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
         {input.genome} {input.fq}/*.fastq.gz | \
         samtools sort -@ {threads} -o {output} - 2> {log}
        """

rule minimap2_pbsv_align:
    input:
        fq = get_samples,
        genome = config["genome"]
    output:
        f"{OUTDIR}/minimap2_pbsv/alignment/{{sample}}.bam"
    threads: get_resource("minimap2_align", "threads")
    resources:
        mem=get_resource("minimap2_align", "mem"),
        walltime=get_resource("minimap2_align", "walltime")
    log:
        f"{OUTDIR}/minimap2_pbsv/{{sample}}.log"
    conda: "../envs/minimap.yaml"
    shell:
        """
        minimap2 -ax map-ont --MD --eqx -L -O 5,56 -E 4,1 -B 5 \
         --secondary=no -z 400,50 -r 2k -Y \
         -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
         -t {threads} {input.genome} {input.fq}/*.fastq.gz | \
         samtools sort -@ {threads} -o {output} - 2> {log}"""

rule ngmlr_align:
    input:
        fq = get_samples,
        genome = config["genome"]
    output:
        protected(f"{OUTDIR}/ngmlr/alignment/{{sample}}.bam")
    threads: get_resource("ngmlr_align", "threads")
    resources:
        mem=get_resource("ngmlr_align", "mem"),
        walltime=get_resource("ngmlr_align", "walltime")
    log:
        f"{LOGDIR}/ngmlr/{{sample}}.log"
    conda: "../envs/ngmlr.yaml"
    shell:
        "zcat {input.fq}/*.fastq.gz | \
         ngmlr --presets ont -t {threads} -r {input.genome} | \
         samtools sort -@ {threads} -o {output} - 2> {log}"

rule samtools_index:
    input:
        f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam"
    output:
        f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai"
    threads: get_resource("samtools_index", "threads")
    resources:
        mem=get_resource("samtools_index", "mem"),
        walltime=get_resource("samtools_index", "walltime")
    conda: "../envs/samtools.yaml"
    log:
        f"{OUTDIR}/{{aligner}}/samtools_index/{{sample}}.log"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"

rule alignment_stats:
    input:
        bam = [f"{OUTDIR}/{{aligner}}/alignment/{sample}.bam" for sample in config["samples"]],
        bai = [f"{OUTDIR}/{{aligner}}/alignment/{sample}.bam.bai" for sample in config["samples"]]
    output:
        f"{OUTDIR}/{{aligner}}/alignment_stats/{{sample}}.txt"
    threads: get_resource("alignment_stats", "threads")
    resources:
        mem=get_resource("alignment_stats", "mem"),
        walltime=get_resource("alignment_stats", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/alignment_stats/{{sample}}.log"
    conda: "../envs/pysam.yaml"
    shell:
        os.path.join(workflow.basedir, "scripts/alignment_stats.py") + \
            " -o {output} {input.bam} 2> {log}"


rule make_last_index:
    input:
        config["genome"]
    output:
        wmstat = f"{OUTDIR}/last/index/genome.wmstat",
        masked_genome = f"{OUTDIR}/last/index/genome-wm.fa",
        indexf = f"{OUTDIR}/last/index/windowmasked-index.bck",
    log:
        f"{LOGDIR}/last/mask_and_build_index/index.log"
    threads: get_resource("make_last_index", "threads")
    resources:
        mem=get_resource("make_last_index", "mem"),
        walltime=get_resource("make_last_index", "walltime")
    params:
        index_base = f"{OUTDIR}/last/index/windowmasked-index"
    conda: "../envs/last.yaml"
    shell:
        """
        bin/windowmasker -mk_counts -in {input} > {output.wmstat} && \
        bin/windowmasker -ustat {output.wmstat} -outfmt fasta -in {input} > {output.masked_genome} && \
        lastdb -P{threads} -uNEAR -R11 -c {params.index_base} {output.masked_genome}
        """

rule last_train:
    input:
        fq = get_all_samples,
        indexf = f"{OUTDIR}/last/index/windowmasked-index.bck",
    output:
        params = f"{OUTDIR}/last/index/last-train.params",
        fqs = temp(f"{OUTDIR}/last/index/fqs.fofn"),
        fas = temp(f"{OUTDIR}/last/index/allreads.fas"),
    log:
        f"{LOGDIR}/last/last-train/train.log"
    conda: "../envs/last.yaml"
    params:
        index_base = f"{OUTDIR}/last/index/windowmasked-index"
    threads: get_resource("last_train", "threads")
    resources:
        mem=get_resource("last_train", "mem"),
        walltime=get_resource("last_train", "walltime")
    shell:
        """
        for i in {input.fq}; do echo ${{i}}/*.fastq.gz >> {output.fqs}; done 2>> {log}
        zcat $(cat {output.fqs} | tr '\\n' ' ') \
         | awk 'NR % 4 == 2 {{print ">" ++n "\\n" $0}}' > {output.fas} 2>> {log}
        last-train -P{threads} -Q0 {params.index_base} {output.fas} > {output.params} 2>> {log}
        """

rule last_align:
    input:
        fq = get_samples,
        genome = config["genome"],
    threads: get_resource("last_align", "threads")
    resources:
        mem=get_resource("last_align", "mem"),
        walltime=get_resource("last_align", "walltime")
    params:
        index_base = config['last-index'],
        train = config['last-train'],
    log:
        f"{LOGDIR}/last/last-align/{{sample}}.log"
    output:
        f"{OUTDIR}/last/last-align/{{sample}}.maf.gz"
    conda: "../envs/last.yaml"
    shell:
        """
        lastal -P{threads} -p {params.train} {params.index_base} {input.fq}/*.fastq.gz \
         | last-split -m1e-6 \
         | gzip > {output} 2> {log}
        """
