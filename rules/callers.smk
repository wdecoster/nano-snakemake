rule svim_call:
    input:
        f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam"
    output:
        f"{OUTDIR}/{{aligner}}/svim_calls/{{sample}}/final_results.vcf"
    threads: get_resource("svim_call", "threads")
    resources:
        mem=get_resource("svim_call", "mem"),
        walltime=get_resource("svim_call", "walltime")
    conda: "../envs/svim.yaml"
    params:
        outdir=f"{OUTDIR}/{{aligner}}/svim_calls/{{sample}}"
    log:
        f"{OUTDIR}/{{aligner}}/svim_call/{{sample}}.log"
    shell:
        "svim alignment --sample {wildcards.sample} \
         {params.outdir}/ {input} 2> {log}"

rule filter_svim:
    input:
        f"{OUTDIR}/{{aligner}}/svim_calls/{{sample}}/final_results.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/svim_genotypes/{{sample,[A-Za-z0-9]+}}.vcf"
    threads: get_resource("filter_svim", "threads")
    resources:
        mem=get_resource("filter_svim", "mem"),
        walltime=get_resource("filter_svim", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/svim_call/{{sample}}.filter.log"
    shell:
        "cat {input} | \
         awk '{{ if($1 ~ /^#/) {{ print $0 }} \
         else {{ if($6>40) {{ print $0 }} }} }}' > {output}"

rule sniffles_call:
    input:
        f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam"
    output:
        f"{OUTDIR}/{{aligner}}/sniffles_calls/{{sample}}.vcf"
    threads: get_resource("sniffles_call", "threads")
    resources:
        mem=get_resource("sniffles_call", "mem"),
        walltime=get_resource("sniffles_call", "walltime")
    conda: "../envs/sniffles.yaml"
    log:
        f"{LOGDIR}/{{aligner}}/sniffles_call/{{sample}}.log"
    shell:
        "sniffles --mapped_reads {input} --vcf {output} --threads {threads} 2> {log}"

rule sniffles_genotype:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam",
        ivcf = f"{OUTDIR}/{{aligner}}/sniffles_combined/calls.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/sniffles_genotypes_temp/{{sample}}.vcf"
    threads: get_resource("sniffles_genotype", "threads")
    resources:
        mem=get_resource("sniffles_genotype", "mem"),
        walltime=get_resource("sniffles_genotype", "walltime")
    conda: "../envs/sniffles.yaml"
    log:
        f"{LOGDIR}/{{aligner}}/sniffles_genotype/{{sample}}.log"
    shell:
        "sniffles --mapped_reads {input.bam} \
                  --vcf {output} \
                  --threads {threads} \
                  --report_seq \
                  --cluster \
                  --Ivcf {input.ivcf} 2> {log}"

rule samtools_split:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai",
    output:
        temp(f"{OUTDIR}/{{aligner}}/alignment/{{sample}}-{{chromosome}}.bam")
    threads: get_resource("samtools_split", "threads")
    resources:
        mem=get_resource("samtools_split", "mem"),
        walltime=get_resource("samtools_split", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/samtools_split/{{sample}}-{{chromosome}}.log"
    conda: "../envs/samtools.yaml"
    shell:
        "samtools view {input.bam} {wildcards.chromosome} -o {output} 2> {log}"


rule nanosv_call:
    '''

    call variants using NanoSV on separate chromosomes
    the shell command will first check if there are reads in this chromosome
    and if not, will just touch the output and leave it empty
    without raising an error
    '''
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}-{{chromosome}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}-{{chromosome}}.bam.bai",
        bed = config["annotbed"]
    output:
        temp(f"{OUTDIR}/{{aligner}}/split_nanosv_genotypes/{{sample}}-{{chromosome}}.vcf")
    threads: get_resource("nanosv_call", "threads")
    resources:
        mem=get_resource("nanosv_call", "mem"),
        walltime=get_resource("nanosv_call", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/nanosv/{{sample}}-{{chromosome}}.log"
    conda: "../envs/nanosv.yaml"
    shell:
        """
        reads=$(samtools idxstats {input.bam} | \
          awk 'BEGIN {{FS = "\\t"}} ; {{sum+=$3}} END {{print sum}}')
        if [ "$reads" -eq "0" ]; then
            echo "##fileformat=VCFv4.1" > {output} && \
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> {output}
            echo "NanoSV: No reads in {input.bam}" >> exceptions.txt 2> {log}
        else
            NanoSV --bed {input.bed} \
                    --threads {threads} \
                    --sambamba samtools {input.bam} \
                    -o {output} 2> {log}
        fi
        """

rule npinv:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai",
    output:
        f"{OUTDIR}/{{aligner}}/npinv/{{sample}}.vcf"
    threads: get_resource("npinv", "threads")
    resources:
        mem=get_resource("npinv", "mem"),
        walltime=get_resource("npinv", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/npinv/{{sample}}.log"
    conda: "../envs/npinv.yaml"
    shell:
        "npinv --input {input.bam} --output {output} 2> {log}"

rule pbsv:
    input:
        bam = f"{OUTDIR}/minimap2_pbsv/alignment/{{sample}}.bam",
        bai = f"{OUTDIR}/minimap2_pbsv/alignment/{{sample}}.bam.bai",
        genome = config["genome"],
    output:
        vcf = f"{OUTDIR}/minimap2_pbsv/pbsv_genotypes/{{sample}}.vcf",
        svsig = temp(f"{OUTDIR}/minimap2_pbsv/pbsv_svsig/{{sample}}.svsig.gz"),
    threads: get_resource("pbsv", "threads")
    resources:
        mem=get_resource("pbsv", "mem"),
        walltime=get_resource("pbsv", "walltime")
    log:
        f"{LOGDIR}/minimap2_pbsv/pbsv/{{sample}}.log"
    conda: "../envs/pbsv.yaml"
    shell:
        """
        pbsv discover {input.bam} {output.svsig} 2> {log} && \
        pbsv call {input.genome} {output.svsig} {output.vcf} 2>> {log}
        """

rule tandem_genotypes:
    input:
        f"{OUTDIR}/last/last-align/{{sample}}.maf.gz"
    output:
        f"{OUTDIR}/last/tandem-genotypes/{{sample}}-tg.txt"
    threads: get_resource("tandem_genotypes", "threads")
    resources:
        mem=get_resource("tandem_genotypes", "mem"),
        walltime=get_resource("tandem_genotypes", "walltime")
    params:
        microsat = config["microsat"],
        refgene = config["refgene"]
    log:
        f"{LOGDIR}/tandem_genotypes/{{sample}}.log"
    conda: "../envs/tandem-genotypes.yaml"
    shell:
        """
        tandem-genotypes -g {params.refgene} {params.microsat} {input} > {output} 2> {log}
        """

rule reformat_tandem_genotypes:
    input:
        [f"{OUTDIR}/last/tandem-genotypes/{sample}-tg.txt" for sample in config["samples"]]
    output:
        f"{OUTDIR}/last/tandem_genotypes_reformatted/combined.txt"
    threads: get_resource("reformat_tandem_genotypes", "threads")
    resources:
        mem=get_resource("reformat_tandem_genotypes", "mem"),
        walltime=get_resource("reformat_tandem_genotypes", "walltime")
    params:
        names = ' '.join(config["samples"].keys())
    log:
        f"{LOGDIR}/tandem_genotypes/reformat.log"
    shell:
        os.path.join(workflow.basedir, "scripts/reformat_tandem-genotypes.py") + \
            " --input {input} --names {params.names} > {output} 2> {log}"
