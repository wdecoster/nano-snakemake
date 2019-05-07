rule svim_call:
    input:
        "{aligner}/alignment/{sample}.bam"
    output:
        "{aligner}/svim_calls/{sample}/final_results.vcf"
    threads: 1
    log:
        "logs/{aligner}/svim_call/{sample}.log"
    shell:
        "svim alignment --sample {wildcards.sample} \
         {wildcards.aligner}/svim_calls/{wildcards.sample}/ {input} 2> {log}"

rule filter_svim:
    input:
        "{aligner}/svim_calls/{sample}/final_results.vcf"
    output:
        "{aligner}/svim_genotypes/{sample,[A-Za-z0-9]+}.vcf"
    threads: 1
    log:
        "logs/{aligner}/svim_call/{sample}.filter.log"
    shell:
        "cat {input} | \
         awk '{{ if($1 ~ /^#/) {{ print $0 }} \
         else {{ if($6>40) {{ print $0 }} }} }}' > {output}"

rule sniffles_call:
    input:
        "{aligner}/alignment/{sample}.bam"
    output:
        "{aligner}/sniffles_calls/{sample}.vcf"
    threads: 8
    log:
        "logs/{aligner}/sniffles_call/{sample}.log"
    shell:
        "sniffles --mapped_reads {input} --vcf {output} --threads {threads} 2> {log}"

rule sniffles_genotype:
    input:
        bam = "{aligner}/alignment/{sample}.bam",
        ivcf = "{aligner}/sniffles_combined/calls.vcf"
    output:
        "{aligner}/sniffles_genotypes_temp/{sample}.vcf"
    threads: 8
    log:
        "logs/{aligner}/sniffles_genotype/{sample}.log"
    shell:
        "sniffles --mapped_reads {input.bam} \
                  --vcf {output} \
                  --threads {threads} \
                  --report_seq \
                  --cluster \
                  --Ivcf {input.ivcf} 2> {log}"

rule samtools_split:
    input:
        bam = "{aligner}/alignment/{sample}.bam",
        bai = "{aligner}/alignment/{sample}.bam.bai",
    output:
        temp("{aligner}/alignment/{sample}-{chromosome}.bam")
    params:
        chrom = "{chromosome}"
    log:
        "logs/{aligner}/samtools_split/{sample}-{chromosome}.log"
    shell:
        "samtools view {input.bam} {params.chrom} -o {output} 2> {log}"


rule nanosv_call:
    '''

    call variants using NanoSV on separate chromosomes
    the shell command will first check if there are reads in this chromosome
    and if not, will just touch the output and leave it empty
    without raising an error
    '''
    input:
        bam = "{aligner}/alignment/{sample}-{chromosome}.bam",
        bai = "{aligner}/alignment/{sample}-{chromosome}.bam.bai",
        bed = config["annotbed"]
    output:
        temp("{aligner}/split_nanosv_genotypes/{sample}-{chromosome}.vcf")
    params:
        samtools = "samtools"
    threads:
        2
    log:
        "logs/{aligner}/nanosv/{sample}-{chromosome}.log"
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
                    --sambamba {params.samtools} {input.bam} \
                    -o {output} 2> {log}
        fi
        """

rule npinv:
    input:
        bam = "{aligner}/alignment/{sample}.bam",
        bai = "{aligner}/alignment/{sample}.bam.bai",
    output:
        "{aligner}/npinv/{sample}.vcf"
    log:
        "logs/{aligner}/npinv/{sample}.log"
    shell:
        "npinv --input {input.bam} --output {output} 2> {log}"

rule pbsv:
    input:
        bam = "minimap2_pbsv/alignment/{sample}.bam",
        bai = "minimap2_pbsv/alignment/{sample}.bam.bai",
        genome = config["genome"],
    output:
        vcf = "minimap2_pbsv/pbsv_genotypes/{sample}.vcf",
        svsig = temp("minimap2_pbsv/pbsv_svsig/{sample}.svsig.gz"),
    log:
        "logs/minimap2_pbsv/pbsv/{sample}.log"
    shell:
        """
        pbsv discover {input.bam} {output.svsig} 2> {log} && \
        pbsv call {input.genome} {output.svsig} {output.vcf} 2>> {log}
        """

rule tandem_genotypes:
    input:
        "last/last-align/{sample}.maf.gz"
    output:
        "last/tandem-genotypes/{sample}-tg.txt"
    params:
        microsat = config["microsat"],
        refgene = config["refgene"]
    log:
        "logs/tandem_genotypes/{sample}.log"
    shell:
        """
        tandem-genotypes -g {params.refgene} {params.microsat} {input} > {output} 2> {log}
        """

rule reformat_tandem_genotypes:
    input:
        expand("last/tandem-genotypes/{sample}-tg.txt", sample=config["samples"])
    output:
        "last/tandem_genotypes_reformatted/combined.txt"
    params:
        names = ' '.join(config["samples"].keys())
    log:
        'last/tandem_genotypes/reformat.log'
    shell:
        os.path.join(workflow.basedir, "scripts/reformat_tandem-genotypes.py") + \
            " --input {input} --names {params.names} > {output} 2> {log}"
