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
        "npinv --input {input} --output {output}"

rule pbsv:
    input:
        bam = "minimap2_pbsv/alignment/{sample}.bam",
        bai = "minimap2_pbsv/alignment/{sample}.bam.bai",
        genome = genome = config["genome"]
    output:
        vcf = "minimap2_pbsv/pbsv/{sample}.vcf",
        svsig = temp("minimap2_pbsv/pbsv/{sample}.svsig.gz)
    log:
        "logs/minimap2_pbsv/pbsv/{sample}.log"
    shell:
        """
        pbsv discover {input.bam} {output.svsig} && \
        pbsv call {input.genome} {output.svsig} {output.vcf}
        """
