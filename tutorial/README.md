# SV-snakemake
A snakemake pipeline for SV analysis from nanopore genome sequencing

## Getting my feet wet with snakemake

- rules in a snakefile: steps in the workflow
- dependencies between steps automatically determined by matching file names

```python
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

input and output specify lists of file names (note the comma)

{input} contains multiple values, which will be concatenated and separated by whitespace  
missing output directories are automatically created  


When executing a workflow Snakemake tries to generate the given target files, which is called a job.  
e.g. `snakemake mapped_reads/A.bam`

extra arguments:  
-n/--dryrun    show execution plan
-p             print the shell command for illustration

A job will not be reran if the output already exists, except if one of the input files is newer than one of the output files (or will be updated by another job).

rules can be generalized using *named wildcards*

```python
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

Specifying input files will fill in the wildcards  
`snakemake -np mapped_reads/A.bam mapped_reads/B.bam`
alternatively, bash magic:  
`snakemake -np mapped_reads/{A,B}.bam`


```python
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
```

Accessing the value of the wildcard in the shell command requires the {wildcards.sample} syntax.

## Visualizing the directed acyclic graph (DAG)
`snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg`
this requires graphviz to be installed (conda)


## Collecting input files using helper function expand
Create a list of files in which {sample} if replaced by all samples in the list SAMPLES
```python
SAMPLES = ["A", "B"]
expand("sorted_reads/{sample}.bam", sample=SAMPLES)
```


## Specifying names for input/output files
```python
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
```

With multiple input or output files you can use names to specify those  
Note that those those file lists can be python expressions, as long as they return a string or list of strings

## Creating a report
```python
rule report:
    input:
        "calls/all.vcf"
    output:
        "report.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as vcf:
            n_calls = sum(1 for l in vcf if not l.startswith("#"))

        report("""
        An example variant calling workflow
        ===================================

        Reads were mapped to the Yeast
        reference genome and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], T1=input[0])
```

Using the built-in report function, which takes in an rst markup string, the path where the report will be stored (output).  
This rule uses a `run` rather than `shell` directive, which is just python code.

## Adding a target rule
It's also possible to specify a rule name as target, if thise one does not contain any wildcards.  
If no target is given on the command line, the first rule of the Snakefile is the target. Therefore it's best practice to hava `all` rule at the top of the workflow as below:
```python
rule all:
    input:
        "report.html"
```

The order of the other rules does not matter and does not influence the DAG.

## specifying the number of threads for a rule

```python
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

When executing a workflow the scheduler makes sure not more than the available CPU cores are used, or the explicitly set maximum using --cores.  
When less --cores are specified than threads, the threads directive will be reduced to the number of given cores.

## Using a config file
Snakemake can use a YAML or JSON config file to store variables
```python
configfile: "config.yaml"
```

```YAML
samples:
    A: data/samples/A.fastq
    B: data/samples/B.fastq
```

the SAMPLES list can be replaced by the information from the yaml config:

```
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
```

## using an input function
Rather than a string as input, also a function can be specified

```python
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

## Setting arbitrary parameters

Can use the `params` directive to add arbitrary parameters to shell commands

```python
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mapped_reads/{sample}.bam"
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    threads: 8
    shell:
        "bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output}"
```

You can read such a parameter from the config file.

## logging

Can specify a `log` directive containing a file in which logs have to be generated
```python
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mapped_reads/{sample}.bam"
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
```


## Using temporary and protected files
Files can be marked as temporary, to make sure those are cleaned up after being used

```python
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
```

Similarly, files can be protected against accidental modification (setting file system permissions):
```python
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.bam")
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
```
