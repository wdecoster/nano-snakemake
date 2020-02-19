configfile: "config.yaml"

rule mosdepth_get:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai"
    threads: get_resource("mosdepth_get", "threads")
    resources:
        mem=get_resource("mosdepth_get", "mem"),
        walltime=get_resource("mosdepth_get", "walltime")
    output:
        protected(f"{OUTDIR}/{{aligner}}/mosdepth/{{sample}}.mosdepth.global.dist.txt"),
        protected(f"{OUTDIR}/{{aligner}}/mosdepth/{{sample}}.regions.bed.gz"),
    params:
        windowsize = 500,
        outdir = f"{OUTDIR}/{{aligner}}/mosdepth/{{sample}}"
    log:
        f"{LOGDIR}/{{aligner}}/mosdepth/mosdepth_{{sample}}.log"
    conda: "../envs/mosdepth.yaml"
    shell:
        "mosdepth --threads {threads} \
                  -n \
                  --by {params.windowsize} \
                  {params.outdir} {input.bam} 2> {log}"


rule mosdepth_combine:
    input:
        [f"{OUTDIR}/{{aligner}}/mosdepth/{sample}.regions.bed.gz" for sample in config["samples"]]
    output:
        f"{OUTDIR}/{{aligner}}/mosdepth/regions.combined.gz"
    threads: get_resource("mosdepth_combine", "threads")
    resources:
        mem=get_resource("mosdepth_combine", "mem"),
        walltime=get_resource("mosdepth_combine", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/mosdepth/mosdepth_combine.log"
    shell:
        os.path.join(workflow.basedir, "scripts/combine_mosdepth.py") + \
            " {input} -o {output} 2> {log}"


rule mosdepth_global_plot:
    input:
        [f"{OUTDIR}/{{aligner}}/mosdepth/{sample}.mosdepth.global.dist.txt" for sample in config["samples"]]
    output:
        f"{OUTDIR}/{{aligner}}/mosdepth_global_plot/global.html"
    threads: get_resource("mosdepth_global_plot", "threads")
    resources:
        mem=get_resource("mosdepth_global_plot", "mem"),
        walltime=get_resource("mosdepth_global_plot", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/mosdepth/mosdepth_global_plot.log"
    shell:
        os.path.join(workflow.basedir, "scripts/mosdepth_plot-dist.py") + \
            " {input} -o {output} 2> {log}"
