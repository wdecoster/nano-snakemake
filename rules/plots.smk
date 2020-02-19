rule SV_length_plot:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_{{stage}}/{{sample}}.vcf"
    output:
        plot = f"{OUTDIR}/{{aligner}}/SV-plots/SV-length_{{caller}}_{{stage}}_{{sample}}.png",
        counts = f"{OUTDIR}/{{aligner}}/SV-plots/SV-nucleotides_affected_{{caller}}_{{stage}}_{{sample}}.txt",
    threads: get_resource("SV_length_plot", "threads")
    resources:
        mem=get_resource("SV_length_plot", "mem"),
        walltime=get_resource("SV_length_plot", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/svplot/svlength_{{caller}}_{{stage}}_{{sample}}.log"
    conda: "../envs/plots.yaml"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-length-plot.py") + \
            " {input} --output {output.plot} --counts {output.counts} 2> {log}"


rule SV_plot_carriers:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_combined/annot_genotypes.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/SV-plots/SV-{{caller}}_carriers.png"
    threads: get_resource("SV_plot_carriers", "threads")
    resources:
        mem=get_resource("SV_plot_carriers", "mem"),
        walltime=get_resource("SV_plot_carriers", "walltime")
    log:
        f"{LOGDIR}/{{aligner}}/svplot/svcarriers_{{caller}}.log"
    conda: "../envs/plots.yaml"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-carriers-plot.py") + \
            " {input} --output {output} 2> {log}"
