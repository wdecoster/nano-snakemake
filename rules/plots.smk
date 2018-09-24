rule SV_length_plot:
    input:
        "{aligner}/{caller}_{stage}/{sample}.vcf"
    output:
        plot = "{aligner}/SV-plots/SV-length_{caller}_{stage}_{sample}.png",
        counts = "{aligner}/SV-plots/SV-nucleotides_affected_{caller}_{stage}_{sample}.txt",
    log:
        "logs/{aligner}/svplot/svlength_{caller}_{stage}_{sample}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-length-plot.py") + \
            " {input} --output {output.plot} --counts {output.counts} 2> {log}"


rule SV_plot_carriers:
    input:
        "{aligner}/{caller}_combined/annot_genotypes.vcf"
    output:
        "{aligner}/SV-plots/SV-{caller}_carriers.png"
    log:
        "logs/{aligner}/svplot/svcarriers_{caller}.log"
    shell:
        os.path.join(workflow.basedir, "scripts/SV-carriers-plot.py") + \
            " {input} {output} 2> {log}"
