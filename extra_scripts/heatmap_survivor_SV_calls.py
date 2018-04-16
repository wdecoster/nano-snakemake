from cyvcf2 import VCF
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def is_variant(call):
    if call == 1 or call == 3:
        return True
    else:
        return False


def main(vcf):
    variants = VCF(vcf)
    samples = variants.samples
    identifiers = {i: set() for i in samples}
    for v in variants:
        for sample, call in zip(samples, v.gt_types):
            if is_variant(call):
                identifiers[sample].add(v.ID)

    df = pd.DataFrame(index=samples, columns=samples)
    for i_sample, i_values in identifiers.items():
        for j_sample, j_values in identifiers.items():
            df.loc[i_sample, j_sample] = len(i_values & j_values)
    try:
        proper_names = [i.split('/')[1].split('.')[0] for i in samples]
    except IndexError:
        proper_names = samples
    sns.heatmap(data=df,
                annot=True,
                fmt="d",
                linewidths=0.5,
                xticklabels=proper_names,
                yticklabels=proper_names)
    plt.savefig("SV-calls_heatmap.png", bbox_inches="tight")


if __name__ == '__main__':
    main(sys.argv[1])
