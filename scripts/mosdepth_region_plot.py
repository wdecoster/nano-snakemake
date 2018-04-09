#! /usr/bin/env python

import pandas as pd
from argparse import ArgumentParser
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt


def main():
    args = get_args()
    region_df = read_mosdepth(args.mosdepth_combined, args.region)
    fig, ax = plt.subplots(figsize=(16, 8))
    ax = sns.pointplot(x="begin",
                       y="coverage",
                       hue="name",
                       data=region_df,
                       scale=0.1,
                       ax=ax)
    ax.set(xlabel="position",
           ylabel="normalizes coverage")
    plt.xticks([tick for i, tick in enumerate(list(plt.xticks()[0])) if not i % 5],
               rotation=30,
               ha='center')
    plt.savefig("Mosdepth_{}.png".format(args.region), dpi=500, bbox_inches='tight')


def get_args():
    parser = ArgumentParser(description="create a coverage plot across a locus using mosdepth data")
    parser.add_argument("-r", "--region", help="Region of interest", required=True)
    parser.add_argument("mosdepth_combined", help="combined mosdepth file")
    return parser.parse_args()


def read_mosdepth(mosdepth_file, region):
    reg_chrom, reg_begin, reg_end = region.replace('-', ':').split(':')
    df = normalize_reads(pd.read_csv(mosdepth_file, sep="\t"))
    fdf = df.loc[(df["chr"] == reg_chrom) &
                 (df["begin"] > int(reg_begin)) &
                 (df["end"] < int(reg_end))]
    return pd.melt(fdf,
                   id_vars=['chr', 'begin', 'end'],
                   value_vars=[i for i in fdf.columns if i not in ["chr", "begin", "end"]],
                   var_name="name",
                   value_name="coverage")


def normalize_reads(df):
    samples = [i for i in df.columns if i not in ["chr", "begin", "end"]]
    read_sums = df.loc[:, samples].sum()
    for sample, read_sum in zip(samples, read_sums):
        df.loc[:, sample] = df.loc[:, sample] / read_sum * 1e6
    return df


if __name__ == '__main__':
    main()
