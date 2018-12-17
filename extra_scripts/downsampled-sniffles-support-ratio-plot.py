import pandas as pd
from argparse import ArgumentParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    args = get_args()
    df = pd.read_csv(filepath_or_buffer=args.pr,
                     sep="\t",
                     header=None,
                     names=['depth_ratio', 'support', 'precision', 'recall'])
    ratio_to_depth = {0.1: 5,  # hard coded conversion dict depending on your coverage
                      0.2: 11,
                      0.3: 17,
                      0.4: 23,
                      0.5: 29,
                      0.6: 35,
                      0.7: 41,
                      0.8: 47,
                      0.9: 52,
                      1:  59}
    df = df.assign(depth=df.depth_ratio.map(ratio_to_depth))
    plt.figure()
    ax = plt.gca()
    for feat in ["precision", "recall"]:
        ax.scatter(x=df["support"] / df["depth"],
                   y=df[feat],
                   label=feat,
                   s=3)
    plt.legend(loc="center left")
    plt.axvline(x=0.25)
    plt.axvline(x=0.33)
    plt.xlabel("Ratio support/coverage")
    plt.ylim(0, 100)
    plt.title("Sniffles minimal support")
    plt.tight_layout()
    plt.savefig("sniffles_support-vs-pr.png")


def get_args():
    parser = ArgumentParser(
        description="Plot precision and recall with different sniffles support and coverage")
    parser.add_argument("pr", help="file with 4 columns, of which 3 and 4 are precision and recall")
    return parser.parse_args()


if __name__ == '__main__':
    main()
