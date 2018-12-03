#! /usr/bin/env python
import sys
from cyvcf2 import VCF
from collections import defaultdict
from argparse import ArgumentParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    args = get_args()
    len_dict = defaultdict(list)
    for v in VCF(args.vcf):
        if not v.INFO.get('SVTYPE') == 'TRA':
            try:
                if abs(v.INFO.get('SVLEN')) >= 50:
                    len_dict[v.INFO.get('SVTYPE')].append(abs(v.INFO.get('SVLEN')))
            except TypeError:
                sys.stderr.write("Exception when parsing variant:\n{}\n\n".format(v))
    with open(args.counts, 'w') as counts:
        counts.write("Number of nucleotides affected by SV:\n")
        for svtype, lengths in len_dict.items():
            counts.write("{}:\t{} variants\t{}bp\n".format(
                svtype, len(lengths), sum(lengths)))
    make_plot(dict_of_lengths=len_dict,
              output=args.output)


def make_plot(dict_of_lengths, output):
    """Makes two stacked bar charts
    Plotting two bar charts of number of SVs by length split by SV type
    Use a consistent colouring scheme for those in "standard_order" to
    make comparison reasonable

    First bar chart is up to 2kb with bins of 10bp
    Second bar chart is up to 20kb, with bins of 100bp
     and uses log scaling on the y-axis
    """
    standard_order = ['DEL', 'INS', 'DUP', 'INV']
    spec_order = sorted([i for i in dict_of_lengths.keys() if i not in standard_order])
    sorter = standard_order + spec_order
    names, lengths = zip(
        *sorted([(svtype, lengths) for svtype, lengths in dict_of_lengths.items()],
                key=lambda x: sorter.index(x[0])))
    plt.subplot(2, 1, 1)
    plt.hist(x=lengths,
             bins=[i for i in range(0, 2000, 10)],
             stacked=True,
             histtype='bar',
             label=names)
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")

    plt.subplot(2, 1, 2)
    plt.hist(x=lengths,
             bins=[i for i in range(0, 20000, 100)],
             stacked=True,
             histtype='bar',
             label=names,
             log=True)
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")
    plt.tight_layout()
    plt.savefig(output)


def get_args():
    parser = ArgumentParser(description="create stacked bar plot of the SV lengths split by type")
    parser.add_argument("vcf", help="vcf file to parse")
    parser.add_argument("-o", "--output",
                        help="output file to write figure to",
                        default="SV-length.png")
    parser.add_argument("-c", "--counts",
                        help="output file to write counts to",
                        default="SV-length.txt")
    return parser.parse_args()


if __name__ == '__main__':
    main()
