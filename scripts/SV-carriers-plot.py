#! /usr/bin/env python

import matplotlib.pyplot as plt
from cyvcf2 import VCF
from argparse import ArgumentParser
import matplotlib
matplotlib.use('Agg')


def main():
    args = get_args()
    variants = VCF(args.vcf)
    plt.hist(x=[get_carrier_number(v.gt_types) for v in variants],
             bins=range(10),
             align="left")
    plt.xlabel('Number of carriers')
    plt.xticks(range(1, len(variants.samples) + 1), range(1, len(variants.samples) + 1))
    plt.ylabel('Number of variants')
    plt.savefig(args.output)


def get_carrier_number(gt_array):
    return sum(1 for i in gt_array if i in [1, 3])


def get_args():
    parser = ArgumentParser(description="create bar plot of the number of carriers per SV")
    parser.add_argument("vcf", help="vcf file to parse")
    parser.add_argument("-o", "--output",
                        help="output file to write figure to",
                        default="SV-carriers.png")
    return parser.parse_args()


if __name__ == '__main__':
    main()
