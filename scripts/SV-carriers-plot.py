#! /usr/bin/env python

from cyvcf2 import VCF
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(vcf, output):
    variants = VCF(vcf)
    plt.hist(x=[get_carrier_number(v.gt_types) for v in variants],
             bins=range(10),
             align="left")
    plt.xlabel('Number of carriers')
    plt.xticks(range(1, len(variants.samples) + 1), range(1, len(variants.samples) + 1))
    plt.ylabel('Number of variants')
    plt.savefig(output)


def get_carrier_number(gt_array):
    return sum(1 for i in gt_array if i in [1, 3])


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
