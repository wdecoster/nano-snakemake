#! /usr/bin/env python

from cyvcf2 import VCF
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(vcf, output):
    plt.hist(x=[get_carrier_number(v.gt_types) for v in VCF(vcf)],
             bins=range(10),
             align="left")
    plt.xlabel('Number of carriers')
    plt.ylabel('Number of variants')
    plt.savefig(output)


def get_carrier_number(gt_array):
    return sum(gt_array > 0)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
