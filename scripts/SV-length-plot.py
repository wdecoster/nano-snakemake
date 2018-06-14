#! /usr/bin/env python

from cyvcf2 import VCF
import sys
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(vcf_file):
    len_dict = defaultdict(list)
    for v in VCF(vcf_file):
        if not v.INFO.get('SVTYPE') == 'TRA':
            len_dict[v.INFO.get('SVTYPE')].append(v.INFO.get('SVLEN'))
    plt.hist(x=np.array(list(len_dict.values())),
             bins=[i for i in range(0, 2000, 50)],
             stacked=True,
             histtype='bar',
             label=list(len_dict.keys()))
    plt.xlabel('Lenghth of structural variant')
    plt.ylabel('Number of variants')
    plt.legend()
    plt.savefig(sys.argv[2])


if __name__ == '__main__':
    main(sys.argv[1])
