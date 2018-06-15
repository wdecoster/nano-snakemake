from cyvcf2 import VCF
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from argparse import ArgumentParser
import sys


def get_args():
    parser = ArgumentParser(description="Comparison of vcf sets merged by SURVIVOR")
    parser.add_argument("--vcf", help="vcf file merged by SURVIVOR containing two or three sets")
    parser.add_argument("--names", help="names of sets", nargs='*')
    return parser.parse_args()


def is_variant(call):
    """Check if a variant position qualifies as a variant"""
    if call == 1 or call == 3:
        return True
    else:
        return False


def check_vcf(vcf):
    """Check if vcf is suited for this script

    returns
    - appropriate venn function
    - appropriate empty list of lists for identifiers
    """
    if len(vcf.samples) == 2:
        return venn2, [[] for i in range(len(vcf.samples))]
    elif len(vcf.samples) == 3:
        return venn3, [[] for i in range(len(vcf.samples))]
    else:
        sys.exit("Fatal: Script only written for vcf files containing 2 or 3 samples")


def main():
    args = get_args()
    vcf = VCF(args.vcf)
    venn, identifier_list = check_vcf(vcf)
    names = args.names or vcf.samples

    for v in vcf:
        for index, call in enumerate(v.gt_types):
            if is_variant(call):
                identifier_list[index].append(v.ID)
    venn([set(i) for i in identifier_list], set_labels=names)
    plt.savefig(args.vcf.replace('.vcf', '.png'))


if __name__ == '__main__':
    main()
