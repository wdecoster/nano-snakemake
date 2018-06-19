import subprocess
from argparse import ArgumentParser
import shlex
import sys
from cyvcf2 import VCF
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


def main():
    args = get_args()
    combined_vcf = survivor(args)
    truth_set, test_set = make_venn(combined_vcf)
    tp = len(truth_set & test_set)
    print(f"Precision: {round(100 * tp / len(test_set), ndigits=1)}%")
    print(f"Recall: {round(100 * tp / len(truth_set), ndigits=1)}%")


def survivor(args):
    """
    Executes SURVIVOR merge, with parameters:
    -samples.fofn (truth and test)
    -distance between calls (args.distance)
    -number of callers to support call (1)
    -require variants to have sampe type (args.ignore_type)
    -require variants to be on same strand (no)
    -estimate distance between calls (no)
    -specify minimal size of SV event (no)
    """
    samples = "samples.fofn"
    vcf_out = "combined.vcf"
    with open(samples, 'w') as fofn:
        fofn.write(args.truth + "\n")
        fofn.write(args.test + "\n")
    if args.ignore_type:
        ignore_type = "-1"
    else:
        ignore_type = "1"
    survivor_cmd = f"SURVIVOR merge {samples} {args.distance} 1 {ignore_type} -1 -1 -1 {vcf_out}"
    sys.stderr.write("Executing SURVIVOR...\n")
    subprocess.call(shlex.split(survivor_cmd))
    return VCF(vcf_out)


def is_variant(call):
    """Check if a variant position qualifies as a variant"""
    if call == 1 or call == 3:
        return True
    else:
        return False


def make_venn(vcf, outname="venn.png"):
    identifier_list = [[], []]
    for v in vcf:
        for index, call in enumerate(v.gt_types):
            if is_variant(call):
                identifier_list[index].append(v.ID)
    identifier_sets = [set(i) for i in identifier_list]
    venn2(identifier_sets, set_labels=('Truth', 'Test'))
    plt.savefig(outname)
    return identifier_sets


def get_args():
    parser = ArgumentParser(description="Calculate precision-recall metrics from 2 vcf files")
    parser.add_argument("--truth", help="vcf containing truth set", required=True)
    parser.add_argument("--test", help="vcf containing test set", required=True)
    parser.add_argument("-d", "--distance",
                        help="maximum distance between test and truth call",
                        default=500)
    parser.add_argument("-i", "--ignore_type",
                        help="Ignore the type of the structural variant",
                        action="store_true")
    return parser.parse_args()


if __name__ == '__main__':
    main()
