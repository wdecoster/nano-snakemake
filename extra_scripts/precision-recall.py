import subprocess
from argparse import ArgumentParser
import shlex
import sys
from cyvcf2 import VCF
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import tempfile
import numpy as np


def main():
    args = get_args()
    if args.ignore_type:
        ignore_type = "-1"
    else:
        ignore_type = "1"
    combined_vcf = survivor(samples=[args.truth, args.test],
                            distance=args.distance,
                            ignore_type=ignore_type,
                            minlength=args.minlength)
    truth_set, test_set = make_venn(combined_vcf)
    tp = len(truth_set & test_set)
    print(f"Precision: {round(100 * tp / len(test_set), ndigits=1)}%")
    print(f"Recall: {round(100 * tp / len(truth_set), ndigits=1)}%")
    if args.bar:
        bar_chart(combined_vcf)


def survivor(samples, distance, ignore_type, minlength):
    """
    Executes SURVIVOR merge, with parameters:
    -samples.fofn (truth and test)
    -distance between calls (args.distance)
    -number of callers to support call (1)
    -require variants to have sampe type (args.ignore_type)
    -require variants to be on same strand (no)
    -estimate distance between calls (no)
    -specify minimal size of SV event (args.minlength)
    """
    fofn_f = tempfile.mkstemp()[1]
    vcf_out = tempfile.mkstemp()[1]
    with open(fofn_f, 'w') as fofn:
        for s in samples:
            fofn.write(s + "\n")
    survivor_cmd = f"SURVIVOR merge {fofn_f} {distance} 1 {ignore_type} -1 -1 {minlength} {vcf_out}"
    sys.stderr.write("Executing SURVIVOR...\n")
    subprocess.call(shlex.split(survivor_cmd))
    return vcf_out


def is_variant(call):
    """Check if a variant position qualifies as a variant

    0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT"""
    if call == 1 or call == 3:
        return True
    else:
        return False


def make_venn(vcf, outname="venn.png"):
    identifier_list = [[], []]
    for v in VCF(vcf):
        for index, call in enumerate(v.gt_types):
            if is_variant(call):
                identifier_list[index].append(v.ID)
    identifier_sets = [set(i) for i in identifier_list]
    venn2(identifier_sets, set_labels=('Truth', 'Test'))
    plt.savefig(outname)
    plt.close()
    return identifier_sets


def bar_chart(vcf, outname="stacked_bar.png"):
    """
    Make a stacked bar chart for length of the SV split by validation status
    This ignores zygosity.
    """
    len_dict = {"true": [], "false": [], "missed": []}
    for v in VCF(vcf):
        if not v.INFO.get('SVTYPE') == 'TRA':
            calls = [is_variant(call) for call in v.gt_types]
            if calls == [True, True]:
                len_dict['true'].append(v.INFO.get('AVGLEN'))
            elif calls == [False, True]:
                len_dict['false'].append(v.INFO.get('AVGLEN'))
            elif calls == [True, False]:
                len_dict['missed'].append(v.INFO.get('AVGLEN'))
    plt.hist(x=np.array(list(len_dict.values())),
             bins=[i for i in range(0, 2000, 50)],
             stacked=True,
             histtype='bar',
             label=list(len_dict.keys()))
    plt.xlabel('Lenghth of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")
    plt.savefig(outname)
    plt.close()


def get_args():
    parser = ArgumentParser(description="Calculate precision-recall metrics from 2 vcf files")
    parser.add_argument("--truth", help="vcf containing truth set", required=True)
    parser.add_argument("--test", help="vcf containing test set", required=True)
    parser.add_argument("-d", "--distance",
                        help="maximum distance between test and truth call",
                        default=500)
    parser.add_argument("--minlength",
                        help="Minimum length of SVs to be taken into account",
                        default=-1)
    parser.add_argument("-i", "--ignore_type",
                        help="Ignore the type of the structural variant",
                        action="store_true")
    parser.add_argument("--bar",
                        help="Make stacked bar chart of SV lengths coloured by validation status",
                        action="store_true")
    return parser.parse_args()


if __name__ == '__main__':
    main()
