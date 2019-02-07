import subprocess
from argparse import ArgumentParser
import shlex
import sys
import os
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
    combined_vcf = survivor(samples=[normalize_vcf(s) for s in [args.truth, args.test]],
                            distance=args.distance,
                            ignore_type=ignore_type,
                            minlength=args.minlength)
    truth_set, test_set = make_venn(combined_vcf)
    tp = len(truth_set & test_set)
    precision = tp / len(test_set)
    recall = tp / len(truth_set)
    fmeasure = 2*(precision*recall)/(precision + recall)
    print(f"Precision: {round(precision, ndigits=4)}%")
    print(f"Recall: {round(recall, ndigits=4)}%")
    print(f"F-measure: {round(fmeasure, ndigits=4)}")
    if args.bar:
        bar_chart(combined_vcf)


def normalize_vcf(vcff):
    handle, name = tempfile.mkstemp()
    out = open(name, 'w')
    for line in open(vcff):
        out.write(line.replace('DUP', 'INS'))
    os.close(handle)
    return name


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
    fhf, fofn_f = tempfile.mkstemp()
    fhv, vcf_out = tempfile.mkstemp()
    with open(fofn_f, 'w') as fofn:
        for s in samples:
            fofn.write(s + "\n")
    survivor_cmd = f"SURVIVOR merge {fofn_f} {distance} 1 {ignore_type} -1 -1 {minlength} {vcf_out}"
    sys.stderr.write("Executing SURVIVOR...\n")
    subprocess.call(shlex.split(survivor_cmd), stdout=subprocess.DEVNULL)
    os.close(fhf)
    os.close(fhv)
    return vcf_out


def is_variant(call):
    """Check if a variant position qualifies as a variant

    0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT"""
    if call == 1 or call == 3:
        return True
    else:
        return False


def make_venn(vcf, outname="venn.png"):
    positions = [[], []]
    for v in VCF(vcf):
        for index, call in enumerate(v.gt_types):
            if is_variant(call):
                positions[index].append("{}:{}-{}".format(v.CHROM, v.start, v.INFO.get('SVTYPE')))
    identifier_sets = [set(i) for i in positions]
    venn2(identifier_sets, set_labels=('Truth', 'Test'))
    plt.savefig(outname)
    plt.close()
    return identifier_sets


def bar_chart(vcf, outname="stacked_bar.png"):
    """
    Make a stacked bar chart for length of the SV split by validation status
    This ignores zygosity.
    """
    len_dict = {"True": [], "False": [], "Missed": []}
    for v in VCF(vcf):
        if not v.INFO.get('SVTYPE') == 'TRA' and abs(v.INFO.get('AVGLEN')) >= 50:
            calls = [is_variant(call) for call in v.gt_types]
            if calls == [True, True]:
                len_dict['True'].append(v.INFO.get('AVGLEN'))
            elif calls == [False, True]:
                len_dict['False'].append(v.INFO.get('AVGLEN'))
            elif calls == [True, False]:
                len_dict['Missed'].append(v.INFO.get('AVGLEN'))
    plt.subplot(2, 1, 1)
    plt.hist(x=np.array(list(len_dict.values())),
             bins=[i for i in range(0, 2000, 10)],
             stacked=True,
             histtype='bar',
             label=list(len_dict.keys()))
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")
    plt.subplot(2, 1, 2)
    plt.hist(x=np.array(list(len_dict.values())),
             bins=[i for i in range(0, 20000, 100)],
             stacked=True,
             histtype='bar',
             label=list(len_dict.keys()),
             log=True)
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")
    plt.tight_layout()
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
                        default=50)
    parser.add_argument("-i", "--ignore_type",
                        help="Ignore the type of the structural variant",
                        action="store_true")
    parser.add_argument("--bar",
                        help="Make stacked bar chart of SV lengths coloured by validation status",
                        action="store_true")
    return parser.parse_args()


if __name__ == '__main__':
    main()
