import subprocess
from argparse import ArgumentParser
import shlex
import sys
import os
from cyvcf2 import VCF
import tempfile
import pandas as pd
from collections import defaultdict
from upsetplot import plot
import matplotlib.pyplot as plt


def main():
    args = get_args()
    if args.ignore_type:
        ignore_type = "-1"
    else:
        ignore_type = "1"
    combined_vcf = survivor(samples=[normalize_vcf(s) for s in args.variants],
                            distance=args.distance,
                            ignore_type=ignore_type,
                            minlength=args.minlength,
                            save=args.store)
    upsets = make_sets(vcf=combined_vcf,
                       names=args.names or args.variants)
    plot(upsets, sort_by='cardinality')
    plt.savefig("UpSetPlot.png")


def normalize_vcf(vcff):
    handle, name = tempfile.mkstemp()
    out = open(name, 'w')
    for line in open(vcff):
        out.write(line.replace('DUP', 'INS'))
    os.close(handle)
    return name


def survivor(samples, distance, ignore_type, minlength, save=False):
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
    if save:
        vcf_out = "overlapping-variants.vcf"
    else:
        fhv, vcf_out = tempfile.mkstemp()
    with open(fofn_f, 'w') as fofn:
        for s in samples:
            fofn.write(s + "\n")
    survivor_cmd = f"SURVIVOR merge {fofn_f} {distance} 1 {ignore_type} -1 -1 {minlength} {vcf_out}"
    sys.stderr.write("Executing SURVIVOR...\n")
    subprocess.call(shlex.split(survivor_cmd), stdout=subprocess.DEVNULL)
    os.close(fhf)
    if not save:
        os.close(fhv)
    return vcf_out


def gt_types_to_binary_comparison(calls):
    """From an array of calls, check if a variant position qualifies as a variant.

    0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
    Return string of 1s and 0s to represent position"""
    binary_calls = []
    for call in calls:
        if call == 1 or call == 3:
            binary_calls.append(1)
        else:
            binary_calls.append(0)
    return ''.join([str(i) for i in binary_calls])


def make_sets(vcf, names):
    """From the merged SV file, return pd.Series of overlapping sets"""
    calls = defaultdict(int)
    for v in VCF(vcf):
        calls[gt_types_to_binary_comparison(v.gt_types)] += 1
    tf_array = [[True, False]] * len(list(calls.keys())[0])
    index = pd.MultiIndex.from_product(tf_array, names=names)
    values = [calls[''.join([str(int(j)) for j in i])] for i in index]
    return pd.Series(values, index=index)


def get_args():
    parser = ArgumentParser(description="Calculate precision-recall metrics from 2 vcf files")
    parser.add_argument("--variants",
                        help="vcfs containing structural variants",
                        required=True,
                        nargs="*")
    parser.add_argument("--names",
                        help="Names of datasets in --variants",
                        nargs="*")
    parser.add_argument("-d", "--distance",
                        help="maximum distance between test and truth call",
                        default=500)
    parser.add_argument("--minlength",
                        help="Minimum length of SVs to be taken into account",
                        default=50)
    parser.add_argument("-i", "--ignore_type",
                        help="Ignore the type of the structural variant",
                        action="store_true")
    parser.add_argument("--store",
                        help="Save vcf of overlapping calls",
                        action="store_true")
    args = parser.parse_args()
    if args.names:
        if not len(args.variants) == len(args.names):
            sys.exit("INPUT ERROR: Need to have same number of values in --names as --variants!")
    return args


if __name__ == '__main__':
    main()
