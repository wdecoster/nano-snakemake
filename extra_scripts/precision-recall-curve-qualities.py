import subprocess
from argparse import ArgumentParser
import shlex
import sys
import os
from cyvcf2 import VCF
import tempfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


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
    variants = extract_variants(combined_vcf)
    quals = get_qualities_from_test(args.test)
    variants["qual_test"] = match_qualities_with_coords(variants["coords_test"], quals)
    res = pd.DataFrame(data=[precision_recall(variants.loc[variants["qual_test"] > q])
                             for q in variants["qual_test"].unique()],
                       columns=["precision", "recall"])
    res.plot(x="precision", y="recall", kind="scatter")
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    plt.show()


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


def extract_variants(vcf):
    return pd.DataFrame(
        data=[tuple(v.format('CO')) + tuple([is_variant(call) for call in v.gt_types])
              + parse_qualities(v.format('QV')) for v in VCF(vcf)],
        columns=["coords_truth", "coords_test", "truth", "test", "quals_truth", "quals_test"])


def parse_qualities(quals):
    return tuple([np.median([int(q) for q in s.split(',')]) if s != 'NaN' else 0 for s in quals])


def precision_recall(df, qual):
    tests = (df["test"] & (df["quals_test"] > qual))
    tp = sum(df["truth"] & tests)
    if tp == 0:
        return 0, 0
    else:


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
    return parser.parse_args()


if __name__ == '__main__':
    main()
