import subprocess
from argparse import ArgumentParser
import shlex
import sys
import os
from cyvcf2 import VCF
import tempfile
import pandas as pd


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
    zygosity_accuracy(combined_vcf)


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


def get_zygosity(call):
    """Check if a variant position qualifies as a variant

    0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT"""
    zygdict = dict([(0, "nocall"), (1, "het"), (2, "nocall"), (3, "hom")])
    return zygdict[call]


def zygosity_accuracy(vcff):
    """
    First level of the dict is the "truth" call, second level is the "test"
    """
    zygosities = dict(het=dict(het=0,
                               hom=0,
                               nocall=0
                               ),
                      hom=dict(het=0,
                               hom=0,
                               nocall=0
                               ),
                      nocall=dict(het=0,
                                  hom=0,
                                  nocall=0
                                  ),
                      )
    for v in VCF(vcff):
        zyg_truth, zyg_test = [get_zygosity(gt) for gt in v.gt_types]
        zygosities[zyg_truth][zyg_test] += 1
    zygs = ["nocall", "het", "hom"]
    df = pd.DataFrame(index=zygs, columns=zygs)
    for tr in zygs:
        for te in zygs:
            df.loc[tr, te] = zygosities[tr][te]
    print(df)


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
