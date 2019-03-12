from cyvcf2 import VCF, Writer
from argparse import ArgumentParser


def main():
    args = get_args()
    vcf = VCF(args.vcf)
    w = Writer(args.output, vcf)
    for v in vcf:
        if v.INFO["SVTYPE"] == "DEL":
            if not v.INFO["SVLEN"] < 0:
                v.INFO["SVLEN"] = - v.INFO["SVLEN"]
        w.write_record(v)


def get_args():
    parser = ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="input vcf file", required=True)
    parser.add_argument("-o", "--output", help="output vcf file", required=True)
    return parser.parse_args()


if __name__ == '__main__':
    main()
