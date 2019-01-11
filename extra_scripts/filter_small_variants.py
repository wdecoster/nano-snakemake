from cyvcf2 import VCF, Writer
from argparse import ArgumentParser


def main():
    args = get_args()
    vcf_in = VCF(args.vcf)
    vcf_out = Writer(args.output, vcf_in)
    for v in vcf_in:
        if v.INFO["SVLEN"] > 49:
            vcf_out.write_record(v)
    vcf_in.close()
    vcf_out.close()


def get_args():
    parser = ArgumentParser(
        description="Remove SVs smaller than 50bp from a vcf")
    parser.add_argument("vcf", help="vcf to start from")
    parser.add_argument("output", help="vcf to produce")
    return parser.parse_args()


if __name__ == '__main__':
    main()
