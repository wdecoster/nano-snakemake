from cyvcf2 import VCF, Writer
from argparse import ArgumentParser


def main():
    args = get_args()
    vcf = VCF(args.vcf)
    output = Writer(args.output, vcf)
    incorrect = 0
    for v in vcf:
        if v.REF == v.ALT[0] and v.INFO["SVTYPE"] == "DEL":
            v.ALT = "<DEL>"
            incorrect += 1
        output.write_record(v)
    print("Fixed {} positions".format(incorrect))
    output.close()
    vcf.close()


def get_reference_nucleotide(chrom, position, fai):
    return fai[chrom][position].seq


def get_args():
    parser = ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="input vcf file", required=True)
    parser.add_argument("-o", "--output", help="output vcf file", required=True)
    return parser.parse_args()


if __name__ == '__main__':
    main()
