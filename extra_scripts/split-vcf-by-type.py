from argparse import ArgumentParser
from cyvcf2 import VCF, Writer
from collections import defaultdict


def main():
    args = get_args()
    vars = defaultdict(list)
    vcf = VCF(args.vcf)
    output = args.vcf.replace('.vcf', '')
    for v in vcf:
        vars[v.INFO.get('SVTYPE')].append(v)
    for k, varlist in vars.items():
        w = Writer(output + '_' + k.replace('/', '') + '.vcf', vcf)
        for v in varlist:
            w.write_record(v)


def get_args():
    parser = ArgumentParser(description="Split a vcf by SV type")
    parser.add_argument("vcf", help="vcf file to split")
    return parser.parse_args()


if __name__ == '__main__':
    main()
