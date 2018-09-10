from cyvcf2 import VCF
from argparse import ArgumentParser


def main():
    args = get_args()
    overlapping_gene = 0
    overlapping_exon = 0
    total = 0
    for v in VCF(args.vcf):
        total += 1
        if v.INFO.get('CODING'):
            overlapping_exon += 1
        if v.INFO.get('GENE'):
            overlapping_gene += 1
    print("Structural variants detected {}".format(total))
    print("Structural variants overlapping a gene {}".format(overlapping_gene))
    print("Structural variants overlapping an exon {}".format(overlapping_exon))


def get_args():
    parser = ArgumentParser(description="calculate metrics based on an annotated vcf of SVs")
    parser.add_argument("vcf", help="vcf to calculate metrics on")
    return parser.parse_args()


if __name__ == '__main__':
    main()
