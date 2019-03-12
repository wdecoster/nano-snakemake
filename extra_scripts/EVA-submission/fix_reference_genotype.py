from cyvcf2 import VCF, Writer
from pyfaidx import Fasta
from argparse import ArgumentParser


def main():
    args = get_args()
    genome = Fasta(args.genome)
    vcf = VCF(args.vcf)
    output = Writer(args.output, vcf)
    incorrect_reference = 0
    for v in vcf:
        ref_nucl = get_reference_nucleotide(v.CHROM, v.start, genome)
        if v.REF != ref_nucl:
            v.REF = ref_nucl
            incorrect_reference += 1
        output.write_record(v)
    print("Fixed {} positions".format(incorrect_reference))
    output.close()
    vcf.close()


def get_reference_nucleotide(chrom, position, fai):
    return fai[chrom][position].seq


def get_args():
    parser = ArgumentParser(description="")
    parser.add_argument("-g", "--genome", help="genome fasta to match", required=True)
    parser.add_argument("-v", "--vcf", help="input vcf file", required=True)
    parser.add_argument("-o", "--output", help="output vcf file", required=True)
    return parser.parse_args()


if __name__ == '__main__':
    main()
