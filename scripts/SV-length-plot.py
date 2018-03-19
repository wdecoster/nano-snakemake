from cyvcf2 import VCF
import sys
import matplotlib.pyplot as plt


def main(vcf):
    plt.hist(x=[v.INFO.get('SVLEN') for v in VCF(vcf) if not v.INFO.get('SVTYPE') == 'TRA'],
             bins=[i for i in range(0, 2000, 50)])
    plt.xlabel('Lenghth of structural variant')
    plt.ylabel('Number of variants')
    plt.savefig(sys.argv[2])


if __name__ == '__main__':
    main(sys.argv[1])
