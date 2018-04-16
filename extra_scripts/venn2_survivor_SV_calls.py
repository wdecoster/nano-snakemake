from cyvcf2 import VCF
import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


def is_variant(call):
    if call == 1 or call == 3:
        return True
    else:
        return False


def main(vcf):
    identifiers = [[], []]
    for v in VCF(vcf):
        for index, call in enumerate(v.gt_types):
            if is_variant(call):
                identifiers[index].append(v.ID)
    venn2([set(i) for i in identifiers], set_labels=('NanoSV', 'Sniffles'))
    plt.savefig(sys.argv[1] + ".png")


if __name__ == '__main__':
    main(sys.argv[1])
