from argparse import ArgumentParser
import csv
import matplotlib.pyplot as plt


class Callset(object):
    def __init__(self, aligner, caller, precision, recall, shape, colour):
        self.aligner = aligner
        self.caller = caller
        self.precision = int(precision)
        self.recall = int(recall)
        self.shape = shape
        self.colour = colour

    def plot(self, axes):
        axes.scatter(self.precision, self.recall, c=self.colour, marker=self.shape)


def main():
    args = get_args()
    calls = [Callset(*line) for line in csv.reader(open(args.results), delimiter="\t")]
    fig, axes = plt.subplots(1, 1)
    for c in calls:
        c.plot(axes)
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    plt.xlabel('Precision')
    plt.ylabel('Recall')
    # markers = plt.legend()
    plt.show()


def get_args():
    parser = ArgumentParser(description="Plot precision and recall for multiple callsets")
    parser.add_argument(
        "results", help="A tsv file with aligner caller precision recall shape colour")
    return parser.parse_args()


if __name__ == '__main__':
    main()
