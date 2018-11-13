from argparse import ArgumentParser
import csv
import matplotlib.pyplot as plt


class Callset(object):
    def __init__(self, aligner, caller, precision, recall, opacity):
        self.aligner = aligner
        self.caller = caller
        self.precision = float(precision)
        self.recall = float(recall)
        self.opacity = float(opacity)
        self.shape = None
        self.colour = None

    def plot(self):
        return plt.scatter(x=self.precision,
                           y=self.recall,
                           c=self.colour,
                           marker=self.shape,
                           alpha=self.opacity)


def main():
    args = get_args()
    calls = [Callset(*line) for line in csv.reader(open(args.results), delimiter="\t")]
    aligner_index, aligners = aligner_to_symbol(calls)
    caller_index, callers = caller_to_colour(calls)
    lines = [c.plot() for c in calls]
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    plt.xlabel('Precision')
    plt.ylabel('Recall')
    if len(aligners) > 1:
        plt.gca().add_artist(plt.legend(handles=[lines[i] for i in aligner_index],
                                        labels=aligners,
                                        loc='upper left',
                                        frameon=False))
    if len(callers) > 1:
        plt.gca().add_artist(plt.legend(handles=[lines[i] for i in caller_index],
                                        labels=callers,
                                        loc='lower left',
                                        frameon=False))
    plt.show()


def aligner_to_symbol(calls):
    """
    Assign symbols to different aligners in the input file
    Set the attribute of the class instances

    return a list of indices for which each aligner is found uniquely and all aligners
    sorted by aligners
    """
    symbols = ['o', '+', 'x', 'v', '*', 'D', 's', 'p', '8', 'X']
    aligners = sorted(set([c.aligner for c in calls]), reverse=True)
    aligner_to_symbol_dict = {a: s for a, s in zip(aligners, symbols)}
    for c in calls:
        c.shape = aligner_to_symbol_dict[c.aligner]
    index_and_aligners = zip([[c.aligner for c in calls].index(i) for i in aligners], aligners)
    return zip(*sorted(index_and_aligners, key=lambda x: x[1]))


def caller_to_colour(calls):
    """
    Assign colours to different callers in the input file
    Set the attribute of the class instances

    return a list of indices for which each caller is found uniquely and all callers
    sorted by callers
    """
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color'] * 10
    callers = sorted(set([c.caller for c in calls]), reverse=True)
    caller_to_colour_dict = {ca: co for ca, co in zip(callers, colours)}
    for c in calls:
        c.colour = caller_to_colour_dict[c.caller]
    index_and_callers = zip([[c.caller for c in calls].index(i) for i in callers], callers)
    return zip(*sorted(index_and_callers, key=lambda x: x[1]))


def get_args():
    parser = ArgumentParser(description="Plot precision and recall for multiple callsets")
    parser.add_argument(
        "results", help="A tsv file with aligner caller precision recall")
    return parser.parse_args()


if __name__ == '__main__':
    main()
