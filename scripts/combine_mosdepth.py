#! /usr/bin/env python

import pandas as pd
from argparse import ArgumentParser


def main():
    args = get_args()
    d = pd.concat([pd.read_csv(filepath_or_buffer=f,
                               sep="\t",
                               header=None,
                               names=['chr', 'begin', 'end', f.replace(".regions.bed.gz", "")])
                   for f in args.inputfiles],
                  axis="columns")
    d.loc[:, ~d.columns.duplicated()] \
        .to_csv(path_or_buf=args.output,
                sep="\t",
                header=True,
                index=False,
                compression="gzip")


def get_args():
    parser = ArgumentParser(
        description="combine gzipped mosdepth region results, assuming same bins")
    parser.add_argument("-o", "--output", help="outputfile to use", required=True)
    parser.add_argument("inputfiles", help="mosdepth region results", nargs='+')
    return parser.parse_args()


if __name__ == '__main__':
    main()
