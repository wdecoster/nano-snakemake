#! /usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import numpy as np
import sys


def main():
    args = get_args()
    res = pd.concat([read_tg(f, n) for f, n in zip(args.input, args.names or args.input)]) \
        .drop_duplicates(subset=['identifier', 'sample'])
    pivot = res.pivot(index='identifier',
                      columns='sample',
                      values='copy_number')
    pivot.loc[:, "max_change"] = pivot.abs().max(axis="columns")
    if args.sort == 'location':
        sorter = ["chromosome", "begin", "end"]
        ascending = True
    else:
        sorter = 'max_change'
        ascending = False
    print(res[["identifier", "chromosome", "begin", "end", "motif", "gene", "element"]]
          .drop_duplicates()
          .set_index('identifier')
          .join(pivot)
          .sort_values(sorter, ascending=ascending)
          .to_csv(sep="\t", na_rep='.')
          )


def read_tg(tgfile, name):
    df = pd.read_csv(tgfile,
                     sep="\t",
                     comment='#',
                     header=None,
                     converters={'motif': truncate_motif},
                     names=['chromosome', 'begin', 'end', 'motif', 'gene', 'element',
                            'copy_numbers_fwd', 'copy_numbers_rev'])
    df["copy_number"] = (df["copy_numbers_fwd"] + ',' + df["copy_numbers_rev"]
                         ).apply(parse_copy_number_changes)
    df["identifier"] = df["chromosome"] + ":" + df["begin"].astype(str) + "-" + df["motif"]
    df["sample"] = name
    return df


def truncate_motif(motif):
    if len(motif) > 8:
        return motif[:8] + '~'
    else:
        return motif


def parse_copy_number_changes(value):
    return np.median([int(i) if i != '.' else 0 for i in value.split(',')])


def get_args():
    parser = ArgumentParser(description="reformat and combine output of tandem-genotypes")
    parser.add_argument("--input", help="tandem-genotypes output files", nargs='+')
    parser.add_argument("--names", help="name of each dataset", nargs='+')
    parser.add_argument("--sort",
                        help="property to sort on. Default: chromosomal location",
                        choices=['location', 'max'],
                        default='location')
    args = parser.parse_args()
    if not args.names or not len(args.input) == len(args.names):
        sys.exit("ERROR: Need same number of values in --names as tandem-genotypes output files")
    return args


if __name__ == '__main__':
    main()
