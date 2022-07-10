#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)

    fh.close()


def add_metavirome_tax(file):

    r = []
    data = {}

    for line in read_tsv(file, "\t"):
        r.append(line)
        if len(line) <=6:
            continue
        if "no support" in line[6]:
            continue
        temp = line[2:7]
        #print(line[6])
        #print(temp)
        if line[6] in data:
            LOG.info("\t".join(temp))
            continue
        if "no support" in temp:
            continue
        data[line[6]] = temp

    #print(data)
    print("#Contig\tSoftware\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies")
    for line in r:
        if len(line) <=6:
            print("\t".join(line))
            continue
        if line[6] in data:
            n = 1 
            for i in data[line[6]]:
                n += 1
                if "no support" in i:
                    continue
                if "no support" not in line[n]:
                    continue
                line[n] = i.strip("*")
        print("\t".join(line))

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="STR", type=str,
        help="Input metavirus species annotation summary file, unique_noncRNA_tax.tsv.")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
URL: https://github.com/zxgsy520/biodb
name:
    add_metavirome_tax: Merge the annotation results of each software virus classification

attention:
    add_metavirome_tax.py unique_noncRNA_tax.tsv >unique_noncRNA_tax_new.tsv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    add_metavirome_tax(args.input)


if __name__ == "__main__":

    main()
