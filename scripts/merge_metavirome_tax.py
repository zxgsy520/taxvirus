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


def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".faa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    r = ""
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            if r:
                yield r.split("\n", 1)
            r = "%s\n" % line.strip(">")
            continue
        r += line.upper()

    if r:
        yield r.split("\n", 1)
    fp.close()


def read_phagcn2(file):

    r = {}
    for line in open(file):
        line = line.strip()
        if not line:
            continue
        line = line.split(",")
        r[line[0]] = line[2]

    return r


def read_vpf_class(files):

    r = {}
    for file in files:
        level = file.split("/")[-1].split(".")[0]

        for line in read_tsv(file, "\t"):
            if line[0] not in r:
                r[line[0]] = ["", ""]
            if "family" in level:
                r[line[0]][0] = line[1]
            if "genus" in level:
                r[line[0]][1] = line[1]

    return r


def read_blast(file):

    r = {}

    for line in read_tsv(file, "\t"):
        r[line[0]] = line[3::]

    return r


def read_cat(file):

    r = {}

    for line in read_tsv(file, "\t"):
        r[line[0]] = line[1::]

    return r


def merge_metavirome_tax(genome, cat, vpf_class, phagcn2, blast):

    dcat = read_cat(cat)
    dvpf = read_vpf_class(vpf_class)
    dphagcn = read_phagcn2(phagcn2)
    dblast = read_blast(blast)
    
    print("#Contig\tSoftware\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies")
    for seqid, seq in read_fasta(genome):
        seqid = seqid.split()[0]
        soft = ""
        family = ""
        genus = ""
        temp = ["no support", "no support", "no support", "no support", "no support", "", ""]
        if seqid in dcat:
            #print("%s\t%s" % (seqid, "\t".join(dcat[seqid])))
            if dcat[seqid]:
                if "Viruses" in dcat[seqid][0]:
                    temp = dcat[seqid]
                    soft = "CAT"

        if seqid in dvpf:
            #print("%s\t%s" % (seqid, "\t".join(dvpf[seqid])))
            family, genus = dvpf[seqid]
            soft = "%s;vpf-class" % soft
        if seqid in dphagcn:
            #print("%s\t%s" % (seqid, dphagcn[seqid]))
            soft = "%s;PhaGCN2" % soft
            family = dphagcn[seqid]
        if family:
            if "no support" in temp[0]:
                temp[0] = "Viruses"
                temp[4] = family
            else:
                if "no support" in temp[4]:
                    temp[4] = family

        if seqid in dblast:
            if dblast[seqid]:
                #print("%s\t%s\t%s" % (seqid, "blastn", "\t".join(dblast[seqid])))
                if ("Viruses" not in temp[0]) and ("Viruses" in dblast[seqid][0]):
                    soft = "blastn"
                    temp = dblast[seqid]
        soft = soft.strip(";")
        print("%s\t%s\t%s" % (seqid, soft, "\t".join(temp)))
               
                
    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="STR", type=str,
        help="Input Metaviral Genome Sequence(fasta).")
    parser.add_argument("-x", "--cat", metavar="FILE", type=str, required=True,
        help="Input cat virus classification results(contig2tax.tsv).")
    parser.add_argument("-v", "--vpf_class", nargs="+", metavar="FILE", type=str, required=True,
        help="Input vpf-class virus classification results(family.tsv genus.tsv).")
    parser.add_argument("-p", "--phagcn2", metavar="FILE", type=str, required=True,
        help="Input phagcn2 virus classification results(PhaGCN2.prediction.csv).")
    parser.add_argument("-b", "--blast", metavar="FILE", type=str, required=True,
        help="Input blast virus classification results(nt.tsv).")

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
    merge_metavirome_tax: Merge the annotation results of each software virus classification

attention:
    merge_metavirome_tax.py vog.proteins.all.fa --members vog.members.tsv.gz --annotations vog.annotations.tsv.gz >vogdb.fasta

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    merge_metavirome_tax(args.input, args.cat, args.vpf_class, args.phagcn2, args.blast)


if __name__ == "__main__":

    main()
