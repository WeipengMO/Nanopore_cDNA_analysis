#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-18 12:13:33
@LastEditTime : 2020-04-18 15:03:28
@Description  : Writes a bed file from a file containing isoforms
'''


import click


@click.command()
@click.option('--infile', required=True)
@click.option('--outfile', required=True)
def convert_isoforms_to_bed(infile, outfile):
    with open(infile, "r") as inFile, open(outfile, "w") as outFile:
        for count, line in enumerate(inFile, 1):
            parts = line.rstrip().split("\t")
            chro = parts[0]
            start = parts[1]
            end = parts[2]
            # readCount=parts[3]
            spliceSites = parts[4].split(',')
            #Translate the splice sites into exons
            exons = []
            exon = [int(start),None]
            for spliceSite in spliceSites:
                if exon is None:
                    exon = [int(spliceSite), None]
                else:
                    exon[1] = int(spliceSite)
                    exons.append(exon)
                    exon=None
            exon[1] = int(end)
            exons.append(exon)
            #Create the blocksizes and blockstarts fields from the exons
            blockSizes=""
            blockStarts=""
            for exon in exons:
                blockSizes += f'{exon[1]-exon[0]},'
                blockStarts += f'{exon[0]-int(start)},'

            outFile.write(f'{chro}\t{start}\t{end}\tisoform_{count}\t0\t.\t{start}\t{end}\t'
                          f'0\t{len(exon)}\t{blockSizes}\t{blockStarts}\n')


if __name__ == "__main__":
    convert_isoforms_to_bed()