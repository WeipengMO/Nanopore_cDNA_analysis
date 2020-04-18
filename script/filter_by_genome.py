#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-17 22:30:41
@LastEditTime : 2020-04-18 10:18:10
@Description  : #Filters isoforms, removing the ones that are not represented enough compared to the coverage
'''

import click
from Bio import SeqIO


@click.command()
@click.option('--raw_isoform_file', required=True)
@click.option('--genome_file', required=True)
@click.option('--filtered_isoform_file', required=True)
@click.option('--max_polya_length', required=False, default=8)  # 可选参数
@click.option('--min_exon_length', required=False, default=100)  # 可选参数
@click.option('--max_single_base_to_length_ratio', required=False, default=0.7)  # 可选参数
def filter_by_genome(raw_isoform_file, genome_file, filtered_isoform_file,
                     max_polya_length, min_exon_length, max_single_base_to_length_ratio):
    genomeMap = {}
    #Read the fill genome file for quick access
    for read in SeqIO.parse(genome_file, 'fasta'):
        chro = read.id
        seq = str(read.seq)
        genomeMap[chro] = seq

    polyA = "A" * max_polya_length
    polyT = "T" * max_polya_length

    #Check the genome for poly-A and poly-T content within each exon and delete the exon if it is found
    '''
    @raw_isoform_file format:
        chr, start, end, cov, splice_site
        1  999887.4  1000413.8  29  1000026,1000111
    '''
    with open(raw_isoform_file, 'r') as rawIsomophs, open(filtered_isoform_file,'w') as filteredIsoforms:
        for line in rawIsomophs:
            parts = line.strip().split("\t")
            chro = parts[0]
            start = parts[1]
            end = parts[2]
            readCount = parts[3]
            spliceSites = parts[4].split(',')
            #Translate the splice sites into exon addresses
            exons=[]
            # exon: 5'ss和3'ss的位置
            exon=[start, None]
            for spliceSite in spliceSites:
                if exon is None:
                    exon = [spliceSite, None]
                else:
                    exon[1] = spliceSite
                    exons.append(exon)
                    exon = None
            exon[1] = end
            exons.append(exon)

            #Look for poly-A, poly-T or suspiciously high A or T count in the exon
            spliceSites=[]
            start = None
            for exon in exons:
                sequence=genomeMap[chro][int(exon[0]): int(exon[1])]
                # TODO 经验值，按需调整
                # 过滤短exon，已经含有多个A/T的exon
                if (len(sequence) < min_exon_length
                    and (polyA in sequence or polyT in sequence 
                         or sequence.count('A')/len(sequence) > max_single_base_to_length_ratio 
                         or sequence.count('T')/len(sequence) > max_single_base_to_length_ratio)):
                    pass
                else:
                    if start is None:
                        start = exon[0]
                    else:
                        spliceSites.append(exon[0])
                    spliceSites.append(exon[1])
            if len(spliceSites) > 2:
                end = spliceSites.pop()
                filteredIsoforms.write("\t".join((chro, start, end, readCount, ",".join(spliceSites)))+"\n")


if __name__ == "__main__":
    filter_by_genome()