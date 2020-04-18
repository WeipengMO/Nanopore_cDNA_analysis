#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-14 10:46:37
@LastEditTime : 2020-04-14 17:35:31
@Description  : add polya length and remove wrong strandness reads
@usage        : python add_tag_to_bam -i <in.bam> -o <out.bam> -g <genome.gff> -p <polya.tsv>
'''


import pysam
import pandas as pd
import numpy as np
import click


def parse_attribute(attr):
    '''get gene_id from attribute
    '''
    for item in attr.split(';'):
        key, value = item.split('=')
        if key == 'ID':
            return value.replace('gene:', '')


def get_gene_is_reverse(gff_path):
    GFF_COLUMNS = ['seqid', 'source', 'type', 'start', 'end', 
                   'score', 'strand', 'phase', 'attribute']
    gff_df = pd.read_csv(gff_path, sep='\t', comment='#', names=GFF_COLUMNS)
    # keep type = 'gene'
    gff_df = gff_df[gff_df['type'].isin(['gene'])]
    # get gene_id
    gff_df['gene_id'] = gff_df['attribute'].map(parse_attribute)
    gff_df['is_reverse'] = gff_df['strand'].map(lambda x: False if x=='+' else True)
    
    return gff_df.set_index('gene_id')['is_reverse'].T.to_dict()


def load_polya_data(polya_data):
    # polya_data = '/public/home/mowp/test/fast5_api/fhh.polyacaller.ann.tsv'
    df = pd.read_csv(polya_data, sep='\t', index_col=0)
    df['read_is_reverse'] = df['read_type'].map(lambda x: True if x == 'polyT' else False)

    # a dict like this: {'read_id': True or False}
    read_is_reverse = df['read_is_reverse'].to_dict()
    # a dict like this: {'read_id': 'gene_id'}
    read_gene_id = df['gene_id'].to_dict()
    # a dict like this: {'read_id': 100.0}
    read_tail_length = df['tail_length'].to_dict()
    return read_is_reverse, read_gene_id, read_tail_length


@click.command()
@click.option('-i', '--infile', required=True)
@click.option('-o', '--outfile', required=True)
@click.option('-p', '--polya_file', required=True)
@click.option('-g', '--gff_file', required=True)
def main(infile, outfile, polya_file, gff_file):
    #gff_file = '/public/home/mowp/db/Arabidopsis_thaliana/gff3/Arabidopsis_thaliana.TAIR10.46.gff3'
    gene_is_reverse_dict = get_gene_is_reverse(gff_file)

    #polya_file = '/public/home/mowp/test/fast5_api/fhh.polyacaller.ann.tsv'
    read_is_reverse_dict, read_gene_id_dict, read_tail_length_dict = load_polya_data(polya_file)

    #infile = '/public/home/mowp/test/nanopore_cdna/aligned_data/fhh.adjust.mm2.sorted.bam'
    #outfile = '/public/home/mowp/test/nanopore_cdna/aligned_data/fhh.tagged.mm2.sorted.bam'
    with pysam.AlignmentFile(infile, 'rb') as inbam, pysam.AlignmentFile(outfile, 'wb', template=inbam) as outbam:
        # add tags: pa for polya_length
        #           gi for gene_id
        for read in inbam.fetch():
            # only output mapped reads
            if not read.is_unmapped:
                read_is_reverse = read_is_reverse_dict[read.query_name]
                read_gene_id = read_gene_id_dict[read.query_name]
                read_tail_length = read_tail_length_dict[read.query_name]
                # 能被注释到gene_id的执行如下操作
                try:
                    gene_is_reverse = gene_is_reverse_dict[read_gene_id]
                    read.set_tag('gi', read_gene_id)
                    if gene_is_reverse is read.is_reverse:
                        read.set_tag('pa', read_tail_length)
                    else:
                        continue
                # 不能注释的直接输出，有可能是新的转录本
                except KeyError:
                    read.set_tag('gi', 'None')
                    read.set_tag('pa', read_tail_length)
                outbam.write(read)


if __name__ == "__main__":
    main()