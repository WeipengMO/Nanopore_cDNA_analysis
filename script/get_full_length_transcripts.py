#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-15 15:26:26
@LastEditTime : 2020-04-15 15:32:29
@Description  : get full length transcripts from bam
'''


import pysam
import pandas as pd
import click


@click.command()
@click.option('-i', '--infile', required=True)
@click.option('--full_len', required=True)
@click.option('--non_full_len', required=True)
@click.option('--first_exon_path', required=True)
def main(infile, full_len, non_full_len, first_exon_path):
    #first_exon_path = '/public/home/mowp/test/nanopore_cdna/supplementary_data/representative_gene_model/representative_gene_first_exon.tsv'
    first_exon_df = pd.read_csv(first_exon_path, sep='\t', index_col=0)
    first_exon_dict = first_exon_df.to_dict(orient='index')

    #infile = '/public/home/mowp/test/nanopore_cdna/aligned_data/fhh.tagged.mm2.sorted.bam'
    with pysam.AlignmentFile(infile, 'rb') as inbam:
        full_len_bam = pysam.AlignmentFile(full_len, 'wb', template=inbam)
        non_full_len_bam = pysam.AlignmentFile(non_full_len, 'wb', template=inbam)
        for read in inbam:
            read_gene_id = read.get_tags()[-2][1]
            if read_gene_id in first_exon_dict:
                if (first_exon_dict[read_gene_id]['strand'] == '+' and 
                    read.reference_start <= first_exon_dict[read_gene_id]['end']):
                    full_len_bam.write(read)
                elif (first_exon_dict[read_gene_id]['strand'] == '-' and
                    read.reference_end >= first_exon_dict[read_gene_id]['start']):
                    full_len_bam.write(read)
                else:
                    non_full_len_bam.write(read)
    
    full_len_bam.close()
    non_full_len_bam.close()

    
if __name__ == "__main__":
    main()