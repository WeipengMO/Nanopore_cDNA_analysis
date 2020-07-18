#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-05-24 19:07:01
@LastEditTime : 2020-07-11 19:05:12
@Description  : 
'''


import pickle
from collections import defaultdict, OrderedDict
import numpy as np
import pysam
from scipy import stats
import click
import pyranges as pr


class Polya_Cluster:
    def __init__(self, polya_sites, gene_id, chrom, is_reverse):
        polya_sites.sort()
        self.polya_sites = np.array(polya_sites)
        self.gene_id = gene_id
        self.chrom = chrom
        self.strand = '-' if is_reverse else '+'
        self.pac_list = None
        self._get_pac()

    def _get_pac(self):
        # 将24nt以内的点合并
        for polya_site in self.polya_sites:
            if self.pac_list is None:
                self.pac_list = [[polya_site]]
            else:
                if polya_site-self.pac_list[-1][-1] <= 24:
                    self.pac_list[-1].append(polya_site)
                else:
                    self.pac_list.append([polya_site])
        # 去除read数少于3的pac
        for i in range(len(self.pac_list)-1, -1, -1):
            if len(self.pac_list[i]) < 3:
                self.pac_list.pop(i)
        # 计算pac的poly(A) site,
        # a cluster supported by the greatest number
        major = float('-inf')
        self.polya_cluster = []
        self.polya_cluster_summit = []
        for n, pac in enumerate(self.pac_list, 1):
            start = pac[0]
            end = pac[-1]
            summit = stats.mode(pac)[0][0]
            self.polya_cluster.append(f'{self.chrom}\t{start-1}\t{end}\t{self.gene_id}_{n}\t{len(pac)}\t{self.strand}')
            self.polya_cluster_summit.append(f'{self.chrom}\t{summit-1}\t{summit}\t{self.gene_id}_{n}\t.\t{self.strand}')
            if len(pac) > major:
                major = len(pac)
                self.major_pac = f'{self.chrom}\t{start-1}\t{end}\t{self.gene_id}\t{major}\t{self.strand}'
                self.major_summit = f'{self.chrom}\t{summit-1}\t{summit}\t{self.gene_id}\t.\t{self.strand}'

        if len(self.polya_cluster) >= 1:
            if self.strand == '+':
                self.last_pac = self.polya_cluster[-1]
                self.last_pac_summit = self.polya_cluster_summit[-1]
            else:
                self.last_pac = self.polya_cluster[-0]
                self.last_pac_summit = self.polya_cluster_summit[-0]
            
            self.polya_cluster = '\n'.join(self.polya_cluster)
            self.polya_cluster_summit = '\n'.join(self.polya_cluster_summit)
        else:
            self.polya_cluster = None
            self.polya_cluster_summit = None
            self.major_pac = None
            self.major_summit = None
            self.last_pac = None
            self.last_pac_summit = None


@click.command()
@click.option('--infile', required=True)
@click.option('--gene_bed', required=True)
@click.option('--out_suffix', required=True)
def main(infile, gene_bed, out_suffix):

    # 1️⃣ 读取gene_bed信息
    gene_bed = pr.read_bed(gene_bed, as_df=True)

    # 2️⃣ 按照gene遍历bam文件，获取每个基因的pA site
    STRAND_TO_BOOL = {'-': True, '+': False}
    gene_polya_site = defaultdict(lambda: [])
    gene_info = OrderedDict()
    inbam = pysam.AlignmentFile(infile, 'rb')
    for item in gene_bed.itertuples():
        _, chrom, start, end, gene_id, _, strand = item
        strand = STRAND_TO_BOOL[strand]
        if chrom in {'Mt', 'Pt'}:
            continue
        read_count = 0
        for read in inbam.fetch(chrom, start, end):
            # 取该位置上与基因方向一致的reads
            if strand is read.is_reverse:
                polya_len = read.get_tag('pa')
                read_gene_id = read.get_tag('gi')
                if polya_len >= 15 and read_gene_id in {'None', gene_id}:
                    read_count += 1
                    if not read.is_reverse:
                        gene_polya_site[gene_id].append(read.reference_end)
                    else:
                        gene_polya_site[gene_id].append(read.reference_start)
        if read_count >= 5:
            gene_info[gene_id] = (chrom, strand)

    inbam.close()

    polya_cluster = open(f'{out_suffix}.polya_cluster.bed', 'w')
    polya_cluster_summit = open(f'{out_suffix}.polya_cluster.summit.bed', 'w')
    major_polya_cluster = open(f'{out_suffix}.major_polya_cluster.bed', 'w')
    major_polya_cluster_summit = open(f'{out_suffix}.major_polya_cluster_summit.bed', 'w')
    last_polya_cluster = open(f'{out_suffix}.last_polya_cluster.bed', 'w')
    last_polya_cluster_summit = open(f'{out_suffix}.last_polya_cluster_summit.bed', 'w')
    for gene_id in gene_info:
        pac = Polya_Cluster(gene_polya_site[gene_id], gene_id, *gene_info[gene_id])
        if pac.polya_cluster is not None and pac.major_pac is not None:
            polya_cluster.write(pac.polya_cluster+'\n')
            polya_cluster_summit.write(pac.polya_cluster_summit+'\n')
            major_polya_cluster.write(pac.major_pac+'\n')
            major_polya_cluster_summit.write(pac.major_summit+'\n')
            last_polya_cluster.write(pac.last_pac+'\n')
            last_polya_cluster_summit.write(pac.last_pac_summit+'\n')
    
    polya_cluster.close()
    polya_cluster_summit.close()
    major_polya_cluster.close()
    major_polya_cluster_summit.close()
    last_polya_cluster.close()
    last_polya_cluster_summit.close()


if __name__ == "__main__":
    main()