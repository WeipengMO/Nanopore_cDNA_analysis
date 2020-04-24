#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-24 09:28:51
@LastEditTime : 2020-04-24 22:28:28
@Description  : 
'''


import numpy as np
from glob import glob
from collections import Counter, namedtuple, defaultdict
from operator import attrgetter
import heapq
import pysam
import pandas as pd


def parse_exons_introns_flank(record, flanksize=200):
    start = int(record[1])
    end = int(record[2])
    exstarts = np.fromstring(record[11], sep=',') + start
    exends = exstarts + np.fromstring(record[10], sep=',')
    exons = np.dstack([exstarts, exends])[0]
    left_flank = np.array([[max(0, start - flanksize), start]])
    right_flank = np.array([[end, end + flanksize]])
    if len(exons) > 1:
        introns = np.dstack([exons[:-1, 1], exons[1:, 0]])[0]
    else:
        introns = np.array([])
    return exons, introns, left_flank, right_flank


def split_intervals(invs, pos, side='left'):
    idx = np.searchsorted(invs.ravel(), pos)
    split = np.insert(invs.ravel(), idx, [pos, pos]).reshape(-1, 2)
    split_idx = (idx + 1) // 2
    return split[:split_idx], split[split_idx:]


def parse_cds_utr_introns_flank(record, flanksize):
    exons, introns, left_flank, right_flank = parse_exons_introns_flank(record, flanksize)
    cds_start = int(record[6])
    cds_end = int(record[7])
    if not cds_start == cds_end:
        utr1, cds = split_intervals(exons, cds_start)
        cds, utr2 = split_intervals(cds, cds_end)
    else:
        utr1 = np.array([])
        cds = np.array([])
        utr2 = np.array([])
    return utr1, cds, utr2, introns, left_flank, right_flank, exons


def parse_features(record, flanksize=500):
    ''' 
    Parse bed format

    flanksize: up/downstream distance

    Return:
        {'chrom': '1',
        'strand': '+',
        'gene_id': 'spliced_9162',
        'invs': {
                'cds': array([[3695., 3695.], [3695., 3913.]]),
                'introns': array([[3913., 3995.], [5326., 5438.]]),
                'exons': array([[3695., 3913.], [5438., 5884.]]),
                '5utr': array([], shape=(0, 2), dtype=float64),
                '3utr': array([[5884., 5884.]]),
                'upstream': array([[3195, 3695]]),
                'downstream': array([[5884, 6384]])
                }
        }
    '''
    features = {}
    invs = {}
    features['chrom'] = record[0].replace('Chr', '')
    features['strand'] = record[5]
    features['gene_id'] = record[3]
    utr1, invs['cds'], utr2, invs['introns'], left_flank, right_flank, invs['exons'] = parse_cds_utr_introns_flank(record, flanksize)
    if features['strand'] == '+':
        invs['5utr'] = utr1
        invs['3utr'] = utr2
        invs['upstream'] = left_flank
        invs['downstream'] = right_flank
    else:
        invs['5utr'] = utr2
        invs['3utr'] = utr1
        invs['upstream'] = right_flank
        invs['downstream'] = left_flank
    features['invs'] = invs
    return features


def get_lengths_for_norm():
    feat_lengths = Counter()
    with open(genes_bed) as bed:
        for record in bed:
            record = parse_features(record.split())
            if record['chrom'] in ['C', 'M']:
                continue
            for feat_type, invs in record['invs'].items():
                for inv in invs:
                    feat_lengths[feat_type] += (inv[1] - inv[0])
    return pd.Series(feat_lengths) / 1000


def intersect(inv_a, inv_b):
    a_start, a_end = inv_a
    b_start, b_end = inv_b
    if a_end < b_start or a_start > b_end:
        return 0
    else:
        s = max(a_start, b_start)
        e = min(a_end, b_end)
        return e - s


def intersect_spliced_invs(invs_a, invs_b):
    '''
    Input:
        invs_a: bam invs
        invs_b: bed invs

    Return:
        overlap bases between invs
    '''
    score = 0
    invs_a = iter(invs_a)
    invs_b = iter(invs_b)
    a_start, a_end = next(invs_a)
    b_start, b_end = next(invs_b)
    while True:
        # 判断有无overlap
        if a_end < b_start:
            try:
                a_start, a_end = next(invs_a)
            except StopIteration:
                break
        elif a_start > b_end:
            try:
                b_start, b_end = next(invs_b)
            except StopIteration:
                break
        else:
            score += intersect([a_start, a_end], [b_start, b_end])
            if a_end > b_end:
                try:
                    b_start, b_end = next(invs_b)
                except StopIteration:
                    break
            else:
                try:
                    a_start, a_end = next(invs_a)
                except StopIteration:
                    break
    return score


class MultiBam(object):

    def __init__(self, bam_fns):
        # 用pysam打开
        # 以字典形式存储 fn -> pysam对象
        self.bam_handles = {bam_fn: pysam.AlignmentFile(bam_fn) for bam_fn in bam_fns}
        self.closed = False

    def fetch(self, *args, **kwargs):
        queries = [bam.fetch(*args, **kwargs) for bam in self.bam_handles.values()]
        yield from heapq.merge(*queries, key=attrgetter('reference_start'))

    def close(self):
        for bam in self.bam_handles.values():
            bam.close()
    
    # 上下文管理特殊方法
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()



def bam_cigar_to_invs(aln, max_allowed_insertion):
    invs = []
    start = aln.reference_start
    end = aln.reference_end
    strand = '-' if aln.is_reverse else '+'
    left = start
    right = left
    has_ins = False
    for op, ln in aln.cigartuples:
        if op in (4, 5):
            # cigartuples: 
            #   S  BAM_CSOFT_CLIP  4
            #   H  BAM_CHARD_CLIP  5
            # does not consume reference
            continue
        elif op == 1 and ln > max_allowed_insertion:
            # cigartuples: 
            #   I  BAM_CINS  1
            has_ins = True
        elif op in (0, 2, 7, 8):
            # cigartuples: 
            #   M  BAM_CMATCH  0
            #   D  BAM_CDEL  2
            #   =  BAM_CEQUAL  7
            #   X  BAM_CDIFF  8
            # consume reference but do not add to invs yet
            right += ln
        elif op == 3:
            # cigartuples: 
            #   N  BAM_CREF_SKIP  3
            invs.append([left, right])
            left = right + ln
            right = left

        # ❗ 以下if语句应该在循环内
        if right > left:
            invs.append([left, right])
    assert invs[0][0] == start
    assert invs[-1][1] == end
    return start, end, strand, np.array(invs), has_ins


PARSED_ALN = namedtuple('Aln', 'chrom start end read_id strand invs')

def parse_pysam_aln(aln, max_allowed_insertion):
    chrom = aln.reference_name
    read_id = aln.query_name
    start, end, strand, invs, has_ins = bam_cigar_to_invs(
        aln, max_allowed_insertion)
    return PARSED_ALN(chrom, start, end, read_id, strand, invs), has_ins


def assign_three_prime_to_feature(three_prime_end, bed_record):
    '''
    判断3'site位于该representative gene model的哪个区间
    '''
    if bed_record['strand'] == '+' and three_prime_end >= bed_record['invs']['exons'][-1][1]:
        return 'downstream'
    elif bed_record['strand'] == '-' and three_prime_end < bed_record['invs']['exons'][0][0]:
        return 'downstream'
    for feature_type in ['3utr', 'cds', 'introns', '5utr']:
        invs = bed_record['invs'][feature_type]
        for start, end in invs:
            if start <= three_prime_end < end:
                return feature_type
    else:
        assert False


def count_three_prime_ends_in_features(annotation_bed_fn, bam_fns):
    feature_pos_counts = defaultdict(Counter)
    feature_read_counts = defaultdict(Counter)
    with open(annotation_bed_fn) as bed, MultiBam(bam_fns) as bam:
        # 💡 遍历每个基因
        for record in bed:
            record = parse_features(record.split())  # record: bed file record
            gene_tpe = Counter()
            if not len(record['invs']['cds']):
                # not protein coding, continue:
                continue
            record_span = (record['invs']['exons'][0][0], record['invs']['exons'][-1][1])
            # 💡 在bam文件中遍历该基因exon区间
            for aln in bam.fetch(record['chrom'], *record_span):
                aln, has_ins = parse_pysam_aln(aln, 30)
                if has_ins:
                    continue
                elif aln.strand != record['strand']:
                    continue
                # aln.infer_query_length ?
                aln_len = sum([e - s for s, e in aln.invs])
                i = intersect_spliced_invs(aln.invs, record['invs']['exons'])
                # overlap rate > 0.2
                if i / aln_len > 0.2:
                    tpe = aln.start if aln.strand == '-' else aln.end
                    # 3'sites position count += 1
                    gene_tpe[tpe] += 1
            # 💡 遍历上述筛选过的3'site
            for tpe, count in gene_tpe.items():
                feat_type = assign_three_prime_to_feature(tpe, record)
                feature_pos_counts[record['gene_id']][feat_type] += 1
                feature_read_counts[record['gene_id']][feat_type] += count
    feature_pos_counts = pd.DataFrame.from_dict(feature_pos_counts, orient='index').fillna(0)
    feature_read_counts = pd.DataFrame.from_dict(feature_read_counts, orient='index').fillna(0)
    return feature_pos_counts, feature_read_counts


ARAPORT = '/cluster/ggs_lab/mtparker/Arabidopsis_annotations/Araport/v11/201606/Araport11_GFF3_genes_transposons.flat_genes.bed'

three_prime_pos_counts, three_prime_read_counts = count_three_prime_ends_in_features(
    ARAPORT,
    glob('../chimeric_transcripts/vir1_vs_col0/aligned_data/201*_col0_*.bam') + \
    ['/cluster/ggs_lab/mtparker/ONT_guppy_pipeline_runs/20180411_1432_20180911_FAH84603_5adapterLIG_Col0_2916/aligned_data/TAIR10/201902_col0_2916_5adapter_exp2.bam',
     '/cluster/ggs_lab/mtparker/ONT_guppy_pipeline_runs/20180508_1522_20180508_FAH82422_5adapt_lig_mRNA_2918/aligned_data/TAIR10/201902_col0_2918_5adapter.bam']
)

three_prime_counts = three_prime_counts.assign(total=three_prime_counts[['downstream', '3utr', 'cds', 'introns', '5utr']].sum(1))
three_prime_counts = three_prime_counts.assign(
    downstream_or_3utr_percent=(three_prime_counts.downstream + three_prime_counts['3utr']) / three_prime_counts.total * 100,
    before_3utr_percent=(three_prime_counts.cds + three_prime_counts.introns + three_prime_counts['5utr']) / three_prime_counts.total * 100,
    utr5_percent=three_prime_counts['5utr'] / three_prime_counts.total * 100,
    cds_percent=three_prime_counts['cds'] / three_prime_counts.total * 100,
    introns_percent=three_prime_counts['introns'] / three_prime_counts.total * 100,
)
three_prime_counts.head()
