#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-17 17:26:26
@LastEditTime : 2020-04-20 22:32:26
@Description  : Find splice sites in a bed file created with the -split option
'''


import logging
from collections import namedtuple
import numpy as np
import click


def get_spliced_transcript(last_line, spliced_sites_list):
    chro = last_line.chro
    start = spliced_sites_list[0]
    end = spliced_sites_list[-1]
    read_id = last_line.read_id
    strand = last_line.strand
    spliced_sites = ','.join(spliced_sites_list[1: -1])

    return chro, start, end, strand, spliced_sites, read_id


def is_cluster(curr_transcript, last_transcript, max_difference):
    # 判断两条reads是否为同一个isoform cluster
    curr_transcript = np.fromstring(curr_transcript, sep=',', dtype=int)
    last_transcript = np.fromstring(last_transcript, sep=',', dtype=int)
    if len(curr_transcript) != len(last_transcript):
        # 如果splice sites数目不一致，认为不是同一个cluster
        return False
    res = (abs(curr_transcript - last_transcript) < max_difference).all()
    return res


def get_isoform(transcripts):
    # 计算这一isoform cluster的 start 和 end 的加权平均值
    start, end = 0, 0
    read_id = []
    spliced_sites = None
    for count, transcript in enumerate(transcripts, 1):
        start += int(transcript[1])
        end += int(transcript[2])
        read_id.append(transcript[5])
        if spliced_sites is None:
            spliced_sites = np.fromstring(transcript[4], sep=',', dtype=int)
        else:
            spliced_sites += np.fromstring(transcript[4], sep=',', dtype=int)

    chro = transcript[0]
    start = round(start/count)
    end = round(end/count)
    strand = transcript[3]
    spliced_sites = np.round(spliced_sites/count).astype(int)
    spliced_sites = ','.join([str(item) for item in spliced_sites])
    read_id = ','.join(read_id)

    return '\t'.join([chro, str(start), str(end), strand, str(count), spliced_sites, read_id])


logging.basicConfig(format='[%(asctime)s]  %(message)s')


@click.command()
@click.option('--infile', required=True)
@click.option('--outfile', required=True)
@click.option('--max_difference', required=False, default = 10)
@click.option('--min_transcript_count', required=False, default = 3)
def find_splice_sites(infile, outfile, max_difference, min_transcript_count):
    # read the file and compare each line with its neighbours
    # infile: bedfile
    # format: 0: chro  1: start  2: end  3: read_id  4: score  5: strand
    Bedformat = namedtuple('Bedformat', 'chro start end read_id score strand')

    last_line = None
    spliced_sites_list = []
    spliced_transcirpts = []
    intronless_transcirpts = []
    # Seach for duplicated instances of reads in the split file (split sites) by comparing their names
    logging.info('start load infile')
    with open(infile, 'r') as bedFile:
        for line in bedFile:
            chro, start, end, read_id, score, strand = line.rstrip().split('\t')
            current_line = Bedformat(chro, start, end, read_id, score, strand)
            # 下一个read_id初始化
            if last_line is not None and current_line.read_id != last_line.read_id:
                if len(spliced_sites_list) > 2:
                    spliced_transcirpts.append(get_spliced_transcript(last_line, spliced_sites_list))
                else:
                    intronless_transcirpts.append(get_spliced_transcript(last_line, spliced_sites_list))
                spliced_sites_list = []
            spliced_sites_list.append(current_line.start)
            spliced_sites_list.append(current_line.end)
            last_line = current_line

        # 处理最后一条read
        spliced_transcirpts.append(get_spliced_transcript(last_line, spliced_sites_list))

    #sort our duplicate list by key, then by chromosome name for later manipulation
    spliced_transcirpts.sort(key=lambda transcript: transcript[3]) # sort by strand
    spliced_transcirpts.sort(key=lambda transcript: transcript[4]) # sort by splice sites
    spliced_transcirpts.sort(key=lambda transcript: transcript[0]) # sort by chr

    logging.info('start output')
    with open(outfile, 'w') as out:
        # for intron-contained transcripts
        last_transcript = None
        current_isoform = []
        for transcript in spliced_transcirpts:
            if last_transcript is not None and not is_cluster(transcript[-2], last_transcript[-2], max_difference):
                # 有3条reads支持的isoform (min_transcript_count)
                if len(current_isoform) >= min_transcript_count:
                    out.write(get_isoform(current_isoform)+'\n')
                # 重新初始化
                current_isoform = [transcript]
            elif last_transcript is not None and is_cluster(transcript[-2], last_transcript[-2], max_difference):
                current_isoform.append(transcript)
            elif last_transcript is None:
                current_isoform.append(transcript)
            last_transcript = transcript
        if len(current_isoform) >= min_transcript_count:
            out.write(get_isoform(current_isoform)+'\n')
        
        # for intronless transcripts
        last_transcript = None
        current_isoform = []
        for transcript in intronless_transcirpts:
            current_position = np.array((int(transcript[1]), int(transcript[2])))
            if (last_transcript is not None and 
                    (last_transcript[3] != transcript[3] 
                        or not(abs(last_position - current_position) < max_difference).all())):
                if len(current_isoform) >= min_transcript_count:
                    out.write(get_isoform(current_isoform)+'\n')
                current_isoform = []
            current_isoform.append(transcript)
            last_position = current_position
            last_transcript = transcript
        if len(current_isoform) >= min_transcript_count:
            out.write(get_isoform(current_isoform)+'\n')


if __name__ == '__main__':
    find_splice_sites()