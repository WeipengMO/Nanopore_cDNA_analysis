#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-17 17:26:26
@LastEditTime : 2020-04-23 19:17:57
@Description  : Find intronless transcript in a bed file created with the -split option
'''


from collections import namedtuple
import numpy as np
from scipy import stats
import click


import logging
logging.basicConfig(level=logging.DEBUG,  
                    format='%(asctime)s %(filename)s: %(message)s',  
                    datefmt='%m/%d/%Y %I:%M:%S %p',  
                    filename='find_intronless_transcripts.log',  
                    filemode='w')


class Cluster:
    '''A cluster of intronless transcript.

    Attributes:
        __chro: chromosome of the cluster
        __start: sum of start sites of transcripts
        __end: sum of start sites of transcripts
        __strand: strand of the cluster
        __count: read counts in the cluster
        __max_difference: max difference between sites
    '''
    def __init__(self, transcript):
        self.__chro = transcript[0]
        self.__start = [int(transcript[1])]
        self.__end = [int(transcript[2])]
        self.__strand = transcript[3] 
        self.__read_id = [transcript[5]] 
        self.__count = 1 
        self.__max_difference = 10

    def add_transcript(self, transcript):
        '''Add transcript to cluster
        ''' 
        self.__start.append(int(transcript[1]))
        self.__end.append(int(transcript[2]))
        self.__count += 1
        self.__read_id.append(transcript[5])

    def in_cluster(self, transcript):
        '''Determine if the transcript belong to cluster
        '''
        # only intronless transcript
        if transcript[3] != self.__strand and transcript[4] != '':
            return False

        start = np.array(self.__start)
        end = np.array(self.__end)
        transcript_start = int(transcript[1])
        transcript_end = int(transcript[2])
        
        #if ((abs(start-transcript_start) < self.__max_difference).any()
        #    and (abs(end-transcript_end) < self.__max_difference).any()):
        #  transcript 3'site within max_difference of the cluster
        #  只考虑3'sites, 因为5'sites不准确
        if ((self.__strand == '-' and (abs(start-transcript_start) < self.__max_difference).any())
            or (self.__strand == '+' and (abs(end-transcript_end) < self.__max_difference).any())):
            return True
        else:
            return False
    
    def __call__(self):
        '''Return isoform info of the cluster
        '''
        chro = self.__chro
        start = str(stats.mode(self.__start)[0][0])
        end = str(stats.mode(self.__end)[0][0])
        strand = self.__strand
        count = str(self.__count)
        spliced_sites = ''
        read_id = ','.join(self.__read_id)

        return chro, start, end, strand, count, spliced_sites, read_id


def get_spliced_transcript(last_line, spliced_sites_list):
    chro = last_line.chro
    start = spliced_sites_list[0]
    end = spliced_sites_list[-1]
    read_id = last_line.read_id
    strand = last_line.strand
    spliced_sites = ','.join(spliced_sites_list[1: -1])

    return chro, start, end, strand, spliced_sites, read_id


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
    intronless_transcirpts = []
    chromosome = set()
    # Seach for duplicated instances of reads in the split file (split sites) by comparing their names
    logging.debug('Loading infile')
    with open(infile, 'r') as bedFile:
        for line in bedFile:
            chro, start, end, read_id, score, strand = line.rstrip().split('\t')
            current_line = Bedformat(chro, start, end, read_id, score, strand)
            # 下一个read_id初始化
            if last_line is not None and current_line.read_id != last_line.read_id:
                if len(spliced_sites_list) <= 2:
                    intronless_transcirpts.append(get_spliced_transcript(last_line, spliced_sites_list))
                spliced_sites_list = []
            spliced_sites_list.append(current_line.start)
            spliced_sites_list.append(current_line.end)
            last_line = current_line
            chromosome.add(chro)

        # 处理最后一条read
        if len(spliced_sites_list) <= 2:
            intronless_transcirpts.append(get_spliced_transcript(last_line, spliced_sites_list))

    total_count = len(intronless_transcirpts)

    
    # Init result dict
    isoform_dict = {}
    for chro in chromosome:
        isoform_dict[chro] = []
    
    logging.debug('Clustering intronless transcripts')
    # TODO 可优化
    max_search = 100
    for count, transcript in enumerate(intronless_transcirpts, 1):
        if count % 10000 == 0:
            logging.debug(f'Treating {count} of {total_count} reads')

        chro = transcript[0]
        need_append = True
        for search, isoform in enumerate(isoform_dict[chro][::-1], 1):
            if search > max_search:
                break
            if isoform.in_cluster(transcript):
                isoform.add_transcript(transcript)
                need_append = False
                break
        if need_append:
            isoform = Cluster(transcript)
            isoform_dict[chro].append(isoform)
    
    logging.info('Starting output')
    with open(outfile, 'w') as out:
        for chro in sorted(isoform_dict):
            for isoform in isoform_dict[chro]:
                isoform = isoform()
                if int(isoform[4]) >= min_transcript_count:
                    out.write('\t'.join(isoform)+'\n')


if __name__ == '__main__':
    find_splice_sites()