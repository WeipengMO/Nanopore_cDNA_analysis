#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-17 17:26:26
@LastEditTime : 2020-04-18 09:32:11
@Description  : Find splice sites in a bed file created with the -split option
'''

import click


def get_spliced_transcript(exons):
    start = None
    lastexon = None
    spliceSites = []
    for exon in exons:
        #the first part of each read should be the same chromosome name
        chro = exon[0]
        #the end position of the previous exon is the start of a splice site
        if lastexon is not None:
            spliceSites.append(lastexon[2])
        #the start position of the first exon should be the start of the transcript
        if start is None:
            start = exon[1]
        #the start position of the other exons are the end of splice sites
        else:
            spliceSites.append(exon[1])
        lastexon = exon
    #the remaining end position at the end of the loop is the end of the transcript
    end = lastexon[2]
    return [chro, start, end, spliceSites]
    

def get_isoform(transcripts):
    #if len(transcripts) <= 0:
    #    return["",0,0,0,[]]
    start = 0
    end = 0
    for transcript in transcripts:
        #the first part of each transcript should be the same chromosome name
        chro = transcript[0]
        #the final start and end positions will be the mean of all start and end positions, respectively
        start += int(transcript[1])
        end += int(transcript[2])
    #the last part of each transcript should be the same splice sites
    spliceSites=transcript[3]
    len_transcript = len(transcripts)
    return [chro, round(start/len_transcript), round(end/len_transcript), len_transcript-1, spliceSites]


@click.command()
@click.option('--infile', required=True)
@click.option('--outfile', required=True)
def find_splice_sites(infile, outfile):
    # read the file and compare each line with its neighbours
    # input: bed file
    # format: 0: chr  1: start  2: end  3: read_id  4: score  5: strand
    with open(infile, "r") as bedFile:
        firstline = True
        # bed文件上一行数据
        lastline = None
        exons = []
        splicedTranscirpts = []
        for line in bedFile:
            #Seach for duplicated instances of reads in the split file (split sites) by comparing their names
            currentline=line.rstrip().split("\t")
            # 如果为同一read_id, 则两行之间为splice区间
            if lastline is not None and lastline[3] == currentline[3]:
                if firstline:
                    exons.append(lastline)
                    firstline=False
                # 该read_id的exon位置
                exons.append(currentline)
            else:
                #If we found a duplicate, then we add it to the list along with its splice sites
                if not firstline:
                    splicedTranscirpts.append(get_spliced_transcript(exons))
                #When a new read name is found, reinitialize the variables
                exons=[]
                firstline=True
            lastline=currentline

        #If the last lines were a duplicate, then we add it to the list along with its splice sites
        # 处理最后一条read
        if not firstline:
            splicedTranscirpts.append(get_spliced_transcript(exons))

    #sort our duplicate list by key, then by chromosome name for later manipulation
    splicedTranscirpts.sort(key=lambda transcript: transcript[3])
    splicedTranscirpts.sort(key=lambda transcript: transcript[0])

    #Keep only the transcripts with at least 3
    isoforms = []
    currentIsoform = []
    firstline = True
    lastSpliceSites = None
    # splicedTranscripts is a list contains [chro, start, end, spliceSites]
    for transcript in splicedTranscirpts:
        spliceSites = " ".join(transcript[3])
        # 完全一样？无考虑测序错误情况
        # TODO(windz) 可优化
        if spliceSites == lastSpliceSites:
            currentIsoform.append(transcript)
        else:
            if firstline:
                firstline=False
            else:
                #Add the isoform to the list if there are at least 3 instances of it
                if len(currentIsoform) > 2:
                    isoforms.append(get_isoform(currentIsoform))

            currentIsoform=[transcript]
        # 上一个splice sites的位置(字符串)
        lastSpliceSites=spliceSites


    with open(outfile, "w") as outFile:
        for transcript in isoforms:
            # chro, start, end, length, spliceSites
            outFile.write(f'{transcript[0]}\t{transcript[1]}\t{transcript[2]}\t{transcript[3]}\t{",".join(transcript[4])}\n')


if __name__ == '__main__':
    find_splice_sites()