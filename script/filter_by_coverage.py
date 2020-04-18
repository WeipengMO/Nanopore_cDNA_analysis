#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-17 21:23:45
@LastEditTime : 2020-04-17 22:29:36
@Description  : Filters isoforms, removing the ones that are not represented enough compared to the coverage
'''


import click


@click.command()
@click.option('--raw_splice_sites', required=True)
@click.option('--coverage_file', required=True)
@click.option('--filtered_splice_sites', required=True)
@click.option('--max_cov', required=False, default=100)
def filterByCoverage(raw_splice_sites, coverage_file, filtered_splice_sites, max_cov):
    '''
    @ raw_splice_sites: chro, start, end, read_count, spliceSites
    @ coverage_file: chro, pos, coverage
    '''
    coverageMap = {}
    currentChr=None
    #Read the fill coverage file for quick access
    with open(coverage_file, 'r') as coverage:
        for line in coverage:
            chro, _, value = line.rstrip().split("\t")
            if chro != currentChr:
                coverageMap[chro]=[None] #the position 0 does not exist in the coverage file, so we fill it with None
                currentChr=chro
            coverageMap[chro].append(int(value))

    #Check the coverage at every splice position
    with open(raw_splice_sites, 'r') as rawIsomophs, open(filtered_splice_sites, 'w') as filteredIsoforms:
        for line in rawIsomophs:
            parts = line.rstrip().split("\t")
            chro=parts[0]
            readCount=int(parts[3])
            spliceSites=parts[4].split(',')
            valid=False
            for spliceSite in spliceSites:
                valid = valid or coverageMap[chro][int(spliceSite)] / readCount < max_cov
            if valid:
                filteredIsoforms.write(line)


if __name__ == '__main__':
    filterByCoverage()