#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-17 11:13:07
@LastEditTime : 2020-04-17 11:43:47
@Description  : 
'''


import click
import pickle


@click.command()
@click.option('--infile', required=True)
@click.option('--outfile', required=True)
@click.option('--min_coverage_threshold', required=False, default=6)
def find_coverage_breaks(infile, outfile, min_coverage_threshold):
    if 'positive' in infile:
        strand = '+'
    elif 'negative' in infile:
        strand = '-'
    else:
        raise NameError('input error')
    
    chromosomes = {}   #The list of all chromosomes with genes position for each one
    currentGenes = []  #The list of gene positions for the currently parsed chromosome
    currentGene = {"start":None, "end":None, "strand":strand} #The current gene for which we search the position
    currentChr = None  #The name of the currently parsed chromosome
    geneNumber = 0  #The number of genes found for information

    with open(infile, "r") as coverage:
        for line in coverage:
            chro, pos, cov=line.rstrip().split("\t")
            #If we reach a new chromosome, we save the current genes positions in an array and begin a new one
            if chro != currentChr:
                if currentChr is not None: 
                    chromosomes[currentChr]=list(currentGenes)
                currentChr=chro
                currentGenes = []
                currentGene={"start":None, "end":None, "strand":strand}

            #In any case try to find the start of a new gene if we don't have one
            if currentGene["start"] is None:
                if int(cov) > min_coverage_threshold: 
                    currentGene["start"]=int(pos)
            #Or the end of the current gene if we already have a start
            else:
                #When we find the end of a gene, we register it in our list and start lookin for a new one
                if int(cov) < min_coverage_threshold:
                    currentGene["end"]=int(pos)
                    currentGenes.append(currentGene)
                    currentGene = {"start":None, "end":None, "strand":strand}
                    geneNumber = geneNumber+1

    #Add the last chromsome list
    if currentChr is not None: 
        chromosomes[currentChr] = list(currentGenes)

    #Write the output as a bed file for later use
    with open(outfile, "w") as bedGenes:
        for chro, genelist in chromosomes.items():
            for gene in genelist:
                bedGenes.write(f'{chro}\t{gene["start"]}\t{gene["end"]}\t{gene["strand"]}\n')

    #return chromosomes
    with open(outfile+'.pydata', 'wb') as o:
        pickle.dump(chromosomes, o)


if __name__ == '__main__':
    find_coverage_breaks()