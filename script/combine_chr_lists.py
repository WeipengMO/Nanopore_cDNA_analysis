#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-17 16:56:35
@LastEditTime : 2020-04-17 17:17:32
@Description  : Combines two chromosomes list together respecting the ordering
'''

import pickle
import click


@click.command()
@click.option('--positive', required=True)
@click.option('--negative', required=True)
@click.option('--outfile', required=True)
def combineChrLists(positive, negative, outfile):
    with open(f'{positive}.pydata', 'rb') as f:
        chrList1 = pickle.load(f)
    with open(f'{negative}.pydata', 'rb') as f:
        chrList2 = pickle.load(f)

    mergedChrList = {}
    geneNumber = 0
    for chro, geneList1 in chrList1.items():
        mergedGeneList = []
        for gene in geneList1:
            #chrList2[chr][0]["start"] is the start position of the first gene left in chrList2 for the current chromosome
            #While there are genes in chrList2 that start before the current gene, add them to the mergedGeneList and remove them from chrList2
            while len(chrList2[chro]) > 0 and chrList2[chro][0]["start"] < gene["start"]:
                mergedGeneList.append(chrList2[chro].pop(0))
                geneNumber += 1
            #When no genes are left that start before the current one, add it to the mergedGeneList
            mergedGeneList.append(gene)
            geneNumber += 1
        #When all genes from the current chromosome have been added for chrList1, add any remaining genes from chrList2 to the mergedGeneList
        for gene in chrList2[chro]:
            mergedGeneList.append(gene)
            geneNumber = geneNumber+1

        #Finaly, add the mergedGeneList to the mergedChrList
        mergedChrList[chro]=mergedGeneList

    #The mergedChrList should now contain every gene for every chromosome sorted by chromosome and position

    #Write the output as a bed file
    with open(outfile, "w") as bedGenes:
        for chro, genelist in mergedChrList.items():
            for gene in genelist:
                bedGenes.write(f'{chro}\t{gene["start"]}\t{gene["end"]}\t{gene["strand"]}\n')


if __name__ == "__main__":
	combineChrLists()