#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-18 10:37:08
@LastEditTime : 2020-04-18 11:55:06
@Description  : Creates a file containing only the most frequent isomoprhs from a file containing a list of them
'''


import click


@click.command()
@click.option('--infile', required=True)
@click.option('--outfile', required=True)
@click.option('--max_splice_difference', required=False, default=10)  # 可选参数
def find_best_isoform(infile, outfile, max_splice_difference):
    #read the file and compare each line with its neighbours
    # infile format: chro, start, end, cov, splice_sites
    #                1  999887  1000414  29  1000026,1000111
    with open(infile, "r") as inFile, open(outfile, "w") as outFile:
            lastline = None
            isoforms = []
            for line in inFile:
                #Seach for instances sharing approximately the same splice site
                currentline = line.rstrip().split("\t")
                currentline += currentline.pop().split(',')
                
                # TODO 只看前两个splice site，是否可以全部看？
                # TODO 如果第一个splice site相同，但后面有不同的转录本就看不了
                if (lastline is not None
                    and (abs(int(lastline[4]) - int(currentline[4])) < max_splice_difference 
                         or abs(int(lastline[5]) - int(currentline[5])) < max_splice_difference)):
                        isoforms.append(currentline)
                else:
                    #Once a new transcript is hit, write the best isoform of the list to the output file
                    if lastline is not None:
                        outFile.write("\t".join(getBest(isoforms))+"\n")
                    #Then reinitialize with the curent line
                    isoforms = [currentline]
                lastline = currentline

            #Write the output for the last batch of isoforms
            outFile.write("\t".join(getBest(isoforms))+"\n")


# Returns the best isoform of a list
# 取cov最大的作为best
def getBest(isoforms):
    best = None
    bestScore = float('-inf')
    for isoform in (isoforms):
        score = int(isoform[3])
        if best is None or score > bestScore:
            best = isoform
            bestScore = score
            
    splice_sites = ','.join(best[4:])
    best = best[:4] + [splice_sites]

    return best


if __name__ == "__main__":
    find_best_isoform()