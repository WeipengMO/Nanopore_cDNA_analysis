#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-17 12:01:25
@LastEditTime : 2020-04-17 16:45:35
@Description  : Mixes a list of 3' or 5' results into a list of chromosomes start and end sites as created by findCoverageBreaks
'''


import pickle
import click


CONFIG = {'min_3dash_threshold': 10, 'min_3dash_peak_separation': 30,
          'min_3dash_end_separation': 100, 'min_5dash_threshold': 10,
          'min_5dash_peak_separation': 30, 'min_5dash_start_separation': 100, 'max_5dash_offset': 200}


def clusterMean(cluster):
    totalSites=0
    totalPosition=0
    #The mean is calculated proportionnaly to the number of sites or each position
    for site in cluster:
        totalSites = totalSites+int(site["coverage"])
        totalPosition = totalPosition+(int(site["position"])*int(site["coverage"]))

    return round(totalPosition/totalSites)


@click.command()
@click.option('--in3dash_file', required=True)
@click.option('--in5dash_file', required=True)
@click.option('--outfile', required=True)
@click.option('--chromosome_data', required=True)
def insert_dash_cuts(in3dash_file, in5dash_file, chromosome_data, outfile):
    '''
    "positive_3dash_file", "positive_5dash_file", positiveChrList, "positive_genelist_file", "+")
    '''
    if 'positive' in in3dash_file and 'positive' in in5dash_file:
        strand = '+'
    elif 'negative' in in3dash_file and 'negative' in in5dash_file:
        strand = '-'
    else:
        raise NameError('input error')
    
    with open(f'{chromosome_data}.pydata', 'rb') as f:
        chromosomeList = pickle.load(f)

    for dash in ('3', '5'):
        #Read the 3' or 5' coverage file o find clusters of start sites or end sites
        currentDashSite = {"chr":None, "position":None, "coverage":None}  #The line read currently
        lastDashSite = None  #The last line that had a coverage over theset threshold
        currentCluster = []  #The current cluster of results
        currentChr = None  #The current chromosome name
        currentChrList = []  #The list 3' and 5' sites for the current chromosome
        dashSiteList = {}  #The full list of relevant 3' or 5' sites
        if dash == "3":
            sourceFile = in3dash_file
        elif dash == "5":
            sourceFile = in5dash_file
        with open(sourceFile, "r") as dashList:
            for line in dashList:
                currentDashSite["chr"],  currentDashSite["position"],  currentDashSite["coverage"] = line.rstrip().split("\t")
                #First check if we are entering a new chromosome
                if currentDashSite["chr"] != currentChr:
                    #if so, store the list of sites for the previous chromosome, then empty the list and reset the variables
                    if currentChr is not None:
                        if currentCluster != []:
                            currentChrList.append(clusterMean(currentCluster))
                        dashSiteList[currentChr] = list(currentChrList)
                    currentCluster = []
                    currentChr=currentDashSite["chr"]
                    currentChrList = []
                    lastDashSite = None

                #lines with a coverage less than a threshold (min_3dash_threshold) are ignored
                if int(currentDashSite["coverage"]) > CONFIG[f"min_{dash}dash_threshold"]:
                    #Peaks within a distance of each other (min_3dash_peak_separation) are clustered together
                    if lastDashSite is None or int(currentDashSite["position"])-int(lastDashSite["position"]) < CONFIG[f"min_{dash}dash_peak_separation"]:
                        currentCluster.append(dict(currentDashSite))
                    #If the peak is far enough from the last cluster, register the mean position for the peak and reset the cluster with our new result
                    else:
                        #The mean position is stored in the current chromosome list
                        if currentCluster != []:
                            currentChrList.append(clusterMean(currentCluster))
                        # 重置，计算下一个cluster
                        currentCluster = []
                        currentCluster.append(dict(currentDashSite))
                    #Once the site has been saved, keep it to compare it to the next
                    lastDashSite = dict(currentDashSite)

        #Add the last chromosome list
        if currentChr is not None:
            currentChrList.append(clusterMean(currentCluster))
            dashSiteList[currentChr] = list(currentChrList)
        #The dashSiteList now contains every cut site given by our 3' or 5' file, organized by chomosome name
        # dashSiteList 为上面结果汇总，为每一挑染色体上，每一个cluster的均值

        #Copy the dash site list to the correct variable
        if dash == "3":
            dashSiteList3=dict(dashSiteList)
        elif dash == "5":
            dashSiteList5=dict(dashSiteList)

    #Compare the dashSiteList with our chromosomeList chromosome by chromosome
    newChrList = {}
    geneNumber = 0
    for chro, geneList in chromosomeList.items():
        #Compare each cut site with each gene site
        newGeneList = []
        for gene in geneList:
            for cut3 in dashSiteList3[chro]:
                #If the 3 dash cut is within the start and end of a gene and separated enough from the end of it,
                #cut the gene accordingly and look for a 5' site to start the next gene
                if gene["start"] < cut3 < gene["end"]:
                    # if dash == "3" and gene["end"]-cut > int(config["min_3dash_end_separation"]):
                    if gene["end"]-cut3 > CONFIG["min_3dash_end_separation"]:
                        #Create a new gene ending at the corresponding 3' position
                        newGeneList.append({"start":gene["start"], "end":cut3, "strand":strand})
                        geneNumber = geneNumber+1
                        gene["start"] = cut3
                        #Look for a 5' site to start the next gene
                        for cut5 in dashSiteList5[chro]:
                            #The 5' site must be between a reasonable length from the 3' cut and the end of the gene
                            if cut3 - CONFIG["max_5dash_offset"] < cut5 < gene["end"]:
                                #Reduce the current gene to start at the 5' just found
                                gene["start"] = cut5
                                break

            #The renaining gene i added to the list at the end
            newGeneList.append(gene)
            geneNumber = geneNumber+1
        #The newGeneList is added to our newChrList
        newChrList[chro]=list(newGeneList)

    #Write the output as a bed file
    with open(outfile, "w") as bedGenes:
        for chro, genelist in newChrList.items():
            for gene in genelist:
                bedGenes.write(f'{chro}\t{gene["start"]}\t{gene["end"]}\t{gene["strand"]}\n')
    
    #return newChrList
    with open(f'{outfile}.pydata', 'wb') as o:
        pickle.dump(newChrList, o)


if __name__ == '__main__':
    insert_dash_cuts()