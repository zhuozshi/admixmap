##########################################################################################
### FileName: getPeaks.py
### Author: Zhuozheng Shi
### Date: 12/07/2023
### Affiliation: UCLA Bogdan Lab
### Description: Generate table of peaks for one trait, containing starting and ending 
###              locations as well as separate chromosomes and base 
### Input: [assocFiles] paths for association file
###        [threshold] p value threshold, e.g. 5e-5
###        [output] output file path for peaks table
### Output: table containing starting and ending locations as well as separate chromosomes
###         and positions
##########################################################################################

import numpy as np
import pandas as pd
import math
import sys
import argparse
import openpyxl
import os
import copy
import time
import re
import structlog

def getBand(chrom, pos, dic):    
    
    chrom = int(chrom)
    pos = int(pos)
    
    temp = dic[dic.iloc[:,0]==f"chr{chrom}"]
    
    temp = temp[temp.iloc[:,1]<pos]

    temp = temp[temp.iloc[:,2]>pos]
    
    
    return f"{chrom}{temp.iloc[0,3]}"
    

logger = structlog.get_logger()
args = sys.argv
args = [0,"interpolation/AoU/AoU_GIA_afreur_height_adm_chrALL_interpolated.csv.gz",5e-5,"0"]
    
    
threshold = float(args[2])
logger.info(
        f"AssocFiles {args[1]} with thresholds {threshold}"
    )

assoc = pd.read_csv(args[1],index_col=0)


peakTable = pd.DataFrame([],columns=["peak start","peak end","chr","pos start","pos end","# regions"])




assoc["chr"] = [int(x.split(":")[0]) for x in assoc.index]
assoc["pos"] = [int(x.split(":")[1]) for x in assoc.index]

flag = 0
toappend = []
for i in range(len(assoc)):

    index = assoc.index[i]
    if assoc.P.iloc[i] <=threshold:
        if flag:
            if int(index.split(sep=":")[0])!=toappend[2]:
                peakTable.loc[len(peakTable)] = toappend
                toappend = [assoc.index[i],assoc.index[i],int(index.split(sep=":")[0]),int(index.split(sep=":")[1]),int(index.split(sep=":")[1]),1]
            else:
                toappend[1] = index
                toappend[4] = int(index.split(sep=":")[1])
                toappend[5] = toappend[5] + 1
                continue
        else:
            toappend = [assoc.index[i],assoc.index[i],int(index.split(sep=":")[0]),int(index.split(sep=":")[1]),int(index.split(sep=":")[1]),1]
            flag = 1
            num = 1
    else:
        if flag:
            flag = 0
            peakTable.loc[len(peakTable)] = toappend
    

    
logger.info(
        f"Table with {peakTable.shape[0]} peaks saved to {args[3]}.csv.gz"
    )
peakTable.to_csv(f"{args[3]}.csv.gz",index=True)
