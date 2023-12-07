##########################################################################################
### FileName: getPeaks.py
### Author: Zhuozheng Shi
### Date: 11/15/2023
### Affiliation: UCLA Bogdan Lab
### Description: Generate table of peaks for association files, containing peak positions
###              and peak bands under strict and loose p value thresholds.
### Input: [assocFiles Directory] directory paths for all association files
###        [bandsDict] dictionary file path for bands and peaks
###        [strict p value threshold] strict p value threshold, e.g. 5e-5
###        [loose p value threshold] loose p value threshold, e.g. 5e-4
###        [output] output file path for peaks table
### Output: table containing peak positions and peak bands under strict and loose p value 
###         thresholds.
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
    
    
threshold_strict = float(args[3])
threshold_loose = float(args[4])
dic = pd.read_csv(args[2])
#ukbb = os.listdir("interpolation/ukbb_v3")
ukbb = os.listdir(args[1])
ukbb.sort()
logger.info(
        f"AssocFiles Directory {args[1]}, band dictionary {args[2]}"
        f"with thresholds {threshold_strict} and {threshold_loose}"
    )



peakTable = pd.DataFrame([],columns=["trait",f"peaks p<{args[3]}",f"pos p<{args[3]}",f"peaks p<{args[4]}",f"pos p<{args[4]}"])


for i in ukbb:
    
    assoc = pd.read_csv(f"{args[1]}/{i}",index_col=0)

    assoc["chr"] = [int(x.split(":")[0]) for x in assoc.index]
    assoc["pos"] = [int(x.split(":")[1]) for x in assoc.index]
    peaks = pd.DataFrame(columns=assoc.columns)
    peakl = pd.DataFrame(columns=assoc.columns)
    
    peak = assoc[assoc.P<threshold_strict]
    bands = []
    while peak.shape[0]!=0:
        
        index = peak.idxmin(axis=0).P
        peaks.loc[index] = peak.loc[index]
        bands.append(getBand(peak.loc[index]["chr"],peak.loc[index]["pos"],dic))

        
        chrom = peak[peak["chr"]==peak.loc[index]["chr"]]
        temp = chrom[chrom["pos"]>chrom.loc[index]["pos"]-1.5e6]
        chrom = temp[temp["pos"]<temp.loc[index]["pos"]+1.5e6]
        indices=list(chrom.index)
        peak = peak.drop(indices)
    peaks["bands"] = bands
    
    peak = assoc[assoc.P<threshold_loose]
    bands = []
    while peak.shape[0]!=0:
        
        index = peak.idxmin(axis=0).P
        peakl.loc[index] = peak.loc[index]
        bands.append(getBand(peak.loc[index]["chr"],peak.loc[index]["pos"],dic))

        
        chrom = peak[peak["chr"]==peak.loc[index]["chr"]]
        temp = chrom[chrom["pos"]>chrom.loc[index]["pos"]-1.5e6]
        chrom = temp[temp["pos"]<temp.loc[index]["pos"]+1.5e6]
        indices=list(chrom.index)
        peak = peak.drop(indices)
    peakl["bands"] = bands
        
    
    
    #peaks = list(peaks["bands"])
    #peakl = list(peakl["bands"])
    sss = pd.unique(np.array(peaks["bands"],dtype="str"))
    lll = pd.unique(np.array(peakl["bands"],dtype="str"))
    spos = np.array(peaks.index)
    lpos = np.array(peakl.index)
    
    peakTable.loc[len(peakTable)] = [i,sss,spos,lll,lpos]
    #peakTable.loc[len(peakTable)] = [i[:-26],sss,spos,lll,lpos]
    
peakTable = peakTable.set_index("trait")
logger.info(
        f"Table with {peakTable.shape[0]} traits saved to {args[5]}.csv.gz"
    )
peakTable.to_csv(f"{args[5]}.csv.gz",index=True)
