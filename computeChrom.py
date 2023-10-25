##########################################################################################
### FileName: computeChrom.py
### Author: Zhuozheng Shi
### Date: 10/25/2023
### Affiliation: UCLA Bogdan Lab
### Description: Compute the admixture mapping using rfmix local ancestry segments.
### Input: [-o, --out] output file name including path prefix, the output file will be a 
###                    zip file containing the pd.DataFrame of result, which should have
###                    the format "{path}/{prefix}_{chromosome_index(1-22)}"
###        [-p, --pheno] input phenotype and covariates pd.DataFrame file path, the first  
###                      column must be in the format of "FID_IID", and the second column
###                      must contain all the phenotypes the rest of column will be all
###                      the covariates
###        [-r, --rfmix] input rfmix file ".msp.tsv" path.
### Output: pd.DataFrame at name and path given as [-o, --out], "{out}.csv.gz", 
###         containning "L1_BETA", "L1_SE", "N", and "P" for local ancestry segments in 
###         rfmix file, with "id" as index.
##########################################################################################

import matplotlib.pyplot as plt
import admix
import numpy as np
from admix.plot import compare_pval
from os.path import join
import pandas as pd
import math
import sys
import argparse
import openpyxl
import os
import statsmodels.api as sm
import statsmodels.formula.api as smf
import copy
import statsmodels
import time
import dask.array as da
import qtl
from qtl.plot import qqplot
import matplotlib.image as mpimg
from scipy import stats 
import re
import structlog

logger = structlog.get_logger()
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--out",dest = "out")
parser.add_argument("-p", "--pheno",dest = "pheno")
parser.add_argument("-r", "--rfmix",dest = "rfmix")
args = parser.parse_args()
logger.info(
        f"Currently computing chromosome {args.rfmix}"
    )

#read phenotypes and covariates
df_pheno = pd.read_csv(args.pheno,index_col="id")
logger.info(
        f"Reading phenotypes and covariates in {args.pheno} with shape {df_pheno.shape}\n"
        f"First 5 lines of data \n{df_pheno.iloc[:5]}"
    )


#read lanc
df_rfmix = pd.read_csv(args.rfmix, sep="\t", skiprows=1)
lanc0 = df_rfmix.loc[:, df_rfmix.columns.str.endswith(".0")].rename(columns=lambda x: x[:-2])
lanc1 = df_rfmix.loc[:, df_rfmix.columns.str.endswith(".1")].rename(columns=lambda x: x[:-2])
pos = list("chr"+df_rfmix["#chm"].astype(str)+":"+df_rfmix["spos"].astype(str)+"-"+df_rfmix["epos"].astype(str))

#find common indiv
notin = [x for x in list(lanc0.columns) if x not in df_pheno.index]
lanc0 = lanc0.drop(notin,axis=1)
lanc1 = lanc1.drop(notin,axis=1)
df_pheno = df_pheno.reindex(list(lanc0.columns))
if len(notin)>0:
    logger.info("Found mismatched phenotype and lanc. Removed")

#remove nan phenotypes
notin = []
phenos = np.array(df_pheno.iloc[:,0])
namess = list(df_pheno.index)
for i in range(len(phenos)):
    if np.isnan(phenos[i]):
        notin.append(namess[i])
lanc0 = lanc0.drop(notin,axis=1)
lanc1 = lanc1.drop(notin,axis=1)
df_pheno = df_pheno.reindex(list(lanc0.columns))
if len(notin)>0:
    logger.info("Found nan phenotypes. Removed")

#compile lanc
l1 = lanc0.to_numpy()
l2 = lanc1.to_numpy()

n_snp,n_indiv = l1.shape

lanc = []
for i in range(n_snp):
    ins = []
    for j in range(n_indiv):
        ins.append([l1[i,j],l2[i,j]])
    lanc.append(ins)

#drop covariates if exclusive\
todrop = []
for c in range(1,df_pheno.shape[1]):
    if len(set(df_pheno.iloc[:,c])) == 1:
        todrop.append(df_pheno.columns[c])
for d in todrop:
    df_pheno = df_pheno.drop([d],axis=1) 
if len(todrop)>0:
    logger.info(f"Found invalid covariates {todrop}. Removed")

    
#association
df_assoc = admix.assoc.marginal(
    geno=lanc,
    lanc=lanc,
    n_anc=len(list(pd.read_csv(args.rfmix,nrows=0,sep="\t"))),
    pheno=np.array(df_pheno.iloc[:,0]),
    cov=df_pheno.iloc[:,1:],
    method="ADM",
    family="linear",
)

#compile segments and reindex
df_assoc["id"] = pos
df_assoc = df_assoc.set_index("id")
df_assoc.to_csv(f"{args.out}.csv.gz", index=True) 

logger.info(
        f"Saving to {args.out}.csv.gz"
    )