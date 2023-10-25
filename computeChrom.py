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

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--out",dest = "out")
parser.add_argument("-p", "--pheno",dest = "pheno")
parser.add_argument("-r", "--rfmix",dest = "rfmix")



args = parser.parse_args()

#read phenotypes and covariates
df_pheno = pd.read_csv(args.pheno,index_col="id")



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


#remove nan phenotypes
notin = []
phenos = np.array(df_pheno["PHENO"])
namess = list(df_pheno.index)
for i in range(len(phenos)):
    if np.isnan(phenos[i]):
        notin.append(namess[i])
lanc0 = lanc0.drop(notin,axis=1)
lanc1 = lanc1.drop(notin,axis=1)
df_pheno = df_pheno.reindex(list(lanc0.columns))

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

    
    

#drop sex if exclusive\
try:
    if sum(np.array(df_pheno["SEX"]))==0 or sum(np.array(df_pheno["SEX"]))==df_pheno.shape[0]:
        df_pheno = df_pheno.drop(["SEX"],axis=1)
except:
    pass
    
    
#association
df_assoc = admix.assoc.marginal(
    geno=lanc,
    lanc=lanc,
    n_anc=len(list(pd.read_csv(args.rfmix,nrows=0,sep="\t"))),
    pheno=np.array(df_pheno["PHENO"]),
    cov=df_pheno.drop(columns=["PHENO"]),
    method="ADM",
    family="linear",
)

    

#compile segments and reindex
df_assoc["id"] = pos
df_assoc = df_assoc.set_index("id")


compression_opts = dict(method='zip',
                        archive_name=f'assoc.csv')  
df_assoc.to_csv(f"{args.out}.zip", index=True, 
          compression=compression_opts) 