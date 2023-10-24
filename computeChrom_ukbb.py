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

#parse all the arguments
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--out",dest = "out")
parser.add_argument("-t", "--trait",dest = "trait")
parser.add_argument("-c", "--chrom",dest = "chrom")
parser.add_argument("-p", "--pheno",dest = "pheno",nargs='?', type=str, default="~/ukbb_project/phenotype-description.csv")
parser.add_argument("-l", "--lanc",dest = "lanc",nargs='?', type=str, default="/u/project/sgss/UKBB/UKB-ADMIXED/01-dataset/out/rfmix-hm3-lanc")
parser.add_argument("-v", "--cov",dest = "cov",nargs='?', type=str, default="/u/project/sgss/UKBB/PRS-RESEARCH/03-compile-pheno/out/")
parser.add_argument("-n", "--nanc",dest = "nanc",nargs='?', type=int, default=2)


args = parser.parse_args()
out = args.out
traitpos = int(args.trait) #traitpos here is the index of trait, not trait name


#read phenotypes and get the trait name
phenos = pd.read_csv(args.pheno)

try:
    temp = float(phenos["id"][traitpos-1])
    traitdir = re.sub(r'\W+', '_',np.array(phenos["description"])[traitpos-1])
except:
    traitdir = np.array(phenos["id"])[traitpos-1]
trait = np.array(phenos["id"])[traitpos-1]


if trait!=traitdir:
    print(f"Working on {traitpos}, {traitdir}_{trait}")
else:
    print(f"Working on {traitpos}, {traitdir}")


#iterate through 22 chromosomes

print(f"Running chromosome {args.chrom}")

#read lanc
path = f"{args.lanc}/chr{args.chrom}.msp.tsv"
df_rfmix = pd.read_csv(path, sep="\t", skiprows=1)
lanc0 = df_rfmix.loc[:, df_rfmix.columns.str.endswith(".0")].rename(columns=lambda x: x[:-2])
lanc1 = df_rfmix.loc[:, df_rfmix.columns.str.endswith(".1")].rename(columns=lambda x: x[:-2])
pos = list("chr"+df_rfmix["#chm"].astype(str)+":"+df_rfmix["spos"].astype(str)+"-"+df_rfmix["epos"].astype(str))

#read phenotypes and covariates
df_dir = args.cov
df_pheno = pd.read_csv(join(df_dir,f"{trait}.tsv"), delim_whitespace=True, low_memory=False,)
df_pheno.index = df_pheno.FID.astype(str) + "_" + df_pheno.IID.astype(str)
df_pheno = df_pheno.drop(columns=["FID", "IID"])
df_covar = pd.read_csv(join(df_dir, "covar.tsv"),delim_whitespace=True,low_memory=False)
df_covar.index = df_covar.FID.astype(str) + "_" + df_covar.IID.astype(str)
df_covar = df_covar.drop(columns=["FID", "IID"])

#find common indiv
notin = [x for x in list(lanc0.columns) if x not in df_pheno.index]
lanc0 = lanc0.drop(notin,axis=1)
lanc1 = lanc1.drop(notin,axis=1)
df_pheno = df_pheno.reindex(list(lanc0.columns))
df_covar = df_covar.reindex(list(lanc0.columns))

#remove nan phenotypes
notin = []
phenos = np.array(df_pheno)
namess = list(df_pheno.index)
for i in range(len(phenos)):
    if np.isnan(phenos[i][0]):
        notin.append(namess[i])
lanc0 = lanc0.drop(notin,axis=1)
lanc1 = lanc1.drop(notin,axis=1)
df_pheno = df_pheno.reindex(list(lanc0.columns))
df_covar = df_covar.reindex(list(lanc0.columns))

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

#drop sex if exclusive
if sum(np.array(df_covar["SEX"]))==0 or sum(np.array(df_covar["SEX"]))==df_covar.shape[0]:
    df_covar = df_covar.drop(["SEX"],axis=1)

#association
df_assoc = admix.assoc.marginal(
    geno=lanc,
    lanc=lanc,
    n_anc=int(args.nanc),
    pheno=np.array(df_pheno).flatten(),
    cov=df_covar,
    method="ADM",
    family="linear",
)

    

#compile segments and reindex
df_assoc["id"] = pos
df_assoc = df_assoc.set_index("id")


if traitdir!=trait:
    name = f"{traitdir}_{trait}_chr{args.chrom}"
else:
    name = f"{trait}_chr{args.chrom}"

compression_opts = dict(method='zip',
                        archive_name=f'{name}.csv')  
df_assoc.to_csv(f"{args.out}/{name}_assoc.zip", index=True, 
          compression=compression_opts) 