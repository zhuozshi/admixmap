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

args = parser.parse_args()
out = args.out
traitpos = int(args.trait)

if False:
    lis = [11, 25, 26, 40]
    traitpos = lis[traitpos-1]
    
if traitpos<=106:
    traitdir = re.sub(r'\W+', '_',np.array(pd.read_csv("~/ukbb_project/phenotype-description.csv")["description"])[traitpos-1])
else:
    traitdir = np.array(pd.read_csv("~/ukbb_project/phenotype-description.csv")["id"])[traitpos-1]
trait = np.array(pd.read_csv("~/ukbb_project/phenotype-description.csv")["id"])[traitpos-1]

if traitpos<=106:
    print(f"Working on {traitpos}, {traitdir}_{trait}")
else:
    print(f"Working on {traitpos}, {traitdir}")



dfs = []
pos = []

for i in range(1,23):
    
    print(f"Running chromosome {i}")
    
    #read lanc
    path = f"/u/project/sgss/UKBB/UKB-ADMIXED/01-dataset/out/rfmix-hm3-lanc/chr{i}.msp.tsv"
    df_rfmix = pd.read_csv(path, sep="\t", skiprows=1)
    lanc0 = df_rfmix.loc[:, df_rfmix.columns.str.endswith(".0")].rename(columns=lambda x: x[:-2])
    lanc1 = df_rfmix.loc[:, df_rfmix.columns.str.endswith(".1")].rename(columns=lambda x: x[:-2])
    pos.append(list("chr"+df_rfmix["#chm"].astype(str)+":"+df_rfmix["spos"].astype(str)+"-"+df_rfmix["epos"].astype(str)))
    
    #read phenotypes and covariates
    df_dir = "/u/project/sgss/UKBB/PRS-RESEARCH/03-compile-pheno/out/"
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
        n_anc=2,
        pheno=np.array(df_pheno).flatten(),
        cov=df_covar,
        method="ADM",
        family="linear",
    )
    
    dfs.append(df_assoc)
    

    
df = pd.concat(dfs)
poss = [item for row in pos for item in row]
df["id"] = poss
df = df.set_index("id")

chrn = []
for i in range(0,22):
    for j in range(dfs[i].shape[0]):
        chrn.append(i+1)

        
        
fig, ax = plt.subplots( figsize=(5, 3), dpi=150)
#axes[0].scatter(np.arange(dset.n_snp), sim_beta[:, 0], s=2)

admix.plot.manhattan(np.array(df.P), ax=ax, s=2,chrom=np.array(chrn),axh_y=-np.log10(5e-5))
plt.tight_layout()
if traitpos<=106:
    plt.title(f"Association of {traitdir}_{trait}")
else:
    plt.title(f"Association of {traitdir}")
plt.savefig(f"/u/home/z/zhuozshi/project-pasaniuc/ukbb_project/v2/out/{traitdir}/{traitdir}_assoc.png",bbox_inches='tight')
plt.close(fig)
    
    


compression_opts = dict(method='zip',
                        archive_name=f'{traitdir}.csv')  
df.to_csv(f"{out}/{traitdir}/{traitdir}_assoc.zip", index=True, index_label='id',
          compression=compression_opts)  



SNPs = pd.DataFrame(columns=['Index_name', 'L1_BETA', 'L1_SE', 'N', 'P'])



locations = [27*x for x in range(1037)]
locations.append(df.shape[0]-1)

if np.isnan(np.array(df.sort_values(by=['P']).P)[0]) and np.isnan(np.array(df.sort_values(by=['P']).P)[-1]):
    with open(f'/u/home/z/zhuozshi/project-pasaniuc/ukbb_project/v2/out/{traitdir}/{traitdir}_qqplot.txt', 'a') as the_file:
        the_file.write('Wrong')
    with open(f'/u/home/z/zhuozshi/project-pasaniuc/ukbb_project/v2/out/{traitdir}/{traitdir}_histo.txt', 'a') as the_file:
        the_file.write('Wrong')
    print("Wrong")
    


    SNPs = pd.DataFrame(columns=['Index_name', 'L1_BETA', 'L1_SE', 'N', 'P'])



    locations = [27*x for x in range(1037)]
    locations.append(df.shape[0]-1)



else:
    for loc in locations:
        SNPs.loc[len(SNPs.index)] = df.iloc[loc]

    fig, ax = plt.subplots(figsize=(5, 5), dpi=150)
    if traitpos<=106:
        ret = qqplot(list(SNPs["P"]),ax=ax,title=F"Q-Q plot of {traitdir}_{trait}")
    else:
        ret = qqplot(list(SNPs["P"]),ax=ax,title=F"Q-Q plot of {traitdir}")

    fig.savefig(f"/u/home/z/zhuozshi/project-pasaniuc/ukbb_project/v2/out/{traitdir}/{traitdir}_qqplot.png",bbox_inches='tight')



    zscore = np.array(df["L1_BETA"])/np.array(df["L1_SE"])

    plt.figure(figsize=(15,10))
    plt.rcParams.update({'font.size': 22})
    plt.hist(zscore, bins=500,density=True) 
    x=np.linspace(-4,4,10000)
    y=stats.norm.pdf(x, 0, 1)
    plt.plot(x,y,color="orange")
    if traitpos<=106:
        plt.title(f"Histogram of {traitdir}_{trait}")
    else:
        plt.title(f"Histogram of {traitdir}")
    plt.savefig(f"/u/home/z/zhuozshi/project-pasaniuc/ukbb_project/v2/out/{traitdir}/{traitdir}_histo.png",bbox_inches='tight')

if traitpos<=106:
    temp = [f"{traitdir}_{trait}"]
else:
    temp = [traitdir]

gc = admix.data.lambda_gc(np.array(list(SNPs.P)))
temp.append(gc)

t1 = df.P[df.P<5e-5]

temp.append(len(t1))

t2 = df.P[df.P<5e-4]

temp.append(len(t2))

np.save(f"/u/home/z/zhuozshi/project-pasaniuc/ukbb_project/v2/out/{traitdir}/{traitdir}",temp)






