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
parser.add_argument("-o", "--output",dest = "output")
parser.add_argument("-i", "--input",dest = "input")
parser.add_argument("-t", "--trait",dest = "trait")
parser.add_argument("-p", "--pheno",dest = "pheno",nargs='?', type=str, default="~/ukbb_project/phenotype-description.csv")


args = parser.parse_args()
out = args.output
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

    
    
chrCount = 22
dfs = []
for i in range(1,chrCount+1):
    if traitdir!=trait:
        name = f"{traitdir}_{trait}_chr{i}"
    else:
        name = f"{trait}_chr{i}"
    
    
    dfs.append(pd.read_csv(f'{args.input}/{name}_assoc.zip',compression='zip',index_col="id"))
    
df = pd.concat(dfs)

    
    
    
    
chrn = []
for i in range(0,chrCount):
    for j in range(dfs[i].shape[0]):
        chrn.append(i+1)   
fig, ax = plt.subplots( figsize=(5, 3), dpi=150)

admix.plot.manhattan(np.array(df.P), ax=ax, s=2,chrom=np.array(chrn),axh_y=-np.log10(5e-5))
plt.tight_layout()

if trait!=traitdir:
    name = f"{traitdir}_{trait}"
else:
    name = traitdir
plt.title(f"Association of {name}")
plt.savefig(f"{out}/{name}_assoc.png",bbox_inches='tight')
plt.close(fig)
    
if np.isnan(np.array(df.sort_values(by=['P']).P)[0]) and np.isnan(np.array(df.sort_values(by=['P']).P)[-1]):
    with open(f'{out}/{name}_qqplot.txt', 'a') as the_file:
        the_file.write('Wrong')
    with open(f'{out}/{name}_histo.txt', 'a') as the_file:
        the_file.write('Wrong')
    print("Wrong")

else:
    #plot qqplot
    SNPs = df
    fig, ax = plt.subplots(figsize=(5, 5), dpi=150)
    ret = qqplot(list(SNPs["P"]),ax=ax,title=F"Q-Q plot of {name}")
    fig.savefig(f"{out}/{name}_qqplot.png",bbox_inches='tight')
    #fig.show()

    
    #plot histogram
    zscore = np.array(df["L1_BETA"])/np.array(df["L1_SE"])

    plt.figure(figsize=(15,10))
    plt.rcParams.update({'font.size': 22})
    plt.hist(zscore, bins=500,density=True) 
    x=np.linspace(-4,4,10000)
    y=stats.norm.pdf(x, 0, 1)
    plt.plot(x,y,color="orange")
    plt.title(f"Histogram of {name}")

    plt.savefig(f"{out}/{name}_histo.png",bbox_inches='tight')
    plt.close(fig)
    
    
temp = [name]

gc = admix.data.lambda_gc(np.array(list(SNPs.P)))
temp.append(gc)

t1 = df.P[df.P<5e-5]

temp.append(len(t1))

t2 = df.P[df.P<5e-4]

temp.append(len(t2))

np.save(f"{out}/{name}",temp)