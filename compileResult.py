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



def compileResult(out,inputPath,name)

    #parse all the arguments




    chrCount = 22
    dfs = []
    for i in inputPath:

        dfs.append(pd.read_csv(i,compression='zip',index_col="id"))

    df = pd.concat(dfs)





    chrn = []
    for i in range(0,chrCount):
        for j in range(dfs[i].shape[0]):
            chrn.append(i+1)
            
    fig, ax = plt.subplots( figsize=(5, 3), dpi=150)
    admix.plot.manhattan(np.array(df.P), ax=ax, s=2,chrom=np.array(chrn),axh_y=-np.log10(5e-5))
    plt.tight_layout()


    plt.title(f"Association of {name}")
    plt.savefig(f"{out}_assoc.png",bbox_inches='tight')
    plt.close(fig)

    if np.isnan(np.array(df.sort_values(by=['P']).P)[0]) and np.isnan(np.array(df.sort_values(by=['P']).P)[-1]):
        with open(f'{out}_qqplot.txt', 'a') as the_file:
            the_file.write('Wrong')
        with open(f'{out}_histo.txt', 'a') as the_file:
            the_file.write('Wrong')
        print("Wrong")

    else:
        #plot qqplot
        SNPs = df
        fig, ax = plt.subplots(figsize=(5, 5), dpi=150)
        ret = qqplot(list(SNPs["P"]),ax=ax,title=F"Q-Q plot of {name}")
        fig.savefig(f"{out}_qqplot.png",bbox_inches='tight')
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

        plt.savefig(f"{out}_histo.png",bbox_inches='tight')
        plt.close(fig)


    temp = [name]

    gc = admix.data.lambda_gc(np.array(list(SNPs.P)))
    temp.append(gc)

    t1 = df.P[df.P<5e-5]

    temp.append(len(t1))

    t2 = df.P[df.P<5e-4]

    temp.append(len(t2))

    np.save(f"{out}",temp)