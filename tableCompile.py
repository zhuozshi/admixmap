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
parser.add_argument("-o", "--output",dest = "output")
parser.add_argument("-n", "--name",dest = "name")
parser.add_argument("-i", "--input",dest = "input",nargs='?', type=str, default="/u/home/z/zhuozshi/project-pasaniuc/ukbb_project/v2/out")
parser.add_argument("-p", "--pheno",dest = "pheno",nargs='?', type=str, default="~/ukbb_project/phenotype-description.csv")


args = parser.parse_args()

table = pd.DataFrame(columns=['Trait', 'lambda_gc', '#hit(threshold 5e-5)', '#hit(threshold 5e-4)'])

phenos = pd.read_csv(args.pheno)


for traitpos in range(1,phenos.shape[0]+1):
    
    try:
        temp = float(phenos["id"][traitpos-1])
        traitdir = re.sub(r'\W+', '_',np.array(phenos["description"])[traitpos-1])
    except:
        traitdir = np.array(phenos["id"])[traitpos-1]
    trait = np.array(phenos["id"])[traitpos-1]
    
    
    temp = np.load(f"{args.input}/{traitdir}/{traitdir}.npy")
    
    
    table.loc[len(table.index)] = temp
    
    
table = table.set_index("Trait")
compression_opts = dict(method='zip',
                        archive_name=f'{args.name}.csv')  
table.to_csv(f'{args.output}/{args.name}.zip', index=True,index_label='Trait',
          compression=compression_opts) 
    
    
    