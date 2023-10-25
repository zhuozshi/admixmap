##########################################################################################
### FileName: tableComile.py
### Author: Zhuozheng Shi
### Date: 10/25/2023
### Affiliation: UCLA Bogdan Lab
### Description: compile table results from all traits in the directory
### Input: [out] output file path and name for table
###        [inputPath] input file path for row of table in  pd.DataFrame file, detecting
###                    files ending in "_tableRow.npy"
###        [name] name of compressed csv table inside the zip file
### Output: table of pd.DataFrame, "{out}.csv.gz", including lambda_gc and #hits for all
###         traits
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
parser.add_argument("-i", "--inputPath",dest = "inputPath")
parser.add_argument("-o", "--out",dest = "out")

args = parser.parse_args()

row = args.inputPath
out = args.out
logger.info(
        f"Currently compiling table in {row}"
    )



table = pd.DataFrame(columns=['Trait', 'lambda_gc', '#hit(threshold 5e-5)', '#hit(threshold 5e-4)'])


rows = []
for file in os.listdir(row):
    if file.endswith("_tableRow.npy"):
        rows.append(os.path.join(row, file))
logger.info(
        f"{len(rows)} found"
    )


for r in rows:
    table.loc[len(table.index)] = np.load(r)


table = table.set_index("Trait")
 
table.to_csv(f'{out}.csv.gz', index=True,index_label='Trait') 
logger.info(
        f"Saving table to {out}.csv.gz"
    )
    
    