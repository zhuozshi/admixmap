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

def tableCompile(rows,out,name='table')


    table = pd.DataFrame(columns=['Trait', 'lambda_gc', '#hit(threshold 5e-5)', '#hit(threshold 5e-4)'])




    for r in rows:
        table.loc[len(table.index)] = np.load(r)


    table = table.set_index("Trait")
    compression_opts = dict(method='zip',
                            archive_name=f'{name}.csv')  
    table.to_csv(f'{out}.zip', index=True,index_label='Trait',
              compression=compression_opts) 
    
    
    