##########################################################################################
### FileName: commonAssoc.py
### Author: Zhuozheng Shi
### Date: 10/30/2023
### Affiliation: UCLA Bogdan Lab
### Description: Generate template for association files, containing all lanc locations
###              Interpolate single association files using all locations in the template
### Input: [assocFiles] file paths for all association files
###        [output] output file path for template
###        or
###        [assocFile] file path for single association file
###        [tempalte] file path for template
###        [output] output file path for interpolated association file
### Output: Template, which contains lanc ending locations in all association files
###         Interpolated association file using template
##########################################################################################

import numpy as np
from os.path import join
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


logger = structlog.get_logger()
args = sys.argv

# flags
if args[1]=="-t":
    logger.info(
            f"Construcing template using {args[2:-1]}"
        )
elif args[1]=="-a":
    assert len(args)==5, f"Number arguments {len(args)-2} does not equal to 3"
    logger.info(
            f"Interpolating association {args[2]} using template {args[3]}"
        )
else:
    logger.error(
            f"False flag given -- {args[1]}"
        )
    
# generate template    
if args[1]=="-t":
    
    # read all association files
    assocs = []
    for file in args[2:-1]:
        assocs.append(pd.read_csv(file,index_col=0)) 
    logger.info(
            f"{len(args)-3} association files read"
        )
    
    # read in all locations
    table = []
    for assoc in assocs:
        for index in assoc.index:
            sep = re.sub(r'\W+', ' ',index).split()
            chrn = re.sub("\D", "", sep[0])
            table.append([int(chrn),int(sep[2])])
            
    # find unique locations
    table = np.array(table)
    nrow = table.shape[0]
    table = np.unique(table[np.lexsort((table[:,1],table[:,0]))],axis=0)
    logger.info(
            f"Final template has {table.shape[0]} lanc segments, before prune: {nrow}"
        )
    
    # save to file
    table = pd.DataFrame(table,columns = ['Chrn','Loc'])
    table.to_csv(f"{args[-1]}.csv.gz",index=False)  
    logger.info(
            f"Template saved to {args[-1]}"
        )

# interpolate association file using template
if args[1]=="-a":
    
    # read template
    table = pd.read_csv(args[3]).to_numpy()
    logger.info(
            f"Read template {args[3]} with shape {table.shape}"
        )
    
    # read association files
    assoc = pd.read_csv(args[2],index_col=0)
    logger.info(
            f"Read assoc {args[2]} with shape {assoc.shape}"
        )
    
    # first format current index
    formatted = []
    for index in assoc.index:

        sep = re.sub(r'\W+', ' ',index).split()
        chrn = re.sub("\D", "", sep[0])

        formatted.append(f"{chrn}:{sep[2]}")
    assoc.index = formatted

    # read in all index
    newid = []
    for t in table:

        newid.append(f"{t[0]}:{t[1]}")

    # reindex
    assoc = assoc.reindex(newid)
    logger.info(
        f"Reindexed assoc with shape {assoc.shape}. Now interpolating"
    )
    
    # interpolation
    assoc = assoc.interpolate(method='bfill', limit_direction='backward', axis=0)

    # saving
    assoc.to_csv(f"{args[-1]}.csv.gz", index=True)
    logger.info(
            f"New assoc file saved to {args[-1]}"
        )
    
    