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


    
    
if args[1]=="-t":
    
    assocs = []
    for file in args[2:-1]:
        assocs.append(pd.read_csv(file,index_col=0)) 
    logger.info(
            f"{len(args)-3} association files read"
        )
    
    
    table = []
    for assoc in assocs:
        for index in assoc.index:
            sep = re.sub(r'\W+', ' ',index).split()
            chrn = re.sub("\D", "", sep[0])
            table.append([int(chrn),int(sep[2])])
            
    
    table = np.array(table)
    nrow = table.shape[0]
    table = np.unique(table[np.lexsort((table[:,1],table[:,0]))],axis=0)
    logger.info(
            f"Final template has {table.shape[0]} lanc segments, before prune: {nrow}"
        )
    
    table = pd.DataFrame(table,columns = ['Chrn','Loc'])
    
    table.to_csv(f"{args[-1]}.csv.gz",index=False)  
    logger.info(
            f"Template saved to {args[-1]}"
        )

if args[1]=="-a":
    
    table = pd.read_csv(args[3]).to_numpy()
    logger.info(
            f"Read template {args[3]} with shape {table.shape}"
        )
    
    assoc = pd.read_csv(args[2],index_col=0)
    logger.info(
            f"Read assoc {args[2]} with shape {assoc.shape}"
        )
    
    formatted = []

    for index in assoc.index:

        sep = re.sub(r'\W+', ' ',index).split()
        chrn = re.sub("\D", "", sep[0])

        formatted.append(f"{chrn}:{sep[2]}")

    assoc.index = formatted

    newid = []
    for t in table:

        newid.append(f"{t[0]}:{t[1]}")


        
    assoc = assoc.reindex(newid)

    logger.info(
        f"Reindexed assoc with shape {assoc.shape}. Now interpolating"
    )
    
    
    assoc = assoc.interpolate(method='bfill', limit_direction='backward', axis=0)

    assoc.to_csv(f"{args[-1]}.csv.gz", index=True)
    
    logger.info(
            f"New assoc file saved to {args[-1]}"
        )
    
    