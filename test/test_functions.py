#!/usr/bin/env python
#
#   merges and analyzes per-barcode dataframe files. 
#   outputs normalized barcode matrix. 
#
#

import argparse
import logging
import os
import sys
import traceback

from configparser import ConfigParser

import pandas as pd

gitpath=os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import *
from mapseq.utils import *
#
#  targets  A,B,C
#  injection D, E
# 
#  S =           A    B    C       D    E   
#                5   20   40      60   50   
#                5   20   40     120   50   
#                5   20   40     180   50  
#              -------------------------------------
#  Sum spikes:  15   60  120  360  150
#
#  Select highest within groups:   
#                 120 for targets
#                 360 for injection
#  
#  Create weight array so that lowest weight = 1.0
#  So weights for targets:
#             120/15  120/60  120/120
#               8.0     2.0     1.0           
#  Weights for injections:
#                          360/360   260/150  
#                              1.0     2.4 
#
#  R =         2     0    0       5   17 
#              0    50    5      10    0 
#             20     5   90      40    0 
#
#  N =       16.0   0.0  0.0     5.0  40.8
#             0.0 100.0  5.0    10.0   0.0
#           160.0  10.0 90.0    40.0   0.0 
# 

def test_normalize():
    col_names = ['A','B','C','D','E']
    target_columns = ['A','B','C']
    injection_columns = ['D','E']

    rdata = [[  2 ,  0 ,  0 ,  5  , 17 ],
             [  0 , 50 ,  5 ,  10 ,  0 ],
             [ 20 ,  5 , 90 ,  40 ,  0 ]]
    
    rdf = pd.DataFrame(rdata, columns=col_names)
    logging.info(f'Real matrix: rdf=\n{rdf}\n')
    
    # for spike-ins it doesn't matter how many per row, just column sum
    sdata = [[ 5   , 20 , 40 ,  60 , 50  ], 
             [ 5   , 20 , 40 , 120 , 50  ], 
             [ 5   , 20 , 40 , 180 , 50  ]]
    sdf = pd.DataFrame(sdata,columns=col_names)
    logging.info(f'Spike-in matrix: sdf=\n{sdf}\n')
    

    ndf = normalize_weight_grouped(rdf, sdf, columns = [target_columns, injection_columns])    
    logging.info(f'Calculated output: ndf=\n{ndf}')    


    cdata = [[  16.0 ,   0.0 ,  0.0 ,   5.0,  40.8 ],
             [   0.0 , 100.0 ,  5.0 ,  10.0,   0.0 ],
             [ 160.0 ,  10.0 , 90.0 ,  40.0,   0.0 ]]
    cdf = pd.DataFrame(cdata, columns=col_names)
    cdf = cdf.astype(float)
    logging.info(f'Hand-checked result:\n{cdf}')


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    test_normalize()    
    