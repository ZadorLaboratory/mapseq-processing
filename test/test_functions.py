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

def test_normalize():
    col_names = ['A' ,'B','C']
    rdata = [[  2 ,  0 ,  0 ],
             [  0 , 50 ,  5 ],
             [ 20 ,  5 , 90 ]]
    rdf = pd.DataFrame(rdata, columns=col_names)
    logging.info(f'rdf=\n{rdf}')
    
    # for spike-ins it doesn't matter how many per row, just column sum
    sdata = [[ 5  , 20 , 40 ], 
             [5   , 20 , 40 ], 
             [5   , 20 , 40 ]]
    sdf = pd.DataFrame(sdata,columns=col_names)
    logging.info(f'sdf=\n{sdf}')
    
    ndf = normalize_weight(rdf, sdf)    
    logging.info(f'ndf=\n{ndf}')    
    cdata = [[  16  ,  0  ,  0 ],
             [  0   , 100 ,  5 ],
             [  160 ,  10 , 90 ]]
    cdf = pd.DataFrame(cdata, columns=col_names)
    cdf = cdf.astype(float)
    logging.info(f'should equal:\n{cdf}')


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    test_normalize()    
    