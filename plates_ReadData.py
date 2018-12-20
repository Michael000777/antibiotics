#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

import platereaderclass as prc

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-o","--outfileprefix",default="",type=str)
parser.add_argument("-E","--error_estimates",default=False,action="store_true")
parser.add_argument("-D","--designassignment",default=None,nargs="*")
args = parser.parse_args()

data = prc.PlateReaderData(**vars(args))

for dataID in range(int(data)):
    title = data.title[dataID]
    print title
    if len(args.outfileprefix) > 0: fn = (args.outfileprefix + '_' + title).replace(' ','_')
    else:                           fn = title.replace(' ','_')
    
    data.write_data_to_file(dataID = dataID,filename = fn,include_error_estimates = args.error_estimates)
