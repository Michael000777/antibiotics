#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

import platereaderclass as prc

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
args = parser.parse_args()


data = prc.PlateReaderData(**vars(args))

for fn,title,transitions,design in data.transitions(threshold = .1):
    print title
    print transitions

#h,b = np.histogram(np.log(data.all_values()),range = (-5,0),bins = 50)
#b = b[:-1] + np.diff(b)

#for x in zip(b,h):
    #print np.exp(x[0]),x[1]
