#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

import platereaderclass as prc

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
args = parser.parse_args()


data = prc.PlateReaderData(**vars(args))

for i in range(data.count_design()):
    print data.get_design(i)

for fn,title,platedata in data:
    print fn,title
    #print platedata
