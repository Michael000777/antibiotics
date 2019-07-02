#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import sys,math

import imp
from sklearn.decomposition import PCA

cmdlineparams = {   'N': 'initialsize',
                    'a': 'growthrate',
                    'y': 'yield',
                    's': 'sigmaB',
                    'p': 'sigmaE',
                    'r': 'rho'}
                    


def extract_parameters(runcmd):
    params   = dict()
    curparam = None
    for value in runcmd.split()[1:]:
        if value[0] == '-':
            if value[1] in cmdlineparams.keys():
                curparam = cmdlineparams[value[1]]
                params[curparam] = list()
            elif value[1] == 'o':
                curparam = 'output'
            else:
                curparam = None
        elif curparam == 'output':
            output = value
        elif not curparam is None:
            params[curparam].append(float(value))
        else:
            curparam = None
    numstrains = len(params['initialsize'])
    try:
        simresults = np.genfromtxt(output)
        params['finalsize'] = simresults[-1,1:1+numstrains]
    except:
        pass
    
    
    
    return pd.DataFrame(params)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--runscript")
    args = parser.parse_args()

    try:
        fp = open(args.runscript,'r')
    except:
        raise IOError('could not open file "{:s}"'.format(args.runscript))
    
    pca = PCA()
    
    for line in fp.readlines():
        if line.strip()[0] != '#' and len(line.strip()) > 0:
            p = extract_parameters(line)
            print(p)
            ft = pca.fit_transform(p)
            f  = pca.fit(p)
            print(ft)
            print(f)
            

if __name__ == "__main__":
    main()
