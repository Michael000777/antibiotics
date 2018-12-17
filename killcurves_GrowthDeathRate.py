#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


def MLSQ(x,y):
    n   = len(x)
    sx  = np.sum(x)
    sy  = np.sum(y)
    sxx = np.dot(x,x)
    sxy = np.dot(x,y)
    syy = np.dot(y,y)
    
    denom  = (n*sxx-sx*sx)
    b      = (n*sxy - sx*sy)/denom
    a      = (sy-b*sx)/n
    estim  = np.array([a,b],dtype=np.float)

    sigma2 = syy + n*a*a + b*b*sxx + 2*a*b*sx - 2*a*sy - 2*b*sxy
    cov    = sigma2 / denom * np.array([[sxx,-sx],[-sx,n]],dtype=np.float)

    return estim,cov

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*")
    parser.add_argument("-t","--timepoints",default=-1,type=int)
    args = parser.parse_args()

    data = list()
    n = 0
    for filename in args.infiles:
        print "# load '{}'".format(filename)
        try:
            data.append(np.genfromtxt(filename))
            n += 1
        except:
            continue
            #raise IOError("could not open file")

    if n > 0:
        gmean = np.power(np.prod(data,axis=0),1./n)
    else:
        raise IOError("could not open any input file")

    for i,conc in enumerate(gmean[:,0]):
        if args.timepoints < 0:
            tp = len(gmean[i,1:])
        else:
            tp = args.timepoints
        time = np.arange(0,10*tp,10)/60.
        celln = gmean[i,1:tp + 2]
        mlsq_data = MLSQ(time,np.log(celln))
        gr    = mlsq_data[0][1]
        grDev = np.sqrt(mlsq_data[1][1,1])
        print '{:7.4f} {:7.4f} {:7.4f}'.format(conc,gr,grDev)

if __name__ == "__main__":
    main()
    
    
