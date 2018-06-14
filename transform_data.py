#!/usr/bin/env python

import argparse
import numpy as np
import sys,math


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    parser.add_argument("-t","--timeinterval",type=float,default=10)
    args = parser.parse_args()

    alldata = dict()

    for fn in args.infiles:
        print "# " + fn
        data = np.genfromtxt(fn)
        
        #print data
        #print data[:,0]
        #exit(1)
        
        labels = data[:,0]
        n = data[:,1:].T
        time = np.arange(len(n[:,0])) * args.timeinterval
        
        #print labels
        #print t
        #print n
        
        for i,label in enumerate(labels):
            if not label in alldata.keys():
                alldata[label] = n[:,i]                

    print "# " + " ".join(["{}".format(key) for key in alldata.keys()])
    for i,t in enumerate(time):
        print t,
        for key in alldata.keys():
            print alldata[key][i],
        print



if __name__ == "__main__":
    main()
    
