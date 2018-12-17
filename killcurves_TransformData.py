#!/usr/bin/env python

import argparse
import numpy as np
import sys,math


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    parser.add_argument("-t","--timeinterval",type=float,default=10)
    parser.add_argument("-v","--verbose",action="store_true",default=False)
    args = parser.parse_args()

    alldata = list()
    alllabels = list()

    for fn in args.infiles:
        if args.verbose:
            print "# " + fn
        data = np.genfromtxt(fn)
        
        labels = data[:,0]
        n = data[:,1:].T
        time = np.arange(len(n[:,0])) * args.timeinterval
        
        for i,label in enumerate(labels):
            alllabels.append(label)
            alldata.append(n[:,i])
    print "# " + " ".join(["{}".format(label) for label in alllabels])
    for i,t in enumerate(time):
        print t," ".join(["{}".format(traj[i]) for traj in alldata])



if __name__ == "__main__":
    main()
    
