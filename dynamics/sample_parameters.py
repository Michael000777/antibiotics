#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math


def read_averages_from_file(filename):
    ret = dict()
    try:
        fp = open(args.AverageFile,'r')
        for line in fp.readlines():
            values = line.split()
            ret[values[0]] = float(values[1])
    except:
        pass
    return ret


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-A", "--AverageFile",  default = None, type = str)
    parser.add_argument("-n", "--numstrains",   default = 3,    type = int)
    parser.add_argument("-N", "--inoculumsize", default = 500,  type = float)
    parser.add_argument("-C", "--CoefficientOfVariation", default = 0.2,  type = float)
    
    args = parser.parse_args()

    avg = dict()
    avg['growthrate'] = 1
    avg['yield']      = 1
    avg['rho']        = 1
    avg['sigmaE']     = 1e3
    avg['sigmaB']     = 1e3

    if not args.AverageFile is None:
        avg.update(read_averages_from_file(args.AverageFile))

    inoc = np.random.multinomial(args.inoculumsize,[1./args.numstrains] * args.numstrains)
    
    growthrate = np.random.normal(loc = avg['growthrate'], scale = args.CoefficientOfVariation * avg['growthrate'], size = args.numstrains)
    yields     = np.random.normal(loc = avg['yield'],      scale = args.CoefficientOfVariation * avg['yield'],      size = args.numstrains)
    rho        = np.random.normal(loc = avg['rho'],        scale = args.CoefficientOfVariation * avg['rho'],        size = args.numstrains)
    sigmaE     = np.random.normal(loc = avg['sigmaE'],     scale = args.CoefficientOfVariation * avg['sigmaE'],     size = args.numstrains)
    sigmaB     = np.random.normal(loc = avg['sigmaB'],     scale = args.CoefficientOfVariation * avg['sigmaB'],     size = args.numstrains)

    retstr  = ' -N ' + ' '.join(['{:d}'.format(n) for n in inoc])
    retstr += ' -a ' + ' '.join(['{:e}'.format(a) for a in growthrate])
    retstr += ' -y ' + ' '.join(['{:e}'.format(y) for y in yields])
    retstr += ' -r ' + ' '.join(['{:e}'.format(r) for r in rho])
    retstr += ' -p ' + ' '.join(['{:e}'.format(s) for s in sigmaE])
    retstr += ' -s ' + ' '.join(['{:e}'.format(s) for s in sigmaB])
    
    print(retstr)

if __name__ == "__main__":
    main()
