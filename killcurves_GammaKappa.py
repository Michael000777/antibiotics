#!/usr/bin/env python

import numpy as np
from scipy.optimize import curve_fit
import argparse
import sys,math


def sigm(logB,logM,kappa,gamma,aa0):
    ekbm  = np.exp(kappa * (logB - logM))
    return aa0*(1. - ekbm)/(1. + ekbm/gamma)


def jac(logB,logM,kappa,gamma,aa0):
    ekbm  = np.exp(kappa * (logB - logM))
    denom = np.power(ekbm + gamma,-2.)
    return np.array([
         aa0 * gamma * (1. + gamma) * kappa * ekbm * denom,
        -aa0 * gamma * (1. + gamma) * (logB - logM) * ekbm * denom,
         aa0 * ekbm * (1. - ekbm) * denom,
        (1. - ekbm) / (1. + ekbm / gamma) ]).T

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*")
    parser.add_argument("-m","--maxConc",default=-1,type=float)
    parser.add_argument("-M","--maxfev",default=5000,type=int)
    parser.add_argument("-I","--starting_guesses",type=float,default=[0,2,2,1],nargs=4)
    parser.add_argument("-w","--useWeights",action="store_true",default=False)
    parser.add_argument("-J","--useJacobian",action="store_true",default=False)
    args = parser.parse_args()

    b  = np.array([])
    gr = np.array([])

    if args.useWeights:     w = np.array([])
    else:                   w = None

    if args.useJacobian:    j = jac
    else:                   j = None

    for fn in args.infiles:
        try:
            data = np.genfromtxt(fn)
        except:
            continue

        b  = np.concatenate([b,data[:,0]])
        gr = np.concatenate([gr,data[:,1]])
        if args.useWeights:
            w = np.concatenate([w,data[:,2]])

    if args.maxConc > 0:
        if args.useWeights:
            w = w[b < args.maxConc]
        gr = gr[b < args.maxConc]
        b  =  b[b < args.maxConc]


    lB  = np.log(b[b > 0])
    grL = gr[b>0]
    if args.useWeights:
        w = w[b>0]

    p0 = np.array(args.starting_guesses)
        
    fit,cov = curve_fit(sigm,lB,grL,p0 = p0,maxfev = args.maxfev, jac = j, sigma = w)
    
    
    # output
    print "Fit parameters"
    print "  lMIC  = {:14.6e}".format(fit[0])
    print "  kappa = {:14.6e}".format(fit[1])
    print "  gamma = {:14.6e}".format(fit[2])
    print "  a0    = {:14.6e}".format(fit[3])
    print
    print "Covariance matrix"
    print cov
    print
    print "Relevant experimental parameters"
    print "  MIC    = {:14.6e}".format(np.exp(fit[0]))
    print "  lambda = {:14.6e}".format(fit[1] * fit[2]/(1.+fit[2]))
    print 
    print "Gnuplot function"
    print "  growthrate(B) = {:e} * (1 - (B/{:e})**{:e})/(1 + (B/{:e})**{:e}/{:e})".format(fit[3],np.exp(fit[0]),fit[1],np.exp(fit[0]),fit[1],fit[2])



if __name__ == "__main__":
    main()
    
    
