#!/usr/bin/env python

import numpy as np
from scipy.optimize import curve_fit
import argparse
import sys,math


def sigm(logB,zMIC,kappa,gamma):
    return a0*(1. - np.power(np.exp(logB)/zMIC,kappa))/(1. + np.power(np.exp(logB)/zMIC,kappa)/gamma)


def sigm_infer_a0(logB,zMIC,kappa,gamma,aa0):
    return aa0*(1. - np.power(np.exp(logB)/zMIC,kappa))/(1. + np.power(np.exp(logB)/zMIC,kappa)/gamma)

def sigm_red(logB,logM,kappa,gamma,aa0):
    return aa0*(1. - np.exp(kappa * (logB-logM)))/(1. + np.exp(kappa * (logB - logM))/gamma)


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-m","--maxConc",default=-1,type=float)
parser.add_argument("-M","--maxfev",default=5000,type=int)
parser.add_argument("-a","--fixA0",default=False,action="store_true")
args = parser.parse_args()

b  = np.array([])
gr = np.array([])
for fn in args.infiles:
    try:
        data = np.genfromtxt(fn)
    except:
        continue

    b  = np.concatenate([b,data[:,0]])
    gr = np.concatenate([gr,data[:,1]])

if args.maxConc > 0:
    gr = gr[b < args.maxConc]
    b  =  b[b < args.maxConc]


lB  = np.log(b[b > 0])
grL = gr[b>0]


#global a0

#if args.fixA0:
    #a0 = np.mean(gr[b==0])

    #p0 = np.array([1e-2,2,2])
    #fit,cov = curve_fit(sigm,lB,grL,p0=p0,maxfev = args.maxfev)

#else:
    #p0 = np.array([1e-2,2,2,2])
    #fit,cov = curve_fit(sigm_infer_a0,lB,grL,p0=p0,maxfev = args.maxfev)
    #a0 = fit[3]


p0 = np.array([1e-2,2,2,2])

fit,cov = curve_fit(sigm_red,lB,grL,p0 = p0,maxfev = args.maxfev)
a0 = fit[2]

print fit
print cov

print 'growthrate(B) = {:e} * (1 - (B/{:e})**{:e})/(1 + (B/{:e})**{:e}/{:e})'.format(a0,fit[0],fit[1],fit[0],fit[1],fit[2])


#for a,b in zip(lB,grL):
    #print a,b


