#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import argparse
import sys,math
import uncertainties

from uncertainties.umath import exp as uexp

import platereaderclass as prc

# *****************************************************************
# ** helper routines
# *****************************************************************

def LMSQ(x,y):
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


# *****************************************************************
# ** process data
# *****************************************************************

def estimate_Tau_sMIC_linearFit(initialconditions,ABlambda = 1):
    # theory predicts: N > 1 + lambda/tau LOG(B/sMIC)
    Nm1 = initialconditions[:,1] - 1
    lB  = np.log(initialconditions[:,0])
    
    fit,cov = LMSQ(lB,Nm1)
    ret     = uncertainties.correlated_values(fit,cov)
    
    u_tau   = ABlambda/ret[1]
    u_sMIC  = uexp(-ret[0]/ret[1])
    
    tau     = np.array([uncertainties.nominal_value(u_tau), uncertainties.std_dev(u_tau)])
    sMIC    = np.array([uncertainties.nominal_value(u_sMIC),uncertainties.std_dev(u_sMIC)])
    
    return tau,sMIC


def estimate_Tau_sMIC_singleParameter(initialconditions,ABlambda = 1):
    mincelldens = np.min(initialconditions[:,1])
    smic_list   = [m[0] for m in initialconditions if m[1] == mincelldens]
    smic        = np.power(np.prod(smic_list),1./len(smic_list))
    
    lBM  = np.log(initialconditions[:,0]/smic)
    n    = initialconditions[:,1]
    
    tau1 = np.sum(lBM*lBM/n)/ ((np.sum(lBM) - np.sum(lBM/n)) * ABlambda)
    tau2 = np.dot(n-1,lBM)/(np.dot(n-1,n-1) * ABlambda)
    return tau1,tau2,smic


# *****************************************************************
# ** main
# *****************************************************************

def main():
    parser = argparse.ArgumentParser()
    
    parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
    parser_io.add_argument("-i", "--infiles",          nargs = "*")
    parser_io.add_argument("-P", "--Images",           default = False, action = "store_true")
    parser_io.add_argument("-T", "--ThresholdFiles",   default = False, action = "store_true")
    parser_io.add_argument("-F", "--DataFiles",        default = False, action = "store_true")
    parser_io.add_argument("-G", "--GnuplotOutput",    default = False, action = "store_true")
    parser_io.add_argument("-g", "--GnuplotColumns",   type = int, default = 3)
    
    parser_alg.add_argument_group(description = "==== Algorithm parameters ====")
    parser_alg.add_argument("-t", "--growthThreshold",  type = float, default = 0.1)
    parser_alg.add_argument("-D", "--designassignment", type = int,   default = [], nargs = "*")
    parser_alg.add_argument("-L", "--AB_lambda",        type = float, default = 1)
    parser_alg.add_argument("-l", "--AB_lambdaStdDev" , type = float, default = 0)
    
    args = parser.parse_args()
    
    # use uncertainties package to compute error propagation
    ABlambda = uncertainties.ufloat(args.AB_lambda,args.AB_lambdaStdDev)
    
    # load all data via the 'PlateReaderData' class
    data = prc.PlateReaderData(**vars(args))
    
    if args.Images:
        data.Export_All_PNGs()
    
    if args.GnuplotOutput:
        sys.stderr.write("set terminal pngcairo enhanced size 1920,1080\n")
        sys.stderr.write("set output \"inoculumeffect.png\"\n")
        sys.stderr.write("set multiplot\n")
        sys.stderr.write("set border 15 lw 2 lc rgb \"#2e3436\"\n")
        sys.stderr.write("set tics front\n")
        sys.stderr.write("set xra [1e-3:1e2]\n")
        sys.stderr.write("set yra [1e2:1e8]\n")
        sys.stderr.write("set logscale\n")
        ysize = 1./args.GnuplotColumns
        if len(data) % args.GnuplotColumns == 0:
            xsize = 1./(len(data)//args.GnuplotColumns)
        else:
            xsize = 1.(len(data)//args.GnuplotColumns + 1.)
        sys.stderr.write("xsize = {:e}\n".format(xsize))
        sys.stderr.write("ysize = {:e}\n".format(ysize))
        sys.stderr.write("xoffset = 0\n")
        sys.stderr.write("yoffset = 0\n")
        sys.stderr.write("set size {:e},{:e}\n".format(xsize,ysize))
        sys.stderr.write("n0(abconc,taulambda,ssmic) = 1 + log(abconc / ssmic) / taulambda\n")
        sys.stderr.write("set label 1 \"empty\" at graph .5,.05 center front\n")
        sys.stderr.write("unset key\n")
        sys.stderr.write("set samples 1001\n")
        sys.stderr.write("\n")

    
    i = 0
    lastfn = ""
    for fn,title,transitions in data.transitions(threshold = args.growthThreshold):
        
        if fn != lastfn:
            print("# data from '{:s}'".format(fn))
        lastfn = fn
        
        tau1,smic1 = estimate_Tau_sMIC_linearFit(transitions)
        tau2,tau3,smic2 = estimate_Tau_sMIC_singleParameter(transitions)
        
        print("{:40s} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}".format(title,tau1[0],tau1[1],smic1[0],smic1[1],tau2,tau3,smic2))
        
        if args.ThresholdFiles or args.GnuplotOutput:
            np.savetxt(title.replace(' ','_') + '.T{:f}.txt'.format(args.growthThreshold),transitions)
        
        if args.GnuplotOutput:
            sys.stderr.write("set origin xoffset + {:d} * xsize, yoffset + {:d} * ysize\n" . format(i//args.GnuplotColumns,i % args.GnuplotColumns))
            sys.stderr.write("set label 1 \"{:s}\"\n".format(title.replace("\"","'")))
            sys.stderr.write("plot \\\n")
            sys.stderr.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#75507b\",\\\n".format(tau2,smic2))
            sys.stderr.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#f57900\",\\\n".format(tau3,smic2))
            sys.stderr.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#cc0000\",\\\n".format(tau1[0],smic1[0]))
            sys.stderr.write("  \"{:s}\" u 1:2 w p pt 7 ps 2 lc rgb \"#3465a4\"\n".format(title.replace(' ','_') + '.T{:f}.txt'.format(args.growthThreshold)))
            sys.stderr.write("\n")
        i += 1
    
    if args.DataFiles:
        for fn,title,platedata in data:
            np.savetxt(title.replace(' ','_') + '.data.txt',platedata)

if __name__ == "__main__":
    main()

