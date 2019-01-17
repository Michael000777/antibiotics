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

def estimate_Tau_sMIC_linearFit_AsFuncOfB(initialconditions,ABlambda = 1,Rsquared = False):
    # theory predicts: N > 1 + lambda/tau LOG(B/sMIC)
    Nm1 = initialconditions[:,1] - 1
    lB  = np.log(initialconditions[:,0])
    
    fit,cov = LMSQ(lB,Nm1)
    ret     = uncertainties.correlated_values(fit,cov)
    
    u_tau   = ABlambda/ret[1]
    u_sMIC  = uexp(-ret[0]/ret[1])
    
    tau     = np.array([uncertainties.nominal_value(u_tau), uncertainties.std_dev(u_tau)])
    sMIC    = np.array([uncertainties.nominal_value(u_sMIC),uncertainties.std_dev(u_sMIC)])
    
    if not Rsquared:
        return tau,sMIC
    else:
        residuals = Nm1 - fit[0] - fit[1] * lB
        ss_res    = np.sum(residuals**2)
        ss_tot    = np.sum((Nm1 - np.mean(Nm1))**2)
        R2        = 1 - ss_res/ss_tot
        return tau,sMIC,R2


def estimate_Tau_sMIC_linearFit_AsFuncOfN(initialconditions,ABlambda = 1,Rsquared = False):
    # theory predicts: N > 1 + lambda/tau LOG(B/sMIC)
    Nm1 = initialconditions[:,1] - 1
    lB  = np.log(initialconditions[:,0])
    
    fit,cov = LMSQ(Nm1,lB)
    ret     = uncertainties.correlated_values(fit,cov)
    
    u_tau   = ret[1]/ABlambda
    u_sMIC  = uexp(ret[0])
    
    tau     = np.array([uncertainties.nominal_value(u_tau), uncertainties.std_dev(u_tau)])
    sMIC    = np.array([uncertainties.nominal_value(u_sMIC),uncertainties.std_dev(u_sMIC)])
    
    if not Rsquared:
        return tau,sMIC
    else:
        residuals = lB - fit[0] - fit[1] * Nm1
        ss_res    = np.sum(residuals**2)
        ss_tot    = np.sum((lB - np.mean(lB))**2)
        R2        = 1 - ss_res/ss_tot
        return tau,sMIC,R2
    


def estimate_Tau_sMIC_singleParameter(initialconditions,ABlambda = 1,Rsquared = False):
    
    # ssMIC as geometric mean of all point with smallest initial population size
    mincelldens = np.min(initialconditions[:,1])
    smic_list   = [m[0] for m in initialconditions if m[1] == mincelldens]
    smic        = np.power(np.prod(smic_list),1./len(smic_list))
    
    # intermediate values for estimation
    Nm1 = initialconditions[:,1] - 1
    lBM = np.log(initialconditions[:,0]/smic)
    
    
    tau1 = np.sum(lBM/Nm1)/ABlambda
    if Rsquared:
        residuals = tau1 - ABlambda * lBM/Nm1
        ss_res    = np.sum(residuals**2)
        ss_tot    = np.sum((ABlambda * lBM/Nm1 - np.mean(ABlambda * lBM/Nm1))**2)
        R2_1      = 1 - ss_res/ss_tot
    
    tau2 = np.dot(Nm1,lBM)/(np.dot(Nm1,Nm1) * ABlambda)
    if Rsquared:
        residuals = lBM - tau2/ABlambda * Nm1
        ss_res    = np.sum(residuals**2)
        ss_tot    = np.sum((lBM - np.mean(lBM))**2)
        R2_2      = 1 - ss_res/ss_tot
    
    if not Rsquared:
        return tau1,tau2,smic
    else:
        return tau1,R2_1,tau2,R2_2,smic


# *****************************************************************
# ** main
# *****************************************************************

def main():
    parser = argparse.ArgumentParser()
    
    parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
    parser_io.add_argument("-i", "--infiles",              nargs = "*")
    parser_io.add_argument("-P", "--Images",               default = False, action = "store_true")
    parser_io.add_argument("-T", "--ThresholdFiles",       default = False, action = "store_true")
    parser_io.add_argument("-F", "--DataFiles",            default = False, action = "store_true")
    parser_io.add_argument("-G", "--GnuplotOutput",        default = False, action = "store_true")
    parser_io.add_argument("-g", "--GnuplotColumns",       default = 3,     type = int)
    parser_io.add_argument("-p", "--GnuplotImageFileName", default = "inoculumeffect", type = str)
    parser_io.add_argument("-H", "--HistogramOD",          default = False, action = "store_true")
    parser_io.add_argument("-b", "--HistogramBins",        default = 20,    type = int)
    parser_io.add_argument("-B", "--HistogramLogscale",    default = False, action = "store_true")
    parser_io.add_argument("-X", "--BasenameExtension",    default = "",    type=str)
    parser_io.add_argument("-d", "--GenerateDesign",       default = [6e6,4,6.25,2], nargs = 4, type = int)
    
    parser_alg = parser.add_argument_group(description = "==== Algorithm parameters ====")
    parser_alg.add_argument("-t", "--growthThreshold",  default = 0.2,   type = float)
    parser_alg.add_argument("-D", "--designassignment", default = [],    type = int, nargs = "*")
    parser_alg.add_argument("-L", "--AB_lambda",        default = 1,     type = float)
    parser_alg.add_argument("-l", "--AB_lambdaStdDev" , default = 0,     type = float)
    parser_alg.add_argument("-R", "--GaussianProcessRegression", default = False, action = "store_true")
    
    args = parser.parse_args()
    
    # use uncertainties package to compute error propagation
    ABlambda = uncertainties.ufloat(args.AB_lambda,args.AB_lambdaStdDev)
    
    # load all data via the 'PlateReaderData' class
    data = prc.PlateReaderData(**vars(args))
    
    if args.Images:
        data.Export_All_PNGs()
    
    if args.GnuplotOutput:
        sys.stderr.write("set terminal pngcairo enhanced size 1920,1080\n")
        sys.stderr.write("set output \"{:s}.png\"\n".format(args.GnuplotImageFileName))
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

    
    
    if data.count_design() == 0:
        data.generate_design(xstart = args.GenerateDesign[0],xdilution = args.GenerateDesign[1], ystart = args.GenerateDesign[2], ydilution = args.GenerateDesign[3])
    
    i = 0
    lastfn = ""
    for fn,title,transitions in data.transitions(threshold = args.growthThreshold, useGPR = args.GaussianProcessRegression):
        basename = (args.BasenameExtension + title).replace(' ','_')
        if fn != lastfn:
            print("# data from '{:s}'".format(fn))
            print("{:40s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s}".format('#','tau lfB','taudev lfB', 'ssmic lfB', 'ssmicdev lfB', 'R2 lfB','tau lfN','taudev lfN', 'ssmic lfN', 'ssmicdev lfN', 'R2 lfN', 'tau Nmin', 'R2 Nmin', 'tau Bmin', 'R2 Bmin', 'ssmic 1'))
        lastfn = fn
        
        tau1,smic1,Rsq1 = estimate_Tau_sMIC_linearFit_AsFuncOfB(transitions,Rsquared = True)
        tau2,smic2,Rsq2 = estimate_Tau_sMIC_linearFit_AsFuncOfN(transitions,Rsquared = True)
        tau3,Rsq3,tau4,Rsq4,smic3 = estimate_Tau_sMIC_singleParameter(transitions,Rsquared = True)
        
        print("{:40s} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}".format(basename,tau1[0],tau1[1],smic1[0],smic1[1],Rsq1,tau2[0],tau2[1],smic2[0],smic2[1],Rsq2,tau3,Rsq3,tau4,Rsq4,smic3))
        
        if args.ThresholdFiles or args.GnuplotOutput:
            np.savetxt(basename + '.threshold{:f}'.format(args.growthThreshold),transitions)
        
        if args.GnuplotOutput:
            sys.stderr.write("set origin xoffset + {:d} * xsize, yoffset + {:d} * ysize\n" . format(i//args.GnuplotColumns,i % args.GnuplotColumns))
            sys.stderr.write("set label 1 \"{:s}\"\n".format(basename.replace('_','-')))
            sys.stderr.write("plot \\\n")
            sys.stderr.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#4e9a06\",\\\n".format(tau3,smic3))
            sys.stderr.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#8ae234\",\\\n".format(tau4,smic3))
            sys.stderr.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#a40000\",\\\n".format(tau1[0],smic1[0]))
            sys.stderr.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#ef2929\",\\\n".format(tau2[0],smic2[0]))
            sys.stderr.write("  \"{:s}\" u 1:2 w p pt 7 ps 2 lc rgb \"#3465a4\"\n".format(basename + '.threshold{:f}'.format(args.growthThreshold)))
            sys.stderr.write("\n")
        i += 1
    
    if args.DataFiles:
        for dataID,plate in enumerate(data):
            fn,title,platedata = plate
            design = data.get_design(dataID = dataID)
            basename = (args.BasenameExtension + title + '.data').replace(' ','_')
            fp = open(basename,'w')
            for i in range(np.shape(platedata)[0]):
                for j in range(np.shape(platedata)[1]):
                    fp.write("{} {} {}\n".format(design[0][i,j],design[1][i,j],platedata[i,j]))
                fp.write("\n")
            fp.close()
    
    if args.HistogramOD:
        if args.HistogramLogscale:
            h,b = np.histogram(np.log(data.all_values()), bins = args.HistogramBins)
            b = np.exp(b[:-1] + np.diff(b))
        else:
            h,b = np.histogram(data.all_values(), range = (0,1), bins = args.HistogramBins)
            b = b[:-1] + np.diff(b)
        np.savetxt('HistogramOD.txt',np.transpose([b,h]))

if __name__ == "__main__":
    main()

