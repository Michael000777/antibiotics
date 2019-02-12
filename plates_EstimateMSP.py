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
    parser_io.add_argument("-X", "--BasenameExtension",    default = "",    type=str)
    parser_io.add_argument("-P", "--Images",               default = False, action = "store_true")
    parser_io.add_argument("-T", "--ThresholdFiles",       default = False, action = "store_true")
    parser_io.add_argument("-F", "--DataFiles",            default = False, action = "store_true")
    parser_io.add_argument("-G", "--GnuplotOutput",        default = None,  type = str)
    parser_io.add_argument("-g", "--GnuplotColumns",       default = 3,     type = int)
    parser_io.add_argument("-d", "--GenerateDesign",       default = [6e6,4,6.25,2], nargs = 4, type = int)
    
    parser_alg = parser.add_argument_group(description = "==== Algorithm parameters ====")
    parser_alg.add_argument("-t", "--growthThreshold",  default = None,  type = float)
    parser_alg.add_argument("-D", "--designassignment", default = [],    type = int, nargs = "*")
    parser_alg.add_argument("-L", "--AB_lambda",        default = 1,     type = float)
    parser_alg.add_argument("-l", "--AB_lambdaStdDev" , default = 0,     type = float)
    parser_alg.add_argument("-R", "--GaussianProcessRegression", default = False, action = "store_true")
    
    args = parser.parse_args()
    
    # use uncertainties package to compute error propagation
    ABlambda = uncertainties.ufloat(args.AB_lambda,args.AB_lambdaStdDev)
    
    # load all data via the 'PlateReaderData' class
    data = prc.PlateReaderData(**vars(args))
    


    if not args.GnuplotOutput is None:
        fGP = open(args.GnuplotOutput, 'w')
        fGP.write("set terminal pngcairo enhanced size 1920,1080\n")
        fGP.write("set output \"{:s}.png\"\n".format(args.GnuplotOutput))
        fGP.write("set multiplot\n")
        fGP.write("set border 15 lw 2 lc rgb \"#2e3436\"\n")
        fGP.write("set tics front\n")
        fGP.write("set xra [1e-3:1e2]\n")
        fGP.write("set yra [1e2:1e8]\n")
        fGP.write("set logscale\n")
        ysize = 1./args.GnuplotColumns
        if len(data) % args.GnuplotColumns == 0:
            xsize = 1./(len(data)//args.GnuplotColumns)
        else:
            xsize = 1.(len(data)//args.GnuplotColumns + 1.)
        fGP.write("xsize = {:e}\n".format(xsize))
        fGP.write("ysize = {:e}\n".format(ysize))
        fGP.write("xoffset = 0\n")
        fGP.write("yoffset = 0\n")
        fGP.write("set size {:e},{:e}\n".format(xsize,ysize))
        fGP.write("n0(abconc,taulambda,ssmic) = 1 + log(abconc / ssmic) / taulambda\n")
        fGP.write("set label 1 \"empty\" at graph .5,.05 center front\n")
        fGP.write("unset key\n")
        fGP.write("set samples 1001\n")
        fGP.write("\n")

    
    
    if data.count_design == 0:
        data.generate_design(xstart = args.GenerateDesign[0],xdilution = args.GenerateDesign[1], ystart = args.GenerateDesign[2], ydilution = args.GenerateDesign[3])
    
    lastfn    = ""
    threshold = data.EstimateGrowthThreshold(dataID = None) # None indicates *ALL* data

    for i in range(len(data)):
        if args.GaussianProcessRegression:  transitions = data.compute_growth_nogrowth_transition_GPR(i,threshold)
        else:                               transitions = data.compute_growth_nogrowth_transition(i,threshold)
            
        title    = data.titles[i]
        fn       = data.filenames[i]
        basename = (args.BasenameExtension + title).replace(' ','_')
        if fn != lastfn:
            print("# data from '{:s}'".format(fn))
            print("{:40s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s} {:>14s}".format('#','tau lfB','taudev lfB', 'ssmic lfB', 'ssmicdev lfB', 'R2 lfB','tau lfN','taudev lfN', 'ssmic lfN', 'ssmicdev lfN', 'R2 lfN', 'tau Nmin', 'R2 Nmin', 'tau Bmin', 'R2 Bmin', 'ssmic 1'))
        lastfn = fn
        
        
        tau1,smic1,Rsq1 = estimate_Tau_sMIC_linearFit_AsFuncOfB(transitions,Rsquared = True)
        tau2,smic2,Rsq2 = estimate_Tau_sMIC_linearFit_AsFuncOfN(transitions,Rsquared = True)
        tau3,Rsq3,tau4,Rsq4,smic3 = estimate_Tau_sMIC_singleParameter(transitions,Rsquared = True)
        
        # main output
        print("{:40s} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e} {:14.6e}".format(basename,tau1[0],tau1[1],smic1[0],smic1[1],Rsq1,tau2[0],tau2[1],smic2[0],smic2[1],Rsq2,tau3,Rsq3,tau4,Rsq4,smic3))
        
        
        if args.ThresholdFiles or not args.GnuplotOutput is None:
            np.savetxt(basename + '.threshold',transitions)
        
        if not args.GnuplotOutput is None:
            fGP.write("set origin xoffset + {:d} * xsize, yoffset + {:d} * ysize\n" . format(i//args.GnuplotColumns,i % args.GnuplotColumns))
            fGP.write("set label 1 \"{:s}\"\n".format(basename.replace('_','-')))
            fGP.write("plot \\\n")
            fGP.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#4e9a06\",\\\n".format(tau3,smic3))
            fGP.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#8ae234\",\\\n".format(tau4,smic3))
            fGP.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#a40000\",\\\n".format(tau1[0],smic1[0]))
            fGP.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#ef2929\",\\\n".format(tau2[0],smic2[0]))
            fGP.write("  \"{:s}\" u 1:2 w p pt 7 ps 2 lc rgb \"#3465a4\"\n".format(basename + '.threshold'))
            fGP.write("\n")
        i += 1
    
        if args.Images:
            prc.PlateImage(data[i],data.titles[i])
    
        if args.DataFiles:
            outfilename = basename + '.data'
            data.write_data_to_file(i,outfilename)
    
    if not args.GnuplotOutput is None:
        fGP.close()

if __name__ == "__main__":
    main()

