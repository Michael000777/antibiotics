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

def estimate_Tau_sMIC_linearFit_AsFuncOfB(initialconditions,Rsquared = False):
    # theory predicts: N > 1 + tau LOG(B/sMIC)
    Nm1 = initialconditions[:,1] - 1
    lB  = np.log(initialconditions[:,0])
    
    fit,cov = LMSQ(lB,Nm1)
    ret     = uncertainties.correlated_values(fit,cov)
    
    u_tau   = 1./ret[1]
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


def estimate_Tau_sMIC_linearFit_AsFuncOfN(initialconditions,Rsquared = False):
    # theory predicts: N > 1 + tau LOG(B/sMIC)
    Nm1 = initialconditions[:,1] - 1
    lB  = np.log(initialconditions[:,0])
    
    fit,cov = LMSQ(Nm1,lB)
    ret     = uncertainties.correlated_values(fit,cov)
    
    u_tau   = ret[1]
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
    


def estimate_Tau_sMIC_singleParameter(initialconditions,Rsquared = False):
    
    # ssMIC as geometric mean of all point with smallest initial population size
    mincelldens = np.min(initialconditions[:,1])
    smic_list   = [m[0] for m in initialconditions if m[1] == mincelldens]
    smic        = np.power(np.prod(smic_list),1./len(smic_list))
    
    # intermediate values for estimation
    Nm1 = initialconditions[:,1] - 1
    lBM = np.log(initialconditions[:,0]/smic)
    
    
    tau1 = np.sum(lBM/Nm1)
    if Rsquared:
        residuals = tau1 - lBM/Nm1
        ss_res    = np.sum(residuals**2)
        ss_tot    = np.sum((lBM/Nm1 - np.mean(lBM/Nm1))**2)
        R2_1      = 1 - ss_res/ss_tot
    
    tau2 = np.dot(Nm1,lBM)/(np.dot(Nm1,Nm1))
    if Rsquared:
        residuals = lBM - tau2 * Nm1
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
    
    parser_alg = parser.add_argument_group(description = "==== Algorithm parameters ====")
    parser_alg.add_argument("-t", "--growthThreshold",           default = None,  type = float)
    parser_alg.add_argument("-D", "--designassignment",          default = [],    type = int, nargs = "*")
    parser_alg.add_argument("-d", "--GenerateDesign",            default = [6e6,4,6.25,2], nargs = 4, type = int)
    parser_alg.add_argument("-R", "--GaussianProcessRegression", default = False, action = "store_true")
    
    args = parser.parse_args()
    
    # load all data via the 'PlateReaderData' class
    data = prc.PlateReaderData(**vars(args))
    


    if not args.GnuplotOutput is None:
        gnuplotoutput = prc.GnuplotMSPOutput(datasize = len(data), outfilename = args.GnuplotOutput, **vars(args))
        gnuplotoutput.write_init()
   
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
            gnuplotoutput.write_plot(i,basename,basename,tau1[0],smic1[0],tau2[0],smic2[0],tau3,tau4,smic3)

        if args.Images:
            prc.PlateImage(data[i],data.titles[i])
    
        if args.DataFiles:
            outfilename = basename + '.data'
            data.write_data_to_file(i,outfilename)

        i += 1
    
if __name__ == "__main__":
    main()

