#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#from __future__ import print_function

import numpy as np
import pandas as pd
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
    Nm1                   = initialconditions[:,1] - 1
    lB                    = np.log(initialconditions[:,0])
    
    fit,cov               = LMSQ(lB,Nm1)
    fit_with_errors       = uncertainties.correlated_values(fit,cov)
    
    ret                   = dict()
    ret['NB_tau']         = uncertainties.nominal_value(1./fit_with_errors[1])
    ret['NB_tau_stddev']  = uncertainties.std_dev(1./fit_with_errors[1])
    ret['NB_sMIC']        = uncertainties.nominal_value(uexp(-fit_with_errors[0]/fit_with_errors[1]))
    ret['NB_sMIC_stddev'] = uncertainties.std_dev(uexp(-fit_with_errors[0]/fit_with_errors[1]))
    
    if Rsquared:
        residuals         = Nm1 - fit[0] - fit[1] * lB
        ss_res            = np.sum(residuals**2)
        ss_tot            = np.sum((Nm1 - np.mean(Nm1))**2)
        ret['NB_R2']      = 1 - ss_res/ss_tot

    return ret


def estimate_Tau_sMIC_linearFit_AsFuncOfN(initialconditions,Rsquared = False):
    # theory predicts: N > 1 + tau LOG(B/sMIC)
    Nm1                   = initialconditions[:,1] - 1
    lB                    = np.log(initialconditions[:,0])
    
    fit,cov               = LMSQ(Nm1,lB)
    fit_with_errors       = uncertainties.correlated_values(fit,cov)
    
    ret                   = dict()
    ret['BN_tau']         = uncertainties.nominal_value(fit_with_errors[1])
    ret['BN_tau_stddev']  = uncertainties.std_dev(fit_with_errors[1])
    ret['BN_sMIC']        = uncertainties.nominal_value(uexp(fit_with_errors[0]))
    ret['BN_sMIC_stddev'] = uncertainties.std_dev(uexp(fit_with_errors[0]))
    
    if Rsquared:
        residuals         = lB - fit[0] - fit[1] * Nm1
        ss_res            = np.sum(residuals**2)
        ss_tot            = np.sum((lB - np.mean(lB))**2)
        ret['BN_R2']      = 1 - ss_res/ss_tot
        
    return ret


def estimate_Tau_sMIC_singleParameter(initialconditions,Rsquared = False):
    ret = dict()
    
    # ssMIC as geometric mean of all point with smallest initial population size
    mincelldens     = np.min(initialconditions[:,1])
    smic_list       = [m[0] for m in initialconditions if m[1] == mincelldens]
    ret['SP_sMIC']  = np.power(np.prod(smic_list),1./len(smic_list))
    
    # intermediate values for estimation
    Nm1 = initialconditions[:,1] - 1
    lBM = np.log(initialconditions[:,0]/smic)
    
    ret['SPNB_tau'] = np.sum(lBM/Nm1)
    ret['SPBN_tau'] = np.dot(Nm1,lBM)/(np.dot(Nm1,Nm1))

    if Rsquared:
        residuals       = ret['SPNB_tau'] - lBM/Nm1
        ss_res          = np.sum(residuals**2)
        ss_tot          = np.sum((lBM/Nm1 - np.mean(lBM/Nm1))**2)
        ret['SPNB_R2']  = 1 - ss_res/ss_tot
    
        residuals       = lBM - ret['SPBN_tau'] * Nm1
        ss_res          = np.sum(residuals**2)
        ss_tot          = np.sum((lBM - np.mean(lBM))**2)
        ret['SPBN_R2']  = 1 - ss_res/ss_tot
    
    return ret


# *****************************************************************
# ** main
# *****************************************************************

def main():
    parser = argparse.ArgumentParser()
    
    parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
    parser_io.add_argument("-i", "--infiles",                  nargs = "*")
    parser_io.add_argument("-X", "--BasenameExtension",        default = "",    type=str)
    parser_io.add_argument("-G", "--GnuplotOutput",            default = None,  type = str)
    parser_io.add_argument("-g", "--GnuplotColumns",           default = 3,     type = int)
    parser_io.add_argument("-r", "--GnuplotRange",             default = [1e-3,1e2,1e2,1e8], nargs = 4, type = float)
    parser_io.add_argument("-P", "--GenerateImages",           default = False, action = "store_true")
    parser_io.add_argument("-I", "--PlotInoculumCombinations", default = False, action = "store_true")
    parser_io.add_argument("-T", "--WriteThresholdFiles",      default = False, action = "store_true")
    parser_io.add_argument("-F", "--WriteDataFiles",           default = False, action = "store_true")
    
    parser_alg = parser.add_argument_group(description = "==== Algorithm parameters ====")
    parser_alg.add_argument("-M", "--InferenceMethods",          default = ["NfuncB"], choices = ["NfuncB", "BfuncN", "SingleParam"], nargs = "*")
    parser_alg.add_argument("-t", "--GrowthThreshold",           default = None,  type = float)
    parser_alg.add_argument("-D", "--DesignAssignment",          default = [],    type = int, nargs = "*")
    parser_alg.add_argument("-d", "--GenerateDesign",            default = [6e6,4,6.25,2], nargs = 4, type = int)
    parser_alg.add_argument("-R", "--GaussianProcessRegression", default = None, action = "store_true")
    parser_alg.add_argument("-n", "--GPRGridsize",               default = 24,    type = int)
    parser_alg.add_argument("-K", "--GPRKernellist",             default = ['white','matern'], type = str, nargs = "*")
    
    args = parser.parse_args()
    
    # load all data via the 'PlateReaderData' class
    data = prc.PlateReaderData(**vars(args))
    


    if not args.GnuplotOutput is None:
        gnuplotoutput = prc.GnuplotMSPOutput(datasize = len(data), outfilename = args.GnuplotOutput, **vars(args))
        gnuplotoutput.write_init()
   
    if data.count_design == 0:
        data.generate_design(xstart = args.GenerateDesign[0], xdilution = args.GenerateDesign[1], ystart = args.GenerateDesign[2], ydilution = args.GenerateDesign[3])
    
    threshold = data.EstimateGrowthThreshold(dataID = None) # None indicates *ALL* data
    
    columnlist = np.array(['Title','Filename'])
    if 'NfuncB'      in args.InferenceMethods:  columnlist = np.concatenate([columnlist,['NB_sMIC','NB_sMIC_stddev','NB_tau','NB_tau_stddev','NB_R2']])
    if 'BfuncN'      in args.InferenceMethods:  columnlist = np.concatenate([columnlist,['BN_sMIC','BN_sMIC_stddev','BN_tau','BN_tau_stddev','BN_R2']])
    if 'SingleParam' in args.InferenceMethods:  columnlist = np.concatenate([columnlist,['SP_sMIC','SPNB_tau','SPNB_R2','SPBN_tau','SPBN_R2']])

    results = pd.DataFrame(columns = columnlist)

    print('{:30s}  '.format('# Title') + '  '.join(['{:>14.14s}'.format(c) for c in columnlist[2:]]))

    for i in range(len(data)):
        if args.GaussianProcessRegression:  transitions = data.compute_growth_nogrowth_transition_GPR(i, threshold, gridsize = args.GPRGridsize, kernellist = args.GPRKernellist, SaveGPRSurfaceToFile = args.WriteDataFiles)
        else:                               transitions = data.compute_growth_nogrowth_transition    (i, threshold)
        
        curdata             = dict()
        curdata['Title']    = data.titles[i].replace(' ','_')
        curdata['Filename'] = data.filenames[i]

        if 'NfuncB'      in args.InferenceMethods:  curdata.update(estimate_Tau_sMIC_linearFit_AsFuncOfB(transitions, Rsquared = True))
        if 'BfuncN'      in args.InferenceMethods:  curdata.update(estimate_Tau_sMIC_linearFit_AsFuncOfN(transitions, Rsquared = True))
        if 'SingleParam' in args.InferenceMethods:  curdata.update(estimate_Tau_sMIC_singleParameter(transitions, Rsquared = True))

        results.append(curdata,ignore_index = True)

        # main output
        print('{:30.30s}  '.format(curdata['Title']) + '  '.join(['{:14.6e}'.format(curdata[c]) for c in columnlist[2:]]))
        
        basename = (args.BasenameExtension + data.titles[i]).replace(' ','_')

        if args.WriteThresholdFiles or not args.GnuplotOutput is None:
            np.savetxt(basename + '.threshold',transitions)
        
        if not args.GnuplotOutput is None:
            if args.PlotInoculumCombinations:
                gnuplotoutput.write_plot(i, basename, basename, curdata, inoculum = data.get_design(dataID = i))
            else:
                gnuplotoutput.write_plot(i, basename, basename, curdata)

        if args.GenerateImages:
            prc.PlateImage(data[i], data.titles[i], growththreshold = threshold)
    
        if args.WriteDataFiles:
            data.WriteData(dataID = i, filename = basename + '.data')
            
    return results

    
if __name__ == "__main__":
    main()

