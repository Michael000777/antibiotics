#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

from itertools import product

import sklearn.gaussian_process as sklgp

import platereaderclass as prc


def generate_kernel(argumets):
    available_kernels = ['CONST', 'WHITE', 'MATERN', 'RBF', 'EXPSINESQUARED', 'DOTPRODUCT', 'RATIONALQUADRATIC']
    def add_kernel(kernel,newkernel):
        if not kernel is None:
            if   k == 'CONST':              kernel += sklgp.kernels.ConstantKernel()
            elif k == 'WHITE':              kernel += sklgp.kernels.WhiteKernel()
            elif k == 'MATERN':             kernel += sklgp.kernels.Matern()
            elif k == 'RBF':                kernel += sklgp.kernels.RBF()
            elif k == 'EXPSINESQUARED':     kernel += sklgp.kernels.ExpSineSquared()
            elif k == 'DOTPRODUCT':         kernel += sklgp.kernels.DotProduct()
            elif k == 'RATIONALQUADRATIC':  kernel += sklgp.kernels.RationalQuadratic()
        else:
            if   k == 'CONST':              kernel  = sklgp.kernels.ConstantKernel()
            elif k == 'WHITE':              kernel  = sklgp.kernels.WhiteKernel()
            elif k == 'MATERN':             kernel  = sklgp.kernels.Matern()
            elif k == 'RBF':                kernel  = sklgp.kernels.RBF()
            elif k == 'EXPSINESQUARED':     kernel  = sklgp.kernels.ExpSineSquared()
            elif k == 'DOTPRODUCT':         kernel  = sklgp.kernels.DotProduct()
            elif k == 'RATIONALQUADRATIC':  kernel  = sklgp.kernels.RationalQuadratic()
        return kernel
    kernel = None
    klist = [ku.upper() for ku in args.kernels if ku.upper() in available_kernels]
    for k in klist:
        kernel = add_kernel(kernel,k)
    if not kernel is None:
        return kernel
    else:
        raise ValueError('could not define kernel')

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*",default=[])
parser.add_argument("-D","--designassignment",nargs="*",default=None)
parser.add_argument("-n","--outputgrid",default=50,type=int)
parser.add_argument("-E","--error_estimates",default=False,action="store_true")
parser.add_argument("-o","--outfileprefix",default="GPR",type=str)
parser.add_argument("-K","--kernels",nargs="*", default = "const white matern")
args = parser.parse_args()

data = prc.PlateReaderData(**vars(args))

for dataID in range(len(data)):
    title = data.title[dataID]
    print title
    
    datagrid0         = data.get_design(dataID = dataID)[0].flatten()
    datagrid1         = data.get_design(dataID = dataID)[1].flatten()
    design            = np.array([[np.log(x),np.log(y)] for x,y in zip(datagrid0,datagrid1)])
    platedata         = np.array([data[dataID].flatten()]).T
    if args.error_estimates:
        errorestimate = data.get_noise_estimates(dataID = dataID).flatten()
        
            

    kernel = generate_kernel(args.kernels)
    if args.error_estimates:    gp = sklgp.GaussianProcessRegressor(kernel = kernel, n_restarts_optimizer = 10, alpha = errorestimate)
    else:                       gp = sklgp.GaussianProcessRegressor(kernel = kernel, n_restarts_optimizer = 10)

    gp.fit(design,platedata)

    grid0                = np.linspace(np.log(np.min(datagrid0)),np.log(np.max(datagrid0)),num=args.outputgrid)
    grid1                = np.linspace(np.log(np.min(datagrid1)),np.log(np.max(datagrid1)),num=args.outputgrid)
    grid                 = np.array([[x[0],x[1]] for x in product(grid0,grid1)])
    platedata_prediction = gp.predict(grid)
    
    if len(args.outfileprefix) > 0: outfile = (args.outfileprefix + '_' + title).replace(' ','_')
    else:                           outfile = title.replace(' ','_')
    
    fp = open(outfile,'w')
    lastx = 0
    for x,z in zip(grid,platedata_prediction):
        if x[0] != lastx:fp.write('\n')
        lastx = x[0]
        fp.write('{:.6e} {:.6e} {:.6e}\n'.format(np.exp(x[0]),np.exp(x[1]),z[0]))
    fp.close()
    
