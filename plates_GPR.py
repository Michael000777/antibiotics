#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math
import itertools
import sklearn.gaussian_process as sklgp

import platereaderclass as prc


def add_kernel(kernel,newkernel):
    if not kernel is None:
        if   newkernel == 'CONST':              kernel += sklgp.kernels.ConstantKernel()
        elif newkernel == 'WHITE':              kernel += sklgp.kernels.WhiteKernel()
        elif newkernel == 'MATERN':             kernel += sklgp.kernels.Matern()
        elif newkernel == 'RBF':                kernel += sklgp.kernels.RBF()
        elif newkernel == 'EXPSINESQUARED':     kernel += sklgp.kernels.ExpSineSquared()
        elif newkernel == 'DOTPRODUCT':         kernel += sklgp.kernels.DotProduct()
        elif newkernel == 'RATIONALQUADRATIC':  kernel += sklgp.kernels.RationalQuadratic()
    else:
        if   newkernel == 'CONST':              kernel  = sklgp.kernels.ConstantKernel()
        elif newkernel == 'WHITE':              kernel  = sklgp.kernels.WhiteKernel()
        elif newkernel == 'MATERN':             kernel  = sklgp.kernels.Matern()
        elif newkernel == 'RBF':                kernel  = sklgp.kernels.RBF()
        elif newkernel == 'EXPSINESQUARED':     kernel  = sklgp.kernels.ExpSineSquared()
        elif newkernel == 'DOTPRODUCT':         kernel  = sklgp.kernels.DotProduct()
        elif newkernel == 'RATIONALQUADRATIC':  kernel  = sklgp.kernels.RationalQuadratic()
    return kernel


def generate_kernel(argumets):
    available_kernels = ['CONST', 'WHITE', 'MATERN', 'RBF', 'EXPSINESQUARED', 'DOTPRODUCT', 'RATIONALQUADRATIC']
    kernel = None

    klist = [ku.upper() for ku in args.kernels if ku.upper() in available_kernels]
    for k in klist:
        kernel = add_kernel(kernel,k)
    if not kernel is None:
        return kernel
    else:
        raise ValueError('could not define kernel')


def main():
    parser = argparse.ArgumentParser()
    
    parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
    parser_io.add_argument("-i","--infiles",nargs="*",default=[])
    parser_io.add_argument("-D","--designassignment",nargs="*",default=None)
    parser_io.add_argument("-o","--OutfilePrefix",default="GPR",type=str)
    parser_io.add_argument("-n","--OutputGrid",default=50,type=int)
    
    parser_alg = parser.add_argument_group(description = "==== Algorithm parameters ====")
    parser_alg.add_argument("-K","--Kernels",nargs="*", default = "const white matern")
    parser_alg.add_argument("-E","--ErrorEstimates",default=False,action="store_true")
    parser_alg.add_argument("-N","--RestartsOptimizer",default = 10, type = int)
    
    args = parser.parse_args()

    # initiate platereader object and load all data from provided XLSX files
    data = prc.PlateReaderData(**vars(args))


    # iterate over all sheets loaded from all files
    for dataID in range(len(data)):
        title = data.titles[dataID]
        print(title)
        
        # load cellcount, AB conc (axes) and data
        datagrid0         = data.get_design(dataID = dataID)[0].flatten()
        datagrid1         = data.get_design(dataID = dataID)[1].flatten()
        design            = np.array([[np.log(x),np.log(y)] for x,y in zip(datagrid0,datagrid1)])
        platedata         = np.array([data[dataID].flatten()]).T


        # define kernels for Gaussian Process
        kernel = generate_kernel(args.Kernels)
        if args.error_estimates:
            if 'WHITE' not in [k.upper() for k in args.Kernels]:
                kernel = add_kernel(kernel,'WHITE')


        # initiate GPR
        gp = sklgp.GaussianProcessRegressor(kernel = kernel, n_restarts_optimizer = args.RestartsOptimizer)


        # estimate hyperparamters for kernels
        gp.fit(design,platedata)


        # use GPR to estimate values on (fine) grid
        grid0                = np.linspace(np.log(np.min(datagrid0)),np.log(np.max(datagrid0)),num=args.OutputGrid)
        grid1                = np.linspace(np.log(np.min(datagrid1)),np.log(np.max(datagrid1)),num=args.OutputGrid)
        grid                 = np.array([[x[0],x[1]] for x in itertools.product(grid0,grid1)])
        platedata_prediction = gp.predict(grid)

        
        # output
        if len(args.OutfilePrefix) > 0: outfile = (args.OutfilePrefix.strip(' _') + '_' + title).replace(' ','_')
        else:                           outfile = title.replace(' ','_')
        fp = open(outfile,'w')
        lastx = 0
        for x,z in zip(grid,platedata_prediction):
            if x[0] != lastx:fp.write('\n')
            lastx = x[0]
            fp.write('{:.6e} {:.6e} {:.6e}\n'.format(np.exp(x[0]),np.exp(x[1]),z[0]))
        fp.close()



if __name__ == "__main__":
    main()
    
    
