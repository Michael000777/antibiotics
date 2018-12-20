#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

from itertools import product

import sklearn.gaussian_process as sklgp

import platereaderclass as prc

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*",default=[])
parser.add_argument("-D","--designassignment",nargs="*",default=None)
parser.add_argument("-n","--outputgrid",default=50,type=int)
args = parser.parse_args()

data = prc.PlateReaderData(**vars(args))

for fn,title,platedata,designID in data:
    print title
    #print platedata
    
    kernel = sklgp.kernels.Matern(1,(1e-5,1e5),1.5)
    gp = sklgp.GaussianProcessRegressor(kernel = kernel, n_restarts_optimizer = 10)
    
    datagrid0 = data.get_design(designid = designID)[0].flatten()
    datagrid1 = data.get_design(designid = designID)[1].flatten()
    
    design = np.array([[np.log(x),np.log(y)] for x,y in zip(datagrid0,datagrid1)])
    pd = np.array([platedata.flatten()]).T

    #print design
    #print np.shape(design)
    
    gp.fit(design,pd)

    grid0 = np.linspace(np.log(np.min(datagrid0)),np.log(np.max(datagrid0)),num=args.outputgrid)
    grid1 = np.linspace(np.log(np.min(datagrid1)),np.log(np.max(datagrid1)),num=args.outputgrid)
    
    grid  = np.array([[x[0],x[1]] for x in product(grid0,grid1)])
    
    pdpred = gp.predict(grid)
    
    #print pdpred
    #print grid
    
    basename = title.replace(' ','_')
    
    fp = open('GPR_{:s}'.format(basename),'w')
    for x,z in zip(grid,pdpred):
        fp.write('{:.6e} {:.6e} {:.6e}\n'.format(np.exp(x[0]),np.exp(x[1]),z[0]))
    fp.close()
    
    
    #break
    
