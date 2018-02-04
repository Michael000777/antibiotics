#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


def iterate(i,imax,dx,errormax):
    # check it to continue iterating
    r = False
    if imax is None:
        if np.dot(dx,dx) < errormax*errormax:
            r = True
    else:
        if i <= imax:
            r = True
    return r


def main():

    # set command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infile")
    parser.add_argument("-C","--maxconcentration",default=None,type=float)
    parser.add_argument("-A","--NRalpha",default=1,type=float)
    parser.add_argument("-I","--initialguess",type=float,nargs=4,default=None)
    parser.add_argument("-v","--verbose",default=False,action="store_true")
    parser_stop = parser.add_mutually_exclusive_group()
    parser_stop.add_argument("-M","--maxiterations",type=int,default=None)
    parser_stop.add_argument("-E","--maxerror",type=float,default=1e-10)
    args = parser.parse_args()

    # load data
    try:
        data = np.genfromtxt(args.infile)
    except:
        raise IOError("could not open file '{}'".format(args.infile))



    # data preprocessing
    a = data[:,1] # growthrate
    b = data[:,0] # antibiotic concentration

    if not args.maxconcentration is None:
        a = a[b<args.maxconcentration]
        b = b[b<args.maxconcentration]

    a = a[b>0]
    b = b[b>0]
    b = np.log(b)


    # initialize variables
    i  = 0
    dx = np.ones(4)
    x  = np.ones(4) # parameters (A,Gamma,mu,kappa)
    f  = np.zeros(4) # system of equations to solve for f(x) = 0
    j  = np.zeros((4,4)) # jacobian of 

    if not args.initialguess is None:
        x = np.array(args.initialguess,dtype=np.float)

    # main iteration loop
    while iterate(i,args.maxiterations,dx,args.maxerror):

        # define often used expressions
        beta   = np.exp(x[3] * (b - x[2]))
        i1gb   = 1./(1 + x[1] * beta)

        # in the current notation, the growthrate is estimated by
        # a = x[0] * (1-beta) * i1gb

        # evaluate function
        f[0]   = np.sum(a * (1-beta) * i1gb)               - x[0] * np.sum((1-beta)*(1-beta)*i1gb*i1gb)
        f[1]   = np.sum(a * (1-beta) * beta * i1gb * i1gb) - x[0] * np.sum((1-beta)*(1-beta)*beta*i1gb*i1gb*i1gb)
        f[2]   = np.sum(a * beta * i1gb * i1gb)            - x[0] * np.sum((1-beta)*beta*i1gb*i1gb*i1gb)
        f[3]   = np.sum(a * beta * (b-x[2]) * i1gb * i1gb) - x[0] * np.sum((1-beta)*beta*(b-x[2])*i1gb*i1gb*i1gb)


        # evaluate jacobian
        j[0,0] = -np.sum((1-beta)*(1-beta)*i1gb*i1gb)
        j[0,1] = -np.sum(a*(1-beta)*beta*i1gb*i1gb) + 2 * x[0]*np.sum((1-beta)*(1-beta)*beta*i1gb*i1gb*i1gb)
        j[0,2] = (1+x[1])*x[3]*np.sum(a*beta*i1gb*i1gb) - 2*x[0]*(1+x[1])*x[3]*np.sum((1-beta)*beta*i1gb*i1gb*i1gb)
        j[0,3] = -(1+x[1])*np.sum(a*(b-x[2])*beta*i1gb*i1gb) + 2*x[0]*(1+x[1])*np.sum((1-beta)*(b-x[2])*beta*i1gb*i1gb*i1gb)
        
        j[1,0] = -np.sum((1-beta)*(1-beta)*beta*i1gb*i1gb*i1gb)
        j[1,1] = -2*np.sum(a*(1-beta)*beta*beta*i1gb*i1gb*i1gb)-3*x[0]*np.sum((1-beta)*(1-beta)*beta*beta*i1gb*i1gb*i1gb*i1gb)
        j[1,2] = -x[3]*np.sum(a*(1+(2+x[1])*beta)*beta*i1gb*i1gb*i1gb)+x[0]*x[3]*np.sum((1-beta)*(1-(3+2*x[1])*beta)*beta*i1gb*i1gb*i1gb*i1gb)
        j[1,3] = np.sum(a*(1+(2+x[1])*beta)*(b-x[2])*beta*i1gb*i1gb*i1gb)-x[0]*np.sum((1-beta)*(1-(3+2*x[1])*beta)*(b-x[2])*beta*i1gb*i1gb*i1gb*i1gb)

        j[2,0] = -np.sum((1-beta)*beta*i1gb*i1gb*i1gb)
        j[2,1] = -2*np.sum(a*beta*beta*i1gb*i1gb*i1gb)+3*x[0]*np.sum((1-beta)*beta*beta*i1gb*i1gb*i1gb*i1gb)
        j[2,2] = -x[3]*np.sum(a*(1-x[1]*beta)*beta*i1gb*i1gb*i1gb)+x[0]*x[3]*np.sum((1-2*(1+x[1])*beta+x[1]*beta*beta)*beta*i1gb*i1gb*i1gb*i1gb)
        j[2,3] = np.sum(a*(1-x[1]*beta)*(b-x[2])*beta*i1gb*i1gb*i1gb)-x[0]*np.sum((1-2*(1+x[1])*beta+x[1]*beta*beta)*(b-x[2])*beta*i1gb*i1gb*i1gb*i1gb)
        
        j[3,0] = -np.sum((b-x[2])*(1-beta)*beta*i1gb*i1gb*i1gb)
        j[3,1] = -2*np.sum(a*(b-x[2])*beta*beta*i1gb*i1gb*i1gb)+3*x[0]*np.sum((1-2*(1+x[1])*beta+x[1]*beta*beta)*(b-x[2])*beta*i1gb*i1gb*i1gb*i1gb)
        j[3,2] = -np.sum((a*(1+x[1]*beta)*beta + x[3]*(b-x[2])*(1-x[1]*beta)*beta)*i1gb*i1gb*i1gb) + x[0]*np.sum(((1-beta)*beta*(1+x[1]*beta)+x[3]*(b-x[2])*(1-2*(1+x[1])*beta+x[1]*beta*beta)*beta)*i1gb*i1gb*i1gb*i1gb)
        j[3,3] = np.sum(a*(b-x[2])*(b-x[2])*(1-x[1]*beta)*beta*i1gb*i1gb*i1gb) - x[0]*np.sum((b-x[2])*(b-x[2])*(1-2*(1+x[1])*beta+x[1]*beta*beta)*beta*i1gb*i1gb*i1gb*i1gb)

        
        # compute and apply iteration step for NR
        dx     = args.NRalpha * np.dot(np.linalg.inv(j),f)
        x     -= dx
        i     += 1

        # output
        if args.verbose:
            print i,x
    
    # final output:
    print x


# run script
if __name__ == "__main__":
    main()
