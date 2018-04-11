#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import argparse
import openpyxl
import cairo
import sys,math
import uncertainties

from uncertainties.umath import exp as uexp
from skimage import measure
from scipy.optimize import curve_fit


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


def column_string(n):
    div=n
    string=""
    temp=0
    while div>0:
        module=(div-1)%26
        string=chr(65+module)+string
        div=int((div-module)/26)
    return string


def rgb(color):
    r = int(color[0:2],16)
    g = int(color[2:4],16)
    b = int(color[4:6],16)
    return np.array([r/255.,g/255.,b/255.])



# *****************************************************************
# ** read data
# *****************************************************************

def read_sheet(sheet,x0 = 2, y0 = 2, width = 12, height = 8):
    return np.array([[sheet['{:s}{:d}'.format(column_string(i),j)].value for i in range(x0,x0+width)] for j in range(y0,y0+height)],dtype=np.float)


def read_initial_conditions(data,sheetname,xab = 4, yab = 14, xcells = 4, ycells = 3, width = 12, height = 8):
    if sheetname in data.sheetnames:
        return read_sheet(data[sheetname],xab,yab),read_sheet(data[sheetname],xcells,ycells)
    else:
        raise KeyError




# *****************************************************************
# ** process data
# *****************************************************************

def rescale(g,logscale=False,logmin=-20):
    if logscale:
        r = np.log(g)
        r[r<logmin] = 0
    else:
        r = g[:,:]
    r = (r - np.min(r))/(np.max(r) - np.min(r))
    return r


def compute_growth_nogrowth_transition(data,threshold = .5):
    if np.any(data < 0) or np.any(data > 1):
        data = rescale(data)
    transition = measure.find_contours(data,threshold)
    p = list()
    for cont in transition:
        for point in cont:
            p.append(point)
    return np.vstack(p)

    
def flatten_thresholds(contours,filename):
    p = list()
    for cont in contours:
        for c in cont:
            p.append(c)
    return np.vstack(p)


def convert_transitions_to_numbers(points,xdata,ydata):
    ret = list()
    for point in points:
        x_index = int(np.floor(point[0]))
        x_ratio = point[0] - x_index
        
        y_index = int(np.floor(point[1]))
        y_ratio = point[1] - y_index

        if x_index + 1 < len(xdata) and x_ratio > 0:
            x = np.exp(np.log(xdata[x_index,y_index]) * (1-x_ratio) + np.log(xdata[x_index+1,y_index]) * x_ratio)
        else:
            x = xdata[x_index,y_index]
        
        
        if y_index + 1 < len(ydata) and y_ratio > 0:
            y = np.exp(np.log(ydata[x_index,y_index]) * (1-y_ratio) + np.log(ydata[x_index,y_index+1]) * y_ratio)
        else:
            y = ydata[x_index,y_index]
        
        ret.append(np.array([x,y]))
        
    return np.vstack(ret)


def estimate_Tau_sMIC(initialconditions,ABlambda = 1):
    # theory predicts: N > 1 + lambda/tau LOG(B/sMIC)
    Nm1 = initialconditions[:,1] - 1
    lB  = np.log(initialconditions[:,0])
    
    fit,cov = LMSQ(lB,Nm1)
    
    ret = uncertainties.correlated_values(fit,cov)
    
    u_tau = ABlambda/ret[1]
    u_sMIC = uexp(-ret[0]/ret[1])
    
    tau  = np.array([uncertainties.nominal_value(u_tau), uncertainties.std_dev(u_tau)])
    sMIC = np.array([uncertainties.nominal_value(u_sMIC),uncertainties.std_dev(u_sMIC)])
    
    return tau,sMIC



# *****************************************************************
# ** write data
# *****************************************************************

def write_PNG(data,filename,colors = ['3465a4','ffffff','2e3436','eeeeec'],wellsize = 50, wellradius = 20, linewidth = 3):
    cFull   = rgb(colors[0])
    cEmpty  = rgb(colors[1])
    cBorder = rgb(colors[2])
    cBack   = rgb(colors[3])

    CairoImage = cairo.ImageSurface(cairo.FORMAT_ARGB32,data.shape[1] * wellsize,data.shape[0] * wellsize)
    context    = cairo.Context(CairoImage)

    context.rectangle(0,0,data.shape[1] * wellsize,data.shape[0] * wellsize)
    context.set_source_rgb(cBack[0],cBack[1],cBack[2])
    context.fill_preserve()
    context.new_path()

    context.set_line_width(linewidth)
    context.translate(.5 * wellsize,.5 * wellsize)
    
    for x in range(int(data.shape[1])):
        for y in range(int(data.shape[0])):
            context.new_path()
            context.arc(0,0,wellradius,0,2*math.pi)
            c = cFull * data[y,x] + cEmpty * (1 - data[y,x])
            context.set_source_rgb(c[0],c[1],c[2])
            context.fill_preserve()
            context.set_source_rgb(cBorder[0],cBorder[1],cBorder[2])
            context.stroke_preserve()
            context.translate(0,wellsize)
        context.translate(wellsize,-data.shape[0] * wellsize)
    
    CairoImage.write_to_png(filename)


def write_data(filename,data,xdata,ydata):
    fp = open(filename,"w")
    abconc = xdata[0]
    celldens = ydata[:,0]
    for i,x in enumerate(celldens):
        for j,y in enumerate(abconc):
            fp.write('{} {} {}\n'.format(x,y,data[i,j]))
        fp.write('\n')
    fp.close()


# *****************************************************************
# ** main
# *****************************************************************

def main():
    ignoresheets = ["Plate design"]

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles",          nargs = "*")
    parser.add_argument("-t", "--growthThreshold",  type = float, default = 0.1)
    
    parser.add_argument("-P", "--noImages",         default = False, action = "store_true")
    parser.add_argument("-T", "--noThresholdFiles", default = False, action = "store_true")
    parser.add_argument("-F", "--noDataFiles",      default = False, action = "store_true")
    
    parser.add_argument("-L", "--AB_lambda",        type = float, default = 1)
    parser.add_argument("-l", "--AB_lambdaStdDev" , type = float, default = 0)
    args = parser.parse_args()
    
    # use uncertainties package to compute error propagation
    ABlambda = uncertainties.ufloat(args.AB_lambda,args.AB_lambdaStdDev)
    
    for fn in args.infiles:
        try:
            data = openpyxl.load_workbook(fn)
        except:
            raise IOError("could not open file")
    
        print("# File:            {:s}".format(fn))
        print("# Threshold value: {:f}".format(args.growthThreshold))
        
        # read plate design, ie initial conditions for cell numbers and antibiotic concentrations
        abconc,celldens = read_initial_conditions(data,"Plate design")
        
        # iterate over all sheets (except the one with initial conditions)
        for sheet in [s for s in data if not s.title in ignoresheets]:
            basename = sheet.title.replace(' ','_')
            
            # read data into numpy and rescale
            growth  = read_sheet(sheet)
            growth  = rescale(growth)
            
            
            
            # get transitions from growth to nogrowth and get the appropriate initial conditions
            transitions = compute_growth_nogrowth_transition(growth,threshold = args.growthThreshold)
            initialconditions = convert_transitions_to_numbers(transitions,abconc,celldens)
            
            # estimate the curve between growth and nogrowth with the analytical prediction
            # N0 > 1 + lambda/tau LOG(B0/sMIC)
            tau,smic = estimate_Tau_sMIC(initialconditions)

            # output
            print("{:40s} {:14.6e} {:14.6e} {:14.6e} {:14.6e}".format(basename,tau[0],tau[1],smic[0],smic[1]))
            if not args.noDataFiles:
                write_data(basename + '_data.txt',growth,abconc,celldens)
            if not args.noThresholdFiles:
                np.savetxt(basename + '.txt',initialconditions)
            if not args.noImages:
                write_PNG(growth,basename + '.png')
            


if __name__ == "__main__":
    main()

