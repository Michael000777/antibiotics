#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import sys,math

def RungeKutta4(func, xx ,time = 0, step = 1e-3):
  # 4th order Runge-Kutta integration scheme
  k1 = step * func(time,              xx )
  k2 = step * func(time + 0.5 * step, xx + 0.5 * k1)
  k3 = step * func(time + 0.5 * step, xx + 0.5 * k2)
  k4 = step * func(time + step,       xx + k3 )
  return xx + (k1 + 2 * k2 + 2 * k3 + k4)/6.


def main():
    parser = argparse.ArgumentParser()
    
    parser_AB = parser.add_argument_group(description = "==== AB paramters ====")
    parser_AB.add_argument("-k", "--ABkappa",           default = 2, type = float)
    parser_AB.add_argument("-g", "--ABgamma",           default = 2, type = float)
    parser_AB.add_argument("-B", "--ABinitconc",        default = 2, type = float)
    parser_AB.add_argument("-s", "--ABsigma",           default = 1, type = float)
    parser_AB.add_argument("-e", "--ABepsilon",         default = 1, type = float)
    parser_AB.add_argument("-K", "--ABmichaelismenten", default = 1, type = float)
    
    parser_PopDyn = parser.add_argument_group(description = "==== Population dynamics ====")
    parser_PopDyn.add_argument("-N", "--PDinitsize",         default = [500], type = float, nargs = "*")
    parser_PopDyn.add_argument("-S", "--PDinitsubstr",       default = 1e6,   type = float, nargs = "*")
    parser_PopDyn.add_argument("-a", "--PDgrowthrate",       default = [1],   type = float, nargs = "*")
    parser_PopDyn.add_argument("-y", "--PDyield",            default = [1],   type = float, nargs = "*")
    parser_PopDyn.add_argument("-r", "--PDrho",              default = [1],   type = float, nargs = "*")
    parser_PopDyn.add_argument("-p", "--PDsigma",            default = [1],   type = float, nargs = "*"
    parser_PopDyn.add_argument("-v", "--PDvolumeseparation", default = 1e-7,  type = float)
    
    parser_alg = parser.add_argument_group(description = "==== Algorithm parameters ====")
    parser_alg.add_argument("-t", "--ALGintegratorstep", default = 1e-3,  type = float)
    parser_alg.add_argument("-T", "--ALGruntime",        default = 24,    type = float)
    parser_alg.add_argument("-O", "--ALGoutputstep",     default = 100,   type = float)
    parser_alg.add_argument("-o", "--ALGoutfilename",    default = "out", type = str)
    
    args = parser.parse_args()
    
    outdata = pd.DataFrame()
    


    outdata.to_csv(args.ALGoutfilename, sep = ' ', header = False)

if __name__ == "__main__":
    main()

