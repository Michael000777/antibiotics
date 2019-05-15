#!/usr/bin/env python3



import numpy as np
import pandas as pd
import argparse
import sys,math


class TimeIntegratorDynamics(object):
    def __init__(self,**kwargs):
        
        self.__params = dict(kwargs)
        print(self.__params)

        
        na  = len(self.__params.get('PD_growthrate'))
        ny  = len(self.__params.get('PD_yield'))
        ns  = len(self.__params.get('PD_initsize'))
        nr  = len(self.__params.get('PD_rho'))
        nsi = len(self.__params.get('PD_sigma'))
        
        assert na == ny == ns == nr == nsi
        
        self.__numstrains = na
        
        init = np.concatenate([ np.zeros(1),
                                np.array(self.__params.get('PD_initsize'),dtype=np.float),
                                np.array([self.__params.get('PD_initsubstr')],dtype=np.float),
                                np.zeros(self.__numstrains),
                                np.zeros(self.__numstrains),
                                np.zeros(1),
                                np.array([self.__params.get('AB_initconc')]) ])
        
        headers = np.concatenate([  ['time'],
                                    np.array(['N{:d}'.format(i) for i in range(self.__numstrains)]),
                                    ['S'],
                                    np.array(['Ein{:d}'.format(i) for i in range(self.__numstrains)]),
                                    np.array(['Bin{:d}'.format(i) for i in range(self.__numstrains)]),
                                    ['Eout'],
                                    ['Bout']    ])
        
        self.__data = pd.DataFrame([init],columns = headers)


    def RungeKutta4(self,func, xx ,time = 0, step = 1e-3):
        # 4th order Runge-Kutta integration scheme
        k1 = step * func(time,              xx )
        k2 = step * func(time + 0.5 * step, xx + 0.5 * k1)
        k3 = step * func(time + 0.5 * step, xx + 0.5 * k2)
        k4 = step * func(time + step,       xx + k3 )
        return xx + (k1 + 2 * k2 + 2 * k3 + k4)/6.


    def SingleStep(self):
        x    = np.array(self.__data.tail(1)).flatten()
        xnew = self.RungeKutta4(self.dynamics, x)
        self.__data.append(xnew)
    
    
    def dynamics(self, state, time):
        N    = self.__data.loc['N0':'N{:d}'.format(self.__numstrains)].copy()
        S    = self.__data.loc['S']
        Ein  = self.__data.loc['Ein0':'Ein{:d}'.format(self.__numstrains)].copy()
        Bin  = self.__data.loc['Bin0':'Bin{:d}'.format(self.__numstrains)].copy()
        Eout = self.__data.loc['Eout']
        Bout = self.__data.loc['Bout']
        
        return np.concatenate([
                        [self.__params.get('ALG_integratorstep',1e-3)],
                        np.array(self.growthrateEff(Bin) * N),
                        [-np.sum(self.growthrateEff(np.zeros(self.__numstrains))/self.__params.get('yield') * N)],
                        np.array(self.__params.get('PD_rho') + self.__params.get('PD_sigma') * (Ein - Eout)),
                        np.array(self.__params.get('AB_epsilon') * Ein/(Ein + self.__params.get('AB_michaelismenten')) * Bin + self.__params.get('AB_sigma') * (Bin - Bout)),
                        [np.sum(self.__params.get('PD_sigma') * N * self.__params.get('PD_volumeseparation') * (Eout - Ein))],
                        [-self.__params.get('AB_epsilon') * Eout/(Eout + self.__params.get('AB_michaelismenten')) + self.__params.get('AB_sigma') * self.__params.get('PD_volumeseparation') *  np.sum(N * (Bout - Bin))]
                    ])
                                   
        


    def run(self,steps = 100):
        for s in range(steps):
            self.SingleStep()
    
    def save_data(self):
        self.__data.to_csv(self.__params.get('ALG_outfilename','out'), sep=' ', index_label = '#')


def main():
    parser = argparse.ArgumentParser()
    
    parser_AB = parser.add_argument_group(description = "==== AB parameters ====")
    parser_AB.add_argument("-k", "--AB_kappa",           default = 2., type = float)
    parser_AB.add_argument("-g", "--AB_gamma",           default = 2., type = float)
    parser_AB.add_argument("-B", "--AB_initconc",        default = 2., type = float)
    parser_AB.add_argument("-s", "--AB_sigma",           default = 1., type = float)
    parser_AB.add_argument("-e", "--AB_epsilon",         default = 1., type = float)
    parser_AB.add_argument("-K", "--AB_michaelismenten", default = 1., type = float)
    
    parser_PopDyn = parser.add_argument_group(description = "==== Population dynamics ====")
    parser_PopDyn.add_argument("-N", "--PD_initsize",         default = [500], type = float, nargs = "*")
    parser_PopDyn.add_argument("-S", "--PD_initsubstr",       default = 1e6,   type = float, nargs = "*")
    parser_PopDyn.add_argument("-a", "--PD_growthrate",       default = [1],   type = float, nargs = "*")
    parser_PopDyn.add_argument("-y", "--PD_yield",            default = [1],   type = float, nargs = "*")
    parser_PopDyn.add_argument("-r", "--PD_rho",              default = [1],   type = float, nargs = "*")
    parser_PopDyn.add_argument("-p", "--PD_sigma",            default = [1],   type = float, nargs = "*")
    parser_PopDyn.add_argument("-v", "--PD_volumeseparation", default = 1e-7,  type = float)
    
    parser_alg = parser.add_argument_group(description = "==== Algorithm parameters ====")
    parser_alg.add_argument("-t", "--ALG_integratorstep", default = 1e-3,  type = float)
    parser_alg.add_argument("-T", "--ALG_runtime",        default = 24,    type = float)
    parser_alg.add_argument("-O", "--ALG_outputstep",     default = 100,   type = float)
    parser_alg.add_argument("-o", "--ALG_outfilename",    default = "out", type = str)
    
    args = parser.parse_args()
    
    abdyn = TimeIntegratorDynamics(**vars(args))
    
    abdyn.run()
    abdyn.save_data()

if __name__ == "__main__":
    main()

