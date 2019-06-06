#!/usr/bin/env python3


import numpy as np
import pandas as pd
import argparse
import sys,math


class TimeIntegratorDynamics(object):
    def __init__(self,**kwargs):
        
        self.__params = dict(kwargs)
        for key in self.__params:
            if isinstance(self.__params[key],(list,tuple)):
                self.__params[key] = np.array(self.__params[key],dtype=np.float)

        self.__verbose = kwargs.get('verbose',False)
        
        na  = len(self.__params.get('PD_growthrate'))
        ny  = len(self.__params.get('PD_yield'))
        ns  = len(self.__params.get('PD_initsize'))
        nr  = len(self.__params.get('PD_rho'))
        nsE = len(self.__params.get('PD_sigma'))
        nsB = len(self.__params.get('AB_sigma'))
        
        assert na == ny == ns == nr == nsE == nsB
        
        self.__numstrains = na
        
        if self.__verbose:
            print("# ==== cmdline parameters ====")
            for k,v in self.__params.items():
                print("# {:<25s} {}".format(k,v))

        
        # initial conditions and description for dynamics
        # population sizes and substrate
        self.__init = np.concatenate([
                                np.array(self.__params.get('PD_initsize'),dtype=np.float),
                                np.array([self.__params.get('PD_initsubstr')],dtype=np.float)
                                ])
        self.__headers = np.concatenate([
                                np.array(['N{:d}'.format(i) for i in range(self.__numstrains)]),
                                ['S']
                                ])

        self.__ABreduction_saturation = True
        if self.__params['AB_michaelismenten'] is None:
            self.__ABreduction_saturation = False

        # track internal concentrations?
        if self.__params['PD_fastinternaldynamics']:
            self.dynamics = self.dynamics_approximateinternal
            
            # store commonly used parameter combinations
            self.__params['es']   = self.__params['AB_epsilon'] / self.__params['AB_sigma']
            self.__params['rs']   = self.__params['PD_rho'] / self.__params['PD_sigma']
            self.__params['i1es'] = 1./(1.+self.__params['es'])
            
        else:
            # needs additional initial conditions and column names
            self.__init = np.concatenate([
                                self.__init,
                                np.zeros(self.__numstrains),
                                np.zeros(self.__numstrains)
                                ])
            self.__headers = np.concatenate([
                                self.__headers,
                                np.array(['Ein{:d}'.format(i) for i in range(self.__numstrains)]),
                                np.array(['Bin{:d}'.format(i) for i in range(self.__numstrains)])
                                ])
            
            
            self.dynamics = self.dynamics_trackinternal
            
            
        # outside concentrations of enzyme and antibiotiocs
        self.__init = np.concatenate([
                                self.__init,
                                np.zeros(1),
                                np.array([self.__params.get('AB_initconc')],dtype=np.float)
                                ])
        self.__headers = np.concatenate([
                                self.__headers,
                                ['Eout'],
                                ['Bout']
                                ])
            
        self.__data = [self.__init]
        self.__time = 0.


    def growthrateEff(self,Bin,substrate):
        """
        helper routines in dynamics
        effective growthrate incorporating death due to antibiotics
        """
        if substrate <= 0:
            return np.zeros(self.__numstrains)
        else:
            bink = np.power(Bin,self.__params['AB_kappa'])
            return self.__params['PD_growthrate'] * (1. - bink)/(1. + bink/self.__params['AB_gamma'])
    
    
    def ABreduction(self,E,B,internal = False):
        """
        returns dynamics for either Michalis-Menten (?) dynamics of AB reduction, or simple linear dynamics
        """
        if self.__ABreduction_saturation:
            if internal:
                return self.__params['AB_epsilon'] * E / (E + self.__params['AB_michaelismenten']*self.__params['i1es']) * B
            else:
                return self.__params['AB_epsilon'] * E / (E + self.__params['AB_michaelismenten']) * B
        else:
            return self.__params['AB_epsilon'] * E * B
    
    
    def dynamics_approximateinternal(self, x, time):
        """
        assume permeability of membrane is fast enough for internal concentrations to equilibrate
        sigmaB, sigmaE >> 1
        """
        N    = x[0:self.__numstrains]
        S    = x[self.__numstrains]
        Eout = x[1+self.__numstrains]
        Bout = x[2+self.__numstrains]
        Ein  = self.__params['rs'] + Eout

        if self.__ABreduction_saturation:   Bin = Bout / (1. + self.__params['es'] * Ein/(Ein + self.__params['AB_michaelismenten']))
        else:                               Bin = Bout / (1. + self.__params['es'] * Ein)
        
        return np.concatenate([
                        np.array(self.growthrateEff(Bin,S) * N),
                        np.array([-np.sum(self.growthrateEff(np.zeros(self.__numstrains),S)/self.__params.get('PD_yield') * N),
                            self.__params['PD_volumeseparation'] * np.sum(self.__params['PD_rho'] * N),
                            -self.ABreduction(Eout,Bout) + self.__params['PD_volumeseparation'] * np.sum(self.__params['i1es'] * N * self.ABreduction(Ein,Bout, internal = True))])
                        ])
                        
                        
    def dynamics_trackinternal(self, x, time):
        """
        all internal concentrations modelled explicitely
        """
        N    = x[0:self.__numstrains]
        S    = x[self.__numstrains]
        Ein  = x[1+self.__numstrains:1+2*self.__numstrains]
        Bin  = x[1+2*self.__numstrains:1+3*self.__numstrains]
        Eout = x[1+3*self.__numstrains]
        Bout = x[2+3*self.__numstrains]
        
        return np.concatenate([
                        np.array(self.growthrateEff(Bin,S) * N),
                        [-np.sum(self.growthrateEff(np.zeros(self.__numstrains),S)/self.__params.get('PD_yield') * N)],
                        np.array(self.__params.get('PD_rho') - self.__params.get('PD_sigma') * (Ein - Eout)),
                        np.array(self.ABreduction(Ein,Bin) - self.__params.get('AB_sigma') * (Bin - Bout)),
                        [np.sum(self.__params.get('PD_sigma') * N * self.__params.get('PD_volumeseparation') * (Ein - Eout))],
                        [-self.ABreduction(Eout,Bout) - self.__params.get('PD_volumeseparation') *  np.sum(self.__params.get('AB_sigma') * N * (Bout - Bin))]
                    ])


    def RungeKutta4(self, func, xx ,time = 0):
        # 4th order Runge-Kutta integration scheme
        k1 = self.__params.get('ALG_integratorstep',1e-3) * func(xx,            time)
        k2 = self.__params.get('ALG_integratorstep',1e-3) * func(xx + 0.5 * k1, time + 0.5 * self.__params.get('ALG_integratorstep',1e-3))
        k3 = self.__params.get('ALG_integratorstep',1e-3) * func(xx + 0.5 * k2, time + 0.5 * self.__params.get('ALG_integratorstep',1e-3))
        k4 = self.__params.get('ALG_integratorstep',1e-3) * func(xx + k3,       time +       self.__params.get('ALG_integratorstep',1e-3))
        return xx + (k1 + 2. * k2 + 2. * k3 + k4)/6.


    def Step(self):
        x = self.__data[-1]
        for s in range(self.__params.get('ALG_outputstep',20)):
            xnew = self.RungeKutta4(self.dynamics, x)
            xnew[np.where(xnew < 0.)[0]]                     = 0.   # concentrations cannot be smaller than 1
            xnew[np.where(xnew[:self.__numstrains] < 1.)[0]] = 0.   # populations are extinct if less than 1 individual
            x    = xnew[:]
        self.__data  = np.concatenate([self.__data,[x]], axis = 0)
        self.__time += self.__params['ALG_integratorstep'] * self.__params['ALG_outputstep']


    def run(self,runtime = None):
        if self.__verbose and self.__time == 0: print(self.lastdata_str)
        starttime = self.__time
        if runtime is None: runtime = self.__params['ALG_runtime']
        while self.__time - starttime <= runtime:
            self.Step()
            if self.__verbose: print(self.lastdata_str)

            
    def save_data(self):
        t = np.array([np.arange(len(self.__data)) * self.__params.get('ALG_integratorstep',1e-3) * self.__params.get('ALG_outputstep')]).T
        d = pd.DataFrame(data = np.concatenate([t, self.__data],axis=1), columns = np.concatenate([['#time'],self.__headers]))
        d.to_csv(self.__params.get('outfilename','out'), sep=' ', float_format = '%.6e', index = False)


    def __getattr__(self,key):
        if key == "lastdata_str":
            return '{:7.3f} '.format(self.__time) + ' '.join(['{:14.6e}'.format(v) for v in self.__data[-1]])








def main():
    parser = argparse.ArgumentParser()
    
    parser_AB = parser.add_argument_group(description = "==== AB parameters ====")
    parser_AB.add_argument("-k", "--AB_kappa",           default = 2.,   type = float)
    parser_AB.add_argument("-g", "--AB_gamma",           default = 2.,   type = float)
    parser_AB.add_argument("-B", "--AB_initconc",        default = 1.5,  type = float)
    parser_AB.add_argument("-e", "--AB_epsilon",         default = 1.,   type = float)
    parser_AB.add_argument("-s", "--AB_sigma",           default = [1.], type = float, nargs = "*")
    parser_AB.add_argument("-K", "--AB_michaelismenten", default = None, type = float)
    
    parser_PopDyn = parser.add_argument_group(description = "==== Population dynamics ====")
    parser_PopDyn.add_argument("-N", "--PD_initsize",             default = [500], type = float, nargs = "*")
    parser_PopDyn.add_argument("-S", "--PD_initsubstr",           default = 1e6,   type = float)
    parser_PopDyn.add_argument("-a", "--PD_growthrate",           default = [1],   type = float, nargs = "*")
    parser_PopDyn.add_argument("-y", "--PD_yield",                default = [1],   type = float, nargs = "*")
    parser_PopDyn.add_argument("-r", "--PD_rho",                  default = [1],   type = float, nargs = "*")
    parser_PopDyn.add_argument("-p", "--PD_sigma",                default = [1],   type = float, nargs = "*")
    parser_PopDyn.add_argument("-V", "--PD_volumeseparation",     default = 1e-10, type = float)
    parser_PopDyn.add_argument("-F", "--PD_fastinternaldynamics", default = False, action = "store_true")
    
    parser_alg = parser.add_argument_group(description = "==== Algorithm parameters ====")
    parser_alg.add_argument("-t", "--ALG_integratorstep", default = 1e-3,  type = float)
    parser_alg.add_argument("-T", "--ALG_runtime",        default = 24,    type = float)
    parser_alg.add_argument("-O", "--ALG_outputstep",     default = 100,   type = float)
    
    parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
    parser_io.add_argument("-o", "--outfilename",    default = "out", type = str)
    parser_io.add_argument("-v", "--verbose",        default = False, action = 'store_true')
    
    args = parser.parse_args()
    
    
    
    # initialize object and run dynamics, save data to text file
    abdyn = TimeIntegratorDynamics(**vars(args))
    
    abdyn.run()
    abdyn.save_data()



if __name__ == "__main__":
    main()


