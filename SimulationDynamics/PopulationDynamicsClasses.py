#!/usr/bin/env python3


import numpy as np
import pandas as pd
import warnings


class TimeIntegrator(object):
    def __init__(self, **kwargs):
        self.AlgParams         = {'IntegrationStep': kwargs.get('IntegrationStep', 2e-3),
                                  'OutputStep':      kwargs.get('OutputStep', 50),
                                  'Runtime':         kwargs.get('Runtime', 24)}

        self.__StoreRuns       = kwargs.get('StoreRuns', False)

        self.Params            = {}
        self.Restrictions      = []

        self.InitialConditions = {}
        self.TrajectoryNames   = []

        self.Trajectories      = []

        self.Reset()


    def RungeKutta4(self, x):
        # 4th order Runge-Kutta integration scheme
        k1 = self.AlgParams['IntegrationStep'] * self.Dynamics(x,            self.time)
        k2 = self.AlgParams['IntegrationStep'] * self.Dynamics(x + 0.5 * k1, self.time + 0.5 * self.AlgParams['IntegrationStep'])
        k3 = self.AlgParams['IntegrationStep'] * self.Dynamics(x + 0.5 * k2, self.time + 0.5 * self.AlgParams['IntegrationStep'])
        k4 = self.AlgParams['IntegrationStep'] * self.Dynamics(x + k3,       self.time +       self.AlgParams['IntegrationStep'])
        return x + (k1 + 2. * k2 + 2. * k3 + k4)/6.


    def Dynamics(self, x, time):
        # needs to be redefined upon inheritance
        return 0


    def Reset(self, **kwargs):
        # possible parameter or IC update
        for key,value in kwargs.items():
            if key in self.InitialConditions.keys(): self.InitialConditions.update({key:value})
            elif key in self.Params.keys():          self.Params.update({key:value})
            else: warnings.warn('Could not update parameter "{}" to value "{}"'.format(key,value))

        # set everything to starting values
        self.time = 0
        self.x    = np.array([[self.InitialConditions[name] for name in self.TrajectoryNames]])


    def Run(self, runtime = None):
        if runtime is None:
            runtime = self.AlgParams['Runtime']

        i         = 0
        xlast     = self.x[-1]
        starttime = self.time
        while self.time <= runtime + starttime:
            xnew  = self.RungeKutta4(xlast)

            # check for restrictions: concentrations can't go below zero or weird things happen
            for name, value1, value2 in self.Restrictions:
                if xnew[self.TrajectoryNames.index(name)] < value1:
                    xnew[self.TrajectoryNames.index(name)] = value2

            xlast = xnew
            i    += 1

            if i % self.AlgParams['OutputStep'] == 0:
                self.x = np.concatenate([self.x, [xnew]], axis = 0)

            self.time += self.AlgParams['IntegrationStep']

        # prepare pandas dataframe as output
        retdict = {'Time': np.arange(len(self.x[:,0])) * self.AlgParams['IntegrationStep'] * self.AlgParams['OutputStep']}
        retdict.update({name:self.x[:,i] for i,name in enumerate(self.TrajectoryNames)})
        retpd = pd.DataFrame(retdict).set_index('Time')

        if self.__StoreRuns:
            self.Trajectories.append((self.Params.copy(), self.InitialConditions.copy(), retpd.copy(deep=True)))

        return retpd



class EnzymePopulationDynamics(TimeIntegrator):
    def __init__(self, **kwargs):
        super(EnzymePopulationDynamics, self).__init__(**kwargs)

        self.TrajectoryNames   = ['N', 'E', 'B']

        self.InitialConditions = {'N': kwargs.get('N', 1e4),
                                  'E': kwargs.get('E', 0),
                                  'B': kwargs.get('B', 1.5)}

        self.Restrictions      = [('N',1,0), ('E',0,0), ('B',0,0)]

        self.Params            = {'epsilon':    kwargs.get('epsilon',1e-3),
                                  'rho':        kwargs.get('rho',1e-3),
                                  'growthrate': kwargs.get('growthrate',1.),
                                  'kappa':      kwargs.get('kappa',2.),
                                  'gamma':      kwargs.get('gamma',2.)}

        self.Reset()


    def Dynamics(self, x, time):
        bk = np.power(x[2], self.Params['kappa'])
        return np.array([
            self.Params['growthrate'] * (1 - bk)/(1 + bk/self.Params['gamma']) * x[0],
            self.Params['rho'] * x[0],
            np.max([-self.Params['epsilon'] * x[1] * x[2],-x[2]]) # second max option restricts
                                                                  # large changes in integrator step.
                                                                  # enzyme x[1] can explode on exp grow pop x[0]
        ])



class EnzymePopulationDynamicsInternal(TimeIntegrator):
    def __init__(self, **kwargs):
        super(EnzymePopulationDynamicsInternal, self).__init__(**kwargs)

        # for our experiments
        volume_ecoli = 1e-15 # L
        volume_well  = 2e-4  # L
        eta_ecoli    = volume_ecoli/volume_well


        self.TrajectoryNames   = ['N', 'Ein', 'Eout', 'Bin', 'Bout' ]#, 'Bout_Din', 'Bout_Dout']

        self.Restrictions      = [('N',1,0), ('Ein',0,0), ('Eout',0,0), ('Bin',0,0), ('Bout',0,0)]

        self.Params            = {'epsilon':    kwargs.get('epsilon',1e-3),
                                  'rho':        kwargs.get('rho',1e-3),
                                  'growthrate': kwargs.get('growthrate',1.),
                                  'kappa':      kwargs.get('kappa',2.),
                                  'gamma':      kwargs.get('gamma',2.),
                                  'sigmaE':     kwargs.get('sigmaE', 1e-4),
                                  'sigmaB':     kwargs.get('sigmaB', 1e-2),
                                  'eta':        kwargs.get('eta', eta_ecoli)}

        Ein = self.Params['rho']/self.Params['sigmaE']
        Phi = (1. + self.Params['rho'] * self.Params['epsilon'] / self.Params['sigmaE'] / self.Params['sigmaB'])

        self.InitialConditions = {'N':          kwargs.get('N', 1e4),
                                  'Ein':        kwargs.get('Ein', None),
                                  'Eout':       kwargs.get('Eout', 0),
                                  'Bin':        kwargs.get('Bin', None),
                                  'Bout':       kwargs.get('Bout', 1.5)
                                 }
        if self.InitialConditions['Bin'] is None:
            self.InitialConditions['Bin'] = self.InitialConditions['Bout']/Phi
        if self.InitialConditions['Ein'] is None:
            self.InitialConditions['Ein'] = Ein

        self.Reset()


    def Dynamics(self, x, time):
        bk = np.power(x[3], self.Params['kappa'])
        return np.array([
            self.Params['growthrate'] * (1 - bk)/(1 + bk/self.Params['gamma']) * x[0],
            self.Params['rho'] + self.Params['sigmaE'] * (x[2] - x[1]),
            self.Params['eta'] * x[0] * self.Params['sigmaE'] * (x[1] - x[2]),
            -self.Params['epsilon'] * x[1] * x[3] + self.Params['sigmaB'] * (x[4] - x[3]),
            -self.Params['epsilon'] * x[2] * x[4] + self.Params['sigmaB'] * self.Params['eta'] * x[0] * (x[3] - x[4])
        ])


