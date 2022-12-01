import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib

class PlateImage(object):
    def __init__(self, data, design, fig = None, ax = None, panelID = None, **kwargs):

        self.__data             = data
        self.__design           = design
        self.panelID            = panelID

        self.platesize          = self.__data.shape

        self.__real_coordinates = kwargs.get('RealCoordinates', False)
        self.__threshold        = kwargs.get('Threshold', 0)


        self.design =   { 'Nmax': np.max(design[1]),
                          'Nmin': np.min(design[1]),
                          'Bmax': np.max(design[0]),
                          'Bmin': np.min(design[0]) }

        self.boundaries = self.design.copy()

        if not kwargs.get('colors', None) is None:
            self.colors = kwargs.get('colors')
        else:
            self.colors =   { 'growth':        ( 65,101,164),
                              'death':         (255,255,255),
                              'border_growth': '#555753',
                              'border_death':  '#a40000',
                              'background':    '#eeeeec'}

        if not kwargs.get('sizes', None) is None:
            self.sizes = kwargs.get('sizes')
        else:
            self.sizes =    { 'xfigsize':   8,
                              'yfigsize':   5.5,
                              'wellradius': 0.4,
                              'wellborder': 3}

        if not kwargs.get('fontsize', None) is None:
            self.fontsize = kwargs.get('fontsize')
        else:
            self.fontsize = { 'axes':    15,
                              'label':   15,
                              'legend':  12,
                              'panelID': 20}

        if not kwargs.get('axeslabels', None) is None:
            self.axeslabels = kwargs.get('axeslabels')
        else:
            self.axeslabels = [r'Initial Antibiotic Concentration $B_0$ $[\mu g/ml]$',
                               r'Inoculum size $N_0$ $[$cells$/ml]$']


        if fig is None:
            self.fig, self.ax = plt.subplots(1,1, figsize = (self.sizes['xfigsize'],self.sizes['yfigsize']))
        else:
            self.fig, self.ax = fig, ax

        # rescale data to [0,1]
        if np.max(self.__data) - np.min(self.__data) > 0:
            self.data_rescale, self.threshold_rescale = self.rescale(self.__data, self.__threshold)
        else:
            self.data_rescale      = np.zeros_like(self.__data)
            self.threshold_rescale = 0

        # apply and adjust figure parameters
        self.adjust_ax()

        # plot all wells
        for i in range(self.platesize[0]):
            for j in range(self.platesize[1]):
                self.plot_well([j,self.platesize[0] - 1 - i], self.data_rescale[i,j])

        if not self.panelID is None:
            self.set_panelID(self.panelID)



    # switch between (N,B) values from design to coordinates on 8x12 grid on plate. and vice versa below
    def values2grid(self, values):
        if isinstance(values, (list,tuple,np.ndarray)) and len(values) == 2:
            return np.array([
                (self.platesize[0]-1) * (np.log(values[0]/self.design['Nmin']) / np.log(self.design['Nmax']/self.design['Nmin'])),
                (self.platesize[1]-1) * (np.log(values[1]/self.design['Bmin']) / np.log(self.design['Bmax']/self.design['Bmin']))
                ])
        else:
            return np.array([self.values2grid(v) for v in values])


    def grid2values(self, grid):
        if isinstance(grid, (list,tuple,np.ndarray)) and len(grid) == 2:
            return np.array([
                self.design['Nmin'] * np.power(self.design['Nmax']/self.design['Nmin'], grid[0]/(self.platesize[0]-1.)),
                self.design['Bmin'] * np.power(self.design['Bmax']/self.design['Bmin'], grid[1]/(self.platesize[1]-1.))
                ])
        else:
            return np.array([self.grid2values(v) for g in grid])


    def set_panelID(self, panelID):
        self.ax.annotate(panelID ,[0,1.015], weight = 'bold', fontsize = self.fontsize['panelID'], xycoords = 'axes fraction')


    def plot_grid(self):
        # set boundaries
        br = (1-2*self.sizes['wellradius']) + self.sizes['wellradius'] # boundary radius
        self.ax.set_xlim([-br, self.platesize[1] - 1 + br])
        self.ax.set_ylim([-br, self.platesize[0] - 1 + br])

        self.ax.vlines(x = np.arange(self.platesize[1] + 1), ymin = -br, ymax = self.platesize[0] - 1 + br)
        self.ax.hlines(y = np.arange(self.platesize[0] + 1), xmin = -br, xmax = self.platesize[1] - 1 + br)


    def adjust_ax(self, plate_size = None):
        # unset boundary boxes
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['bottom'].set_visible(False)
        self.ax.spines['left'].set_visible(False)


        # set boundaries
        br = (1 - 2*self.sizes['wellradius']) + self.sizes['wellradius'] # boundary radius
        self.ax.set_xlim([-br, self.platesize[1] - 1 + br])
        self.ax.set_ylim([-br, self.platesize[0] - 1 + br])


        # compute and set positions of ticks on axes
        grid_boundaries = np.array([self.ax.get_ylim(), self.ax.get_xlim()]).T
        self.boundaries['Nmin'], self.boundaries['Bmin'] = self.grid2values(grid_boundaries[0])
        self.boundaries['Nmax'], self.boundaries['Bmax'] = self.grid2values(grid_boundaries[1])

        lNmin           = int(np.ceil (np.log10(self.boundaries['Nmin'])))
        lBmin           = int(np.ceil (np.log10(self.boundaries['Bmin'])))
        lNmax           = int(np.floor(np.log10(self.boundaries['Nmax'])))
        lBmax           = int(np.floor(np.log10(self.boundaries['Bmax'])))

        Nticks          = np.power(10., np.arange(lNmin, lNmax + np.max([1,(lNmax-lNmin)//2]), np.max([1,(lNmax-lNmin)//2])))
        Bticks          = np.power(10., np.arange(lBmin, lBmax + np.max([1,(lBmax-lBmin)//2]), np.max([1,(lBmax-lBmin)//2])))
        Nticks_grid     = self.values2grid([(x,1) for x in Nticks]).T[0]
        Bticks_grid     = self.values2grid([(1,x) for x in Bticks]).T[1]

        self.ax.set_xticks(Bticks_grid)
        self.ax.set_yticks(Nticks_grid)
        self.ax.set_xticklabels([r'$10^{{{:d}}}$'.format(int(np.log10(b))) for b in Bticks], fontsize = self.fontsize['label'])
        self.ax.set_yticklabels([r'$10^{{{:d}}}$'.format(int(np.log10(n))) for n in Nticks], fontsize = self.fontsize['label'])


        # add background
        background = matplotlib.patches.Rectangle([self.ax.get_xlim()[0],self.ax.get_ylim()[0]],np.diff(self.ax.get_xlim()),np.diff(self.ax.get_ylim()), fc = self.colors['background'], zorder = -1)
        self.ax.add_patch(background)


        # set axes labels
        self.ax.set_xlabel(self.axeslabels[0], fontsize = self.fontsize['label'])
        self.ax.set_ylabel(self.axeslabels[1], fontsize = self.fontsize['label'])


    def patch_color(self, value, color1 = None, color2 = None, hexoutput = True):
        value = np.max([0,np.min([1,value])])
        if color1 is None: color1 = self.colors['growth']
        if color2 is None: color2 = self.colors['death']

        if hexoutput:
            return '#{:02X}{:02X}{:02X}'.format(int(color1[0] * value + color2[0] * (1-value)),int(color1[1] * value + color2[1] * (1-value)),int(color1[2] * value + color2[2] * (1-value)))
        else:
            return [int(color1[0] * value + color2[0] * (1-value)),int(color1[1] * value + color2[1] * (1-value)),int(color1[2] * value + color2[2] * (1-value))]


    def plot_well(self, pos, value):
        border_color = self.colors['border_growth']
        if value < self.threshold_rescale:
            border_color = self.colors['border_death']
        circle = matplotlib.patches.Circle(pos, self.sizes['wellradius'], facecolor = self.patch_color(value), edgecolor = border_color, linewidth = self.sizes['wellborder'])
        self.ax.add_patch(circle)


    def rescale(self, platedata, threshold):
        return (platedata - np.min(platedata))/(np.max(platedata) - np.min(platedata)), (threshold - np.min(platedata))/(np.max(platedata)-np.min(platedata))



    def plot_curve(self, popsize = None, abconc = None, **kwargs):
        grid = self.values2grid(np.array([popsize,abconc]).T).T
        self.ax.plot(grid[1], grid[0], **kwargs)


    def plot_MSIcurve(self, tau = None, mueff = None, xi = 1, resolution = 100, **kwargs):
        nlist = np.exp(np.linspace(np.log(self.boundaries['Nmin']*0.9), np.log(self.boundaries['Nmax']/0.9), resolution)) # extend by factor 0.9 to remove bulky end of plotted lines
        blist = np.exp(np.power((nlist-1)/tau,1./xi)) * mueff
        self.plot_curve(popsize = nlist, abconc = blist, **kwargs)






