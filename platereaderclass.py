#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import sys,math
import itertools

from skimage import measure

## loaded upon chosing output format below
#import cairo
#import svgwrite


class PlateReaderData(object):
    def __init__(self,**kwargs):

        self.__infilenames = kwargs.get("infiles",[])
        if self.__infilenames is None:
            raise IOError('No input data provided. Use option -i FILENAME1 [FILENAME2 ...]')
        self.__ignoresheetsparameter = kwargs.get("ignoresheets",[])

        self.__designdata = list()
        self.__designtitle = list()

        self.__data = list()
        self.__filenames = list()
        self.__sheetnames = list()
        self.__designassignment = list()
        
        self.__rescale = kwargs.get("DataRescale",False)
        self.__logscale = kwargs.get("DataLogscale",False)
        self.__logmin  = kwargs.get("DataLogscaleMin",-20)
        
        
        # default positions in Excel-file for data and plate design
        # note: indices start at 1 to account for the normal numbering in Excel!
        self.__read_coordinates         = { 'x0':      kwargs.get("xstart",  2),
                                            'y0':      kwargs.get("ystart",  2),
                                            'xwidth':  kwargs.get("xwidth", 12),
                                            'yheight': kwargs.get("yheight", 8) }
        
        self.__read_coordinates_design0 = { 'x0':      kwargs.get("d0_xstart",  4),
                                            'y0':      kwargs.get("d0_ystart", 14),
                                            'xwidth':  kwargs.get("d0_xwidth", 12),
                                            'yheight': kwargs.get("d0_yheight", 8) }
        
        self.__read_coordinates_design1 = { 'x0':      kwargs.get("d1_xstart",  4),
                                            'y0':      kwargs.get("d1_ystart",  3),
                                            'xwidth':  kwargs.get("d1_xwidth", 12),
                                            'yheight': kwargs.get("d1_yheight", 8) }
        
        self.__ignoresheets = ['Plate Design*']
        if len(self.__ignoresheetsparameter) > 0:
            self.__ignoresheets += self.__ignoresheetsparameter
        
        # load all data at __init__()
        if len(self.__infilenames) > 0:
            for fn in self.__infilenames:
                try:
                    tmpExcelFile = pd.ExcelFile(fn)
                except:
                    continue
                
                for designsheet in [s for s in tmpExcelFile.sheet_names if self.ignoresheet(s)]:
                    tmpDesignData0 = self.read_sheet(tmpExcelFile,designsheet,**self.__read_coordinates_design0)
                    tmpDesignData1 = self.read_sheet(tmpExcelFile,designsheet,**self.__read_coordinates_design1)
                    self.__designdata.append((tmpDesignData0,tmpDesignData1))
                    self.__designtitle.append(designsheet)
                i = 0
                for sheet in [s for s in tmpExcelFile.sheet_names if not self.ignoresheet(s)]:
                    tmpSheetData = self.read_sheet(tmpExcelFile, sheet, **self.__read_coordinates)
                    self.__data.append(tmpSheetData)
                    self.__sheetnames.append(sheet)
                    self.__filenames.append(fn)
                    try:
                        design = kwargs.get("DesignAssignment",[])[i]
                        if not 0 <= design < len(self.__designdata):
                            raise KeyError
                        self.__designassignment.append(design)
                    except:
                        self.__designassignment.append(0)
                    i+=1
        else:
            raise IOError('No input data provided. Use option -i FILENAME1 [FILENAME2 ...]')
        
        
    ########################################
    # code for loading and generating data #
    ########################################
    
    def column_string(self,n):
        # needed to convert column numbers to Excel naming scheme (A ... Z AA AB ... AZ BA ...)
        div=n
        string=""
        temp=0
        while div>0:
            module=(div-1)%26
            string=chr(65+module)+string
            div=int((div-module)/26)
        return string


    def ignoresheet(self,sheetname):
        r = False
        for ignore in self.__ignoresheets:
            if len(ignore) > 0:
                if ignore[-1] == '*':
                    if sheetname[:min(len(sheetname),len(ignore)-1)].upper() == ignore[:-1].upper():
                        r = True
                else:
                    if sheetname.upper() == ignore.upper():
                        r = True
            else:
                if sheetname.upper() == ignore.upper():
                    r = True
        return r

    
    def read_sheet(self, excelfile, sheetname, x0 = None, y0 = None, xwidth = None, yheight = None):
        if x0 is None:       x0      = self.__read_coordinates['x0']
        if y0 is None:       y0      = self.__read_coordinates['y0']
        if xwidth is None:   xwidth  = self.__read_coordinates['xwidth']
        if yheight is None:  yheight = self.__read_coordinates['yheight']

        if sheetname in excelfile.sheet_names:
            return np.array(excelfile.parse(sheetname, usecols = '{}:{}'.format(self.column_string(x0),self.column_string(x0 + xwidth - 1)), header = None)[y0-1 : y0 + yheight-1],dtype = np.float)
        else:
            raise KeyError


    def add_design(self,xstart=6e6,xdilution = 4.,ystart=6.25,ydilution = 2,xdirection = 1,ydirection = -1,xsize = 8, ysize = 12, designtitle = 'generated'):
        if ydirection < 0:  ydirection = -1
        else:               ydirection =  1
        if xdirection < 0:  xdirection = -1
        else:               xdirection =  1
        if designtitle == 'generated':
            designtitle += '_X{}_dil{}_Y{}_dil{}'.format(xstart,xdilution,ystart,ydilution)
        design_x = np.array([[xstart * np.power(xdilution,-i) for j in np.arange(ysize)[::ydirection]] for i in np.arange(xsize)[::xdirection]],dtype=np.float)
        design_y = np.array([[ystart * np.power(ydilution,-j) for j in np.arange(ysize)[::ydirection]] for i in np.arange(xsize)[::xdirection]],dtype=np.float)
        
        self.__designdata.append((design_x,design_y))
        self.__designtitle.append(designtitle)


    def get_design(self,designid = 0, designname = None, dataID = None):
        # get the design, depending if either 'designid', 'designname' or 'dataID' is specified
        if not designname is None:
            if designname in self.__designtitle:
                return self.__designdata[self.__designtitle.index(designname)]
        elif not dataID is None:
            return self.__designdata[self.__designassignment[dataID]]
        else:
            return self.__designdata[designid]
    

    ########################
    # math helper routines #
    ########################

    def avg(self,values,geom = True):
        if len(values) > 0:
            if geom:    return np.power(np.product(values),1./len(values))
            else:       return np.mean(values)
        else:
            raise ValueError
    
    
    def rescale(self,g):
        if self.__rescale:
            if self.__logscale:
                r = np.log(g)
                r[r<self.__logmin] = self.__logmin
            else:
                r = g[:,:]
                
            r = (r - np.min(r))/(np.max(r) - np.min(r))
            return r
        else:
            return g

    
    def interpolate_log_log(self,x1,x2,y1,y2,ythreshold):
        if x1 != x2:
            return x1 * np.exp( np.log(ythreshold/y1) * np.log(x2/x1) / np.log(y2/y1) )
        else:
            return x1


    ##################################
    # Export data in various formats #
    ##################################

    def Export_All_PNGs(self):
        for i in range(self.count):
            PlateImage(self.__data[i], outfilename = self.__sheetnames)


    def write_data_to_file(self,dataID,filename = 'out',include_error_estimates = False):
        xlist = self.__designdata[self.__designassignment[dataID]][0]
        ylist = self.__designdata[self.__designassignment[dataID]][1]
        zlist = self.__data[dataID]
        if include_error_estimates:
            elist = self.get_noise_estimates(dataID)
        s = np.shape(zlist)
        
        fp = open(filename,'w')
        
        for i in range(s[0]):
            for j in range(s[1]):
                if include_error_estimates:
                    fp.write('{} {} {} {}\n'.format(xlist[i,j],ylist[i,j],zlist[i,j],elist[i,j]))
                else:
                    fp.write('{} {} {}\n'.format(xlist[i,j],ylist[i,j],zlist[i,j]))
            fp.write('\n')
        
        fp.close()


    ##########################
    # data analysis routines #
    ##########################
    
    def get_noise_estimates(self,dataID):
        # get rough estimate of noise as variance between neighboring wells on plate
        shape = np.shape(self.__data[dataID])
        ne = np.zeros(shape)

        #corners
        ne[0,0]                   = np.std(self.__data[dataID][:2,:2])
        ne[0,shape[1]-1]          = np.std(self.__data[dataID][:2:,-2:])
        ne[shape[0]-1,0]          = np.std(self.__data[dataID][-2:,:2])
        ne[shape[0]-1,shape[1]-1] = np.std(self.__data[dataID][-2:,-2:])

        # edges
        for i in range(1,shape[0]-1):
            ne[i,0]               = np.std(self.__data[dataID][i-1:i+1,:2])
            ne[i,shape[1]-1]      = np.std(self.__data[dataID][i-1:i+1,-2:])
        
        # edges
        for j in range(1,shape[1]-1):
            ne[0,j]               = np.std(self.__data[dataID][:2,j-1:j+1])
            ne[shape[0]-1,j]      = np.std(self.__data[dataID][-2:,j-1:j+1])
            
            # bulk
            for i in range(1,shape[0]-1):
                ne[i,j]           = np.std(self.__data[dataID][i-1:i+1,j-1:j+1])
        
        return ne



    def transitions(self,threshold,useGPR = False):
        if useGPR:
            threshold = self.EstimateGrowthThreshold(dataID = None)
            return [(self.__filenames[i],self.__sheetnames[i],self.compute_growth_nogrowth_transition_GPR(i,threshold)) for i in range(len(self.__data))]
        else:
            return [(self.__filenames[i],self.__sheetnames[i],self.compute_growth_nogrowth_transition(i,threshold)) for i in range(len(self.__data))]


    def compute_growth_nogrowth_transition(self, dataID, threshold = None, geom = True):
        r = list()
        if threshold is None:
            threshold = self.EstimateGrowthThreshold(dataID)
        platesize = np.shape(self.__data[dataID])
        allcont = measure.find_contours(self.__data[dataID],threshold)
        for cont in allcont:
            for pos in cont:
                    # logarithmic interpolation
                    i, j  = int(pos[0]//1), int(pos[1]//1)
                    di,dj = pos[0] - i, pos[1] - j
                    x = self.__designdata[self.__designassignment[dataID]][0][i,j]
                    if j + 1 < platesize[1] and dj > 0:
                        x *= np.power(self.__designdata[self.__designassignment[dataID]][0][i,j+1]/self.__designdata[self.__designassignment[dataID]][0][i,j],dj)
                    
                    y = self.__designdata[self.__designassignment[dataID]][1][i,j]
                    if i + 1 < platesize[0] and di > 0:
                        y *= np.power(self.__designdata[self.__designassignment[dataID]][1][i+1,j]/self.__designdata[self.__designassignment[dataID]][1][i,j],di)
                    
                    r.append(np.array([x,y]))

        if len(r) > 0:
            return np.vstack(r)
        else:
            return None



    def EstimateGrowthThreshold(self, dataID = None, logscale = True):
        # otsu's method
        # described in IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS (1979)
        # usually used to binarize photos into black/white, here we separate the growth/no-growth transition
        # here modified to have uniform distribution, ie in each bin of the histogram we have one value, but spacing between bins is different
        
        x = list()
        if dataID is None:
            # all data is used
            for i in range(len(self)):
                x += list(self.__data[i].flatten())
        elif isinstance(dataID,(list)):
            # pick specific IDs, provided as list
            for i in dataID:
                x += list(self.__data[i].flatten())
        elif isinstance(dataID,int):
            # estimate threshold only from a single plate
            x = list(self.__data[dataID])
        else:
            # something went wrong
            raise TypeError
        x = np.array(x)
        if logscale: x = np.log(x)

        sx = np.sort(x)
        lx = len(x)

        w   = np.arange(lx,dtype=np.float)/lx
        m   = np.array([np.mean(sx[:k]) * w[k] if k > 0 else 0 for k in range(lx)])
        mT  = np.mean(sx)
        
        sB  = np.array([(mT * w[k] - m[k])**2/(w[k]*(1.-w[k])) if w[k]*(1.-w[k]) > 0 else 0 for k in range(lx)])
        idx = np.argmax(sB)
        
        if logscale:
            return np.exp(sx[idx])
        else:
            return sx[idx]
    
    
    def GaussianProcessRegression(self,dataID,kernellist = ['white','matern'], restarts_optimizer = 10, outputgrid = (50,50), AxesLogScale = True):
        
        if not 'sklgp' in sys.modules:  import sklearn.gaussian_process as sklgp

        
        # helper routines to allow arbitrary combinations of kernels.
        # as this is noisy data, 'WHITE' should be among the choices; testing revealed 'RBF' and 'MATERN' both worked reasonably well
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

        def generate_kernel(kernellist):
            available_kernels = ['CONST', 'WHITE', 'MATERN', 'RBF', 'EXPSINESQUARED', 'DOTPRODUCT', 'RATIONALQUADRATIC']
            kernel = None
            klist = [ku.upper() for ku in kernellist if ku.upper() in available_kernels]
            for k in klist:
                kernel = add_kernel(kernel,k)
            if not kernel is None:
                return kernel
            else:
                raise ValueError('could not define kernel')

        
        # main routine of GPR
        
        # reformat input data into correct array sizes
        datagrid0   = self.__designdata[self.__designassignment[dataID]][0].flatten()
        datagrid1   = self.__designdata[self.__designassignment[dataID]][1].flatten()
        if AxesLogScale:
            design  = np.array([[np.log(x),np.log(y)] for x,y in zip(datagrid0,datagrid1)])
        else:
            design  = np.array([[x,y] for x,y in zip(datagrid0,datagrid1)])
        platedata   = np.array([self.__data[dataID].flatten()]).T

        
        # define kernels for Gaussian Process
        kernel = generate_kernel(kernellist)
        
        # initiate GPR
        gp = sklgp.GaussianProcessRegressor(kernel = kernel, n_restarts_optimizer = restarts_optimizer)

        # estimate hyperparamters for kernels
        gp.fit(design,platedata)
        
        # use GPR to estimate values on (fine) grid
        if isinstance(outputgrid,(list,tuple)):
            grid0            = np.linspace(np.log(np.min(datagrid0)),np.log(np.max(datagrid0)),num=outputgrid[0])
            grid1            = np.linspace(np.log(np.min(datagrid1)),np.log(np.max(datagrid1)),num=outputgrid[1])
            outshape         = (outputgrid[0],outputgrid[1])
        elif isinstance(outputgrid,int):
            grid0            = np.linspace(np.log(np.min(datagrid0)),np.log(np.max(datagrid0)),num=outputgrid)
            grid1            = np.linspace(np.log(np.min(datagrid1)),np.log(np.max(datagrid1)),num=outputgrid)
            outshape         = (outputgrid,outputgrid)
        else:
            raise TypeError
        grid                 = np.array([[x[0],x[1]] for x in itertools.product(grid0,grid1)])
        platedata_prediction = gp.predict(grid)
        
        if AxesLogScale:
            return np.exp(grid0), np.exp(grid1), platedata_prediction.reshape(outshape)
        else:
            return grid0, grid1, platedata_prediction.reshape(outshape)
    

    def compute_growth_nogrowth_transition_GPR(self,dataID,threshold):
        grid0,grid1,pdpred = self.GaussianProcessRegression(dataID)
        contours           = measure.find_contours(pdpred,threshold)
        
        finalc    = list()
        for c in contours:
            for i in range(len(c)):
                ix, iy = int(np.floor(c[i,0])), int(np.floor(c[i,1]))
                px, py = c[i,0] - ix, c[i,1] - iy
                try:    cx = (1.-px)*np.log(grid0[ix]) + px*np.log(grid0[ix+1])
                except: cx = np.log(grid0[ix])
                try:    cy = (1.-py)*np.log(grid1[iy]) + py*np.log(grid1[iy+1])
                except: cy = np.log(grid1[iy])
                finalc.append(np.exp([cx,cy]))
        finalc = np.vstack(finalc)
        return finalc
    
    

    ################
    # python stuff #
    ################
    
    def __iter__(self):
        for fn,title,platedata,designassignment in zip(self.__filenames,self.__sheetnames,self.__data,self.__designassignment):
            yield fn,title,self.rescale(platedata),designassignment
    
    
    def __int__(self):
        return len(self.__data)


    def __len__(self):
        return len(self.__data)

    
    def __getattr__(self,key):
        # access internal variables only via __iter__ method!
        if key == "count":
            return len(self.__data)
        elif key == "count_design":
            return len(self.__designdata)
        elif key == "titles":
            return self.__sheetnames
        elif key == "filenames":
            return self.__filenames
        elif key == "all_values":
            return np.concatenate([x.flatten() for x in self.__data])
        elif key == "designassignment":
            return self.__designassignment
    
    
    def __getitem__(self,key):
        return self.__data[key]
























class PlateImage(object):
    def __init__(self, data, outfilename = None,**kwargs):
        self.__data            = data
        self.__outputformat    = kwargs.get('outputformat','svg')
        self.__logscale        = kwargs.get('logscale',False)
        self.__growththreshold = kwargs.get('growththreshold')
        
        self.extract_figure_file_parameters(kwargs)

        if outfilename is None:
            self.__outfilename = 'out.' + self.__outputformat.lower()
        else:
            self.__outfilename = outfilename
            if self.__outfilename[:-4].upper() != '.' + self.__outputformat.upper():
                self.__outfilename += '.' + self.__outputformat
        
        # load appropriate modules only now
        if self.__outputformat.upper() == 'SVG':
            if not 'svgwrite' in sys.modules: import svgwrite
            self.write_SVG()
        elif self.__outputformat.upper() == 'PNG':
            if not 'cairo' in sys.modules: import cairo
            self.write_PNG()
            


    def rescale(self, g, logscale = False):
        if logscale:
            dmin = np.log(np.min(g))
            dmax = np.log(np.max(g))
            return (np.log(g) - dmin)/(dmax - dmin)
        else:
            dmin = np.min(g)
            dmax = np.max(g)
            return (g - dmin)/(dmax - dmin)
        

    def extract_figure_file_parameters(self,kwargs):
        # default values
        self.figureparameters = {   'colors':   ['3465a4','ffffff','eeeeec','2e3436','a40000'],
                                    'wellradius': 20,
                                    'wellsize':50,
                                    'linewidth':3}
        # update default values if part of the argument dictionary
        if 'FigureWellDistance'        in kwargs.keys():    self.figureparameters['wellsize']   = kwargs['FigureWellDistance']
        if 'FigureWellRadius'          in kwargs.keys():    self.figureparameters['wellradius'] = kwargs['FigureWellRadius']
        if 'FigureLinewidth'           in kwargs.keys():    self.figureparameters['linewidth']  = kwargs['FigureLinewidth']
        if 'FigureColorFull'           in kwargs.keys():    self.figureparameters['colors'][0]  = kwargs['FigureColorFull']
        if 'FigureColorEmpty'          in kwargs.keys():    self.figureparameters['colors'][1]  = kwargs['FigureColorEmpty']
        if 'FigureColorBackground'     in kwargs.keys():    self.figureparameters['colors'][2]  = kwargs['FigureColorBackground']
        if 'FigureColorBorder'         in kwargs.keys():    self.figureparameters['colors'][3]  = kwargs['FigureColorBorder']
        if 'FigureColorBorderNoGrowth' in kwargs.keys():    self.figureparameters['colors'][4]  = kwargs['FigureColorBorderNoGrowth']


        
    def rgb(self,color, outformat = 'list'):
        if color is None: return None
        if isinstance(color, str):
            r = int(color[0:2],16)
            g = int(color[2:4],16)
            b = int(color[4:6],16)
        elif isinstance(color,(list,tuple,np.array)):
            if isinstance(color[0],(int,np.int)):
                r = color[0]
                g = color[1]
                b = color[2]
            elif isinstance(color[0],(float,np.float)):
                r = int(255 * color[0])
                g = int(255 * color[1])
                b = int(255 * color[2])
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError
        
        if outformat == 'list':
            return np.array([r/255.,g/255.,b/255.])
        elif outformat == 'xml':
            return 'rgb({:d},{:d},{:d})'.format(r,g,b)


    def interpolate_color(self, rdatavalue, color1, color2, outformat = 'list'):
        ic = rdatavalue * self.rgb(color1, outformat = 'list') + (1-rdatavalue) * self.rgb(color2, outformat = 'list')
        
        if outformat == 'list':
            return ic
        elif outformat == 'xml':
            return 'rgb({:d},{:d},{:d})'.format(int(255 * ic[0]), int(255 * ic[1]), int(255 * ic[2]))


    def write_SVG(self):
        rdata = self.rescale(self.__data, logscale = self.__logscale)
        rthreshold = -1
        if self.__logscale:
            dmin = np.log(np.min(self.__data))
            dmax = np.log(np.max(self.__data))
            if not self.__growththreshold is None:
                rthreshold = (np.log(self.__growththreshold) - dmin)/(dmax - dmin)
        else:
            dmin = np.min(self.__data)
            dmax = np.max(self.__data)
            if not self.__growththreshold is None:
                rthreshold = (self.__growththreshold - dmin)/(dmax - dmin)

        # for some reason 'svgwrite' is not loaded into global namespace with the local import above in __init__
        # need to access it like that:
        img = sys.modules['svgwrite'].Drawing(self.__outfilename, size = (self.figureparameters['wellsize'] * rdata.shape[1], self.figureparameters['wellsize'] * rdata.shape[0]))
        
        if not self.figureparameters['colors'][2] is None:
            img.add(img.rect((0,0), (self.figureparameters['wellsize'] * rdata.shape[1], self.figureparameters['wellsize'] * rdata.shape[0]), fill = self.rgb(self.figureparameters['colors'][2], outformat = 'xml')))
        
        for x in range(rdata.shape[0]):
            for y in range(rdata.shape[1]):
                if rdata[x,y] < rthreshold: bordercolor = self.rgb(self.figureparameters['colors'][4], outformat = 'xml')
                else:                       bordercolor = self.rgb(self.figureparameters['colors'][3], outformat = 'xml')
                fillcolor = self.interpolate_color(rdata[x,y],self.figureparameters['colors'][0], self.figureparameters['colors'][1], outformat = 'xml')
                
                img.add(img.circle( ( (y + .5) * self.figureparameters['wellsize'], (x + .5) * self.figureparameters['wellsize'] ), self.figureparameters['wellradius'], stroke_width = self.figureparameters['linewidth'], fill = fillcolor, stroke = bordercolor))
        
        img.save()


    def write_PNG(self):
        cFull           = self.rgb(self.figureparameters['colors'][0])
        cEmpty          = self.rgb(self.figureparameters['colors'][1])
        cBack           = self.rgb(self.figureparameters['colors'][2])
        cBorder         = self.rgb(self.figureparameters['colors'][3])
        cBorderNoGrowth = self.rgb(self.figureparameters['colors'][4])

        rthreshold = -1
        if self.__logscale:
            dmin = np.log(np.min(self.__data))
            dmax = np.log(np.max(self.__data))
            if not self.__growththreshold is None:
                rthreshold = (np.log(self.__growththreshold) - dmin)/(dmax - dmin)
        else:
            dmin = np.min(self.__data)
            dmax = np.max(self.__data)
            if not self.__growththreshold is None:
                rthreshold = (self.__growththreshold - dmin)/(dmax - dmin)

        CairoImage      = sys.modules['cairo'].ImageSurface(sys.modules['cairo'].FORMAT_ARGB32,self.__data.shape[1] * self.figureparameters['wellsize'],self.__data.shape[0] * self.figureparameters['wellsize'])
        context         = sys.modules['cairo'].Context(CairoImage)

        if not cBack is None:
            context.rectangle(0,0,self.__data.shape[1] * self.figureparameters['wellsize'],self.__data.shape[0] * self.figureparameters['wellsize'])
            context.set_source_rgb(cBack[0],cBack[1],cBack[2])
            context.fill_preserve()
            context.new_path()

        context.set_line_width(self.figureparameters['linewidth'])
        context.translate(.5 * self.figureparameters['wellsize'],.5 * self.figureparameters['wellsize'])
        for x in range(int(self.__data.shape[1])):
            for y in range(int(self.__data.shape[0])):
                context.new_path()
                context.arc(0,0,self.figureparameters['wellradius'],0,2*math.pi)
                r = (self.__data[y,x] - dmin)/(dmax - dmin)
                c = r * cFull + (1 - r) * cEmpty
                context.set_source_rgb(c[0],c[1],c[2])
                context.fill_preserve()
                if r < threshold:
                    context.set_source_rgb(cBorderNoGrowth[0],cBorderNoGrowth[1],cBorderNoGrowth[2])
                else:
                    context.set_source_rgb(cBorder[0],cBorder[1],cBorder[2])
                context.stroke_preserve()
                context.translate(0,self.figureparameters['wellsize'])
            context.translate(self.figureparameters['wellsize'],-self.__data.shape[0] * self.figureparameters['wellsize'])
        
        CairoImage.write_to_png(self.__outfilename)











class GnuplotMSPOutput(object):
    def __init__(self,**kwargs):
        self.__outfilename   = kwargs.get('outfilename','out.gnuplot')
        self.__columns       = kwargs.get('GnuplotColumns',3)
        self.__datasize      = kwargs.get('datasize')
        self.__xrange        = np.array(kwargs.get('GnuplotRange',[1e-3,1e2,0,0])[:2],dtype=np.float)
        self.__yrange        = np.array(kwargs.get('GnuplotRange',[0,0,1e2,1e8])[2:],dtype=np.float)
        
        self.fp              = open(self.__outfilename,'w')
        
    def __del__(self):
        self.fp.close()
    
    def write_init(self):
        self.fp.write("set terminal pngcairo enhanced size 1920,1080\n")
        self.fp.write("set output \"{:s}.png\"\n".format(self.__outfilename))
        self.fp.write("set multiplot\n")
        self.fp.write("set border 15 lw 2 lc rgb \"#2e3436\"\n")
        self.fp.write("set tics front\n")
        self.fp.write("set xra [{:e}:{:e}]\n".format(*self.__xrange))
        self.fp.write("set yra [{:e}:{:e}]\n".format(*self.__yrange))
        self.fp.write("set logscale\n")
        self.fp.write("set format \"10^\{%L\}\"\n")
        self.fp.write("set samples 2001")
        ysize = 1./self.__columns
        if self.__datasize % self.__columns == 0:
            xsize = 1./(self.__datasize//self.__columns)
        else:
            xsize = 1./(self.__datasize//self.__columns + 1.)
        self.fp.write("xsize = {:e}\n".format(xsize))
        self.fp.write("ysize = {:e}\n".format(ysize))
        self.fp.write("xoffset = 0\n")
        self.fp.write("yoffset = 0\n")
        self.fp.write("set size {:e},{:e}\n".format(xsize,ysize))
        self.fp.write("n0(abconc,taulambda,ssmic) = 1 + log(abconc / ssmic) / taulambda\n")
        self.fp.write("set label 1 \"empty\" at graph .5,.9 center front\n")
        self.fp.write("set key at graph .95,.4\n")
        self.fp.write("set samples 1001\n")
        self.fp.write("\n")
        
    def write_plot(self,ID,label,basename,curdata, inoculum = None):
        self.fp.write("set origin xoffset + {:d} * xsize, yoffset + {:d} * ysize\n" . format(ID//self.__columns,ID % self.__columns))
        self.fp.write("set label 1 \"{:s}\"\n".format(label.replace('_','-')))
        self.fp.write("plot \\\n")
        if not inoculum is None:
            for x,y in zip(inoculum[0].flatten(),inoculum[1].flatten()):
                self.fp.write("   '+' u ({:e}):({:e}) w p pt 7 ps 1 lc rgb \"#eeeeec\" title \"\",\\\n".format(x,y))
        if 'SP_sMIC' in curdata.keys(): self.fp.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#4e9a06\" title \"SP B(N)\",\\\n".format(curdata['SPBN_tau'],curdata['SP_sMIC']))
        if 'SP_sMIC' in curdata.keys(): self.fp.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#8ae234\" title \"SP N(B)\",\\\n".format(curdata['SPNB_tau'],curdata['SP_sMIC']))
        if 'BN_sMIC' in curdata.keys(): self.fp.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#a40000\" title \"B(N)\",\\\n".format(curdata['BN_tau'],curdata['BN_sMIC']))
        if 'NB_sMIC' in curdata.keys(): self.fp.write("  n0(x,{:e},{:e}) w l lw 4 lc rgb \"#ef2929\" title \"N(B)\",\\\n".format(curdata['NB_tau'],curdata['NB_sMIC']))
        self.fp.write("  \"{:s}\" u 1:2 w p pt 7 ps 2 lc rgb \"#3465a4\" title \"\"\n".format(basename + '.threshold'))
        self.fp.write("\n")
        
