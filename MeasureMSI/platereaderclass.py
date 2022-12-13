#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
import itertools

from skimage import measure

# getting rid of deprecation warning in sklearn
try:
    from collections import Sequence
except:
    from collections.abc import Sequence


class PlateReaderData(object):
    def __init__(self,**kwargs):

        self.__infilenames              = kwargs.get("infiles",[])
        if self.__infilenames is None:    raise IOError('No input data provided. Use option -i FILENAME1 [FILENAME2 ...]')
        self.__ignoresheetsparameter    = kwargs.get("ignoresheets",[])

        self.__designdata               = list()
        self.__designtitle              = list()

        self.__data                     = list()
        self.__filenames                = list()
        self.__sheetnames               = list()
        self.__designassignment         = list()

        self.__rescale                  = kwargs.get("DataRescale",False)
        self.__logscale                 = kwargs.get("DataLogscale",False)
        self.__logmin                   = kwargs.get("DataLogscaleMin",-20)

        self.__force_orientation        = kwargs.get("ForceOrientation", False)
        self.__force_equaldesignspacing = kwargs.get("ForceEqualDesignSpacing", False)
        self.__threshold_weight         = kwargs.get("ThresholdWeight", 0.5)
        self.__threshold_logscale       = kwargs.get("ThresholdLog", False)

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

        self.__ignoresheets             =  ['Plate Design*']
        if len(self.__ignoresheetsparameter) > 0:
            self.__ignoresheets        += self.__ignoresheetsparameter

        # set threshold to not yet computed
        self.__UseBinarizedData         = kwargs.get('UseBinarizedData',False)
        self.__threshold                = None
        if self.__UseBinarizedData:
            self.__theshold = 0.5

        # load all data at __init__()
        if len(self.__infilenames) > 0:
            i = 0
            for fn in self.__infilenames:
                try:
                    tmpExcelFile = pd.ExcelFile(fn, engine = 'openpyxl')
                except:
                   continue

                for designsheet in [s for s in tmpExcelFile.sheet_names if self.ignoresheet(s)]:
                    tmpDesignData0 = self.read_sheet(tmpExcelFile,designsheet,**self.__read_coordinates_design0)
                    tmpDesignData1 = self.read_sheet(tmpExcelFile,designsheet,**self.__read_coordinates_design1)
                    self.__designdata.append((tmpDesignData0,tmpDesignData1))
                    self.__designtitle.append(designsheet)
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


        # force equal spacing between wells and ignore individual values in design sheets, except for max and min
        if self.__force_equaldesignspacing:
            for i in range(len(self.__designdata)):
                Bstart, Bstop = self.__designdata[i][0][0,0], self.__designdata[i][0][-1,-1]
                Nstart, Nstop = self.__designdata[i][1][0,0], self.__designdata[i][1][-1,-1]
                Ncount, Bcount = self.__designdata[i][0].shape

                Bvalues = np.repeat([np.exp(np.linspace(np.log(Bstart), np.log(Bstop), Bcount))], Ncount, axis = 0)
                Nvalues = np.repeat([np.exp(np.linspace(np.log(Nstart), np.log(Nstop), Ncount))], Bcount, axis = 0).T

                self.__designdata[i] = [Bvalues, Nvalues]


        # force orientation of all data to be the same
        if self.__force_orientation:
            tmp_switch_designs = []
            for i in range(len(self.__data)):
                orientation0 = orientation1 = 1
                if self.__designdata[self.__designassignment[i]][0][0,0] > self.__designdata[self.__designassignment[i]][0][-1,-1]: orientation0 = -1
                if self.__designdata[self.__designassignment[i]][1][0,0] < self.__designdata[self.__designassignment[i]][1][-1,-1]: orientation1 = -1

                if orientation0 != 1 or orientation1 != 1:
                    self.__data[i] = self.__data[i][::orientation1,::orientation0]

                    if self.__designassignment[i] not in tmp_switch_designs:
                        tmp_switch_designs.append(self.__designassignment[i])

            for i in tmp_switch_designs:
                orientation0 = orientation1 = 1
                if self.__designdata[i][0][0,0] > self.__designdata[i][0][-1,-1]: orientation0 = -1
                if self.__designdata[i][1][0,0] < self.__designdata[i][1][-1,-1]: orientation1 = -1

                self.__designdata[i] = (self.__designdata[i][0][::orientation1, ::orientation0],self.__designdata[i][1][::orientation1, ::orientation0])





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
            PlateImage(self[i], outfilename = self.__sheetnames)


    def WriteData(self, dataID, filename = 'out', include_error_estimates = False):
        xlist = self.__designdata[self.__designassignment[dataID]][0]
        ylist = self.__designdata[self.__designassignment[dataID]][1]
        zlist = self[dataID]
        s     = np.shape(zlist)
        fp    = open(filename,'w')

        for i in range(s[0]):
            for j in range(s[1]):
                if include_error_estimates:
                    fp.write('{} {} {} {}\n'.format(xlist[i,j],ylist[i,j],zlist[i,j],elist[i,j]))
                else:
                    fp.write('{} {} {}\n'.format(xlist[i,j],ylist[i,j],zlist[i,j]))
            fp.write('\n')

        fp.close()


    def RelativePosToInoculum(self, position, designID = None, dataID = None):
        if designID is None:
            if dataID is None:
                raise IOError
            else:
                design0 = self.__designdata[self.__designassignment[dataID]][0]
                design1 = self.__designdata[self.__designassignment[dataID]][1]
        else:
            design0 = self.__designdata[designID][0]
            design1 = self.__designdata[designID][1]

        designsize = np.array(np.shape(design0),dtype = np.float) - np.ones(2)

        i = np.array(np.trunc(position * designsize),dtype = np.int)
        f = position * designsize - i

        inoc0 = design0[i[0],i[1]]
        inoc1 = design1[i[0],i[1]]

        if i[0] < designsize[0]:
            inoc0 *= np.power(design0[i[0]+1,i[1]]/design0[i[0],i[1]],f[0])
            inoc1 *= np.power(design1[i[0]+1,i[1]]/design1[i[0],i[1]],f[0])

        if i[1] < designsize[1]:
            inoc0 *= np.power(design0[i[0],i[1]+1]/design0[i[0],i[1]],f[1])
            inoc1 *= np.power(design1[i[0],i[1]+1]/design1[i[0],i[1]],f[1])

        return np.array([inoc0,inoc1], dtype = np.float)


    def WriteArrayWithDesign(self, filename, data, designID):
        fp = open(filename,'w')
        datashape = np.shape(data)
        for i in range(datashape[0]):
            for j in range(datashape[1]):
                fp.write('{:.6e} {:.6e} {:.6e}\n'.format(*self.RelativePosToInoculum(position = [i/(1.*datashape[0]-1),j/(1.*datashape[1]-1)],designID = designID),data[i,j]))
            fp.write('\n')
        fp.close()


    ##########################
    # data analysis routines #
    ##########################

    def EstimateGrowthThreshold(self, dataID = None, logscale = True, weight = None):
        # otsu's method
        # described in IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS (1979)
        # usually used to binarize photos into black/white, here we separate the growth/no-growth transition
        # modified to use 'uniform distribution', ie in each bin of the histogram we have one value, but spacing between bins is different

        if weight is None: weight = self.__threshold_weight

        # prepare data
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

        # need sorted data
        sx = np.sort(x)
        lx = len(x)

        # actual algorithm
        w   = np.arange(lx,dtype=np.float)/lx
        m   = np.array([np.mean(sx[:k]) * w[k] if k > 0 else 0 for k in range(lx)])
        mT  = np.mean(sx)

        sB  = np.array([(mT * w[k] - m[k])**2/(w[k]*(1.-w[k])) if w[k]*(1.-w[k]) > 0 else 0 for k in range(lx)])
        idx = np.argmax(sB)

        try:
            avgthr = (1 - weight)*sx[idx] + weight*sx[idx+1]
        except:
            avgthr = sx[idx]

        if logscale:
            return np.exp(avgthr)
        else:
            return avgthr


    def ComputeThreshold(self, weight = None):
        self.__threshold = self.EstimateGrowthThreshold(dataID = None, weight = weight)


    def get_noise_estimates(self,dataID):
        # get rough estimate of noise as variance between neighboring wells on plate
        shape = np.shape(self[dataID])
        ne = np.zeros(shape)

        #corners
        ne[0,0]                   = np.std(self[dataID][:2,:2])
        ne[0,shape[1]-1]          = np.std(self[dataID][:2:,-2:])
        ne[shape[0]-1,0]          = np.std(self[dataID][-2:,:2])
        ne[shape[0]-1,shape[1]-1] = np.std(self[dataID][-2:,-2:])

        # edges
        for i in range(1,shape[0]-1):
            ne[i,0]               = np.std(self[dataID][i-1:i+1,:2])
            ne[i,shape[1]-1]      = np.std(self[dataID][i-1:i+1,-2:])

        # edges
        for j in range(1,shape[1]-1):
            ne[0,j]               = np.std(self[dataID][:2,j-1:j+1])
            ne[shape[0]-1,j]      = np.std(self[dataID][-2:,j-1:j+1])

            # bulk
            for i in range(1,shape[0]-1):
                ne[i,j]           = np.std(self[dataID][i-1:i+1,j-1:j+1])

        return ne




    def BinarizedData(self, dataID, threshold = None):
        if threshold is None:
            if self.__threshold is None:
                threshold = self.EstimateGrowthThreshold()
            else:
                threshold = self.__threshold
        return np.where(self.__data[dataID] < threshold, 0, 1)



    #########################################################################
    # compute transitions between growth and no-growth with various methods #
    #########################################################################

    def transitions(self,threshold = None, useGPR = False, weight = None):
        if threshold is None:
            if self.__threshold is None:
                threshold = self.EstimateGrowthThreshold(dataID = None, weight = weight)
            else:
                threshold = self.__threshold
        if useGPR:
            return [(self.__filenames[i], self.__sheetnames[i], self.compute_growth_nogrowth_transition_GPR(i,threshold)) for i in range(len(self.__data))]
        else:
            return [(self.__filenames[i], self.__sheetnames[i], self.compute_growth_nogrowth_transition(i,threshold)) for i in range(len(self.__data))]


    def compute_growth_nogrowth_transition(self, dataID, threshold = None, geom = True, onlyLongest = True, weight = None):

        if threshold is None:
            if self.__threshold is None:
                threshold = self.EstimateGrowthThreshold(dataID, weight = weight)
            else:
                threshold = self.__threshold

        if threshold == 0 and self.__threshold_logscale:
            threshold = 0.9 * np.min(self[dataID])

        platesize_m1      = np.array(np.shape(self[dataID]),dtype=np.float) - np.ones(2)
        if not self.__threshold_logscale:
            contours      = measure.find_contours(self[dataID],threshold)
        else:
            contours      = measure.find_contours(np.log(self[dataID]), np.log(threshold))

        threshold_contour = list()
        contour_length    = np.array([len(c) for c in contours],dtype = np.int)
        for i,c in enumerate(contours):
            if i == contour_length.argmax() or not onlyLongest:
                for pos in c:
                    threshold_contour.append(self.RelativePosToInoculum(position = pos/platesize_m1, designID = self.__designassignment[dataID]))

        if len(threshold_contour) > 0:
            return np.vstack(threshold_contour)
        else:
            return None



    ###############################
    # Gaussian Process Regression #
    ###############################

    # first get a much finer grid of interpolated data, not just the 8x12 wells on a plate
    # then estimate curve through this surface on fine grid at threshold

    def GaussianProcessRegression(self,dataID,kernellist = ['white','matern'], restarts_optimizer = 10, outputgrid = (20,20), FitToIndexGrid = False):

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

        # define input grid on which data is defined, define output grid
        if isinstance(outputgrid,(list,tuple,np.ndarray)):
            outshape = (outputgrid[0],outputgrid[1])
        elif isinstance(outputgrid,int):
            outshape = (outputgrid,outputgrid)
        else:
            raise TypeError

        designshape = np.shape(self.__designdata[self.__designassignment[dataID]][0])
        fitshape    = (designshape[0] * designshape[1],1)
        predshape   = (outshape[0] * outshape[1],1)

        if FitToIndexGrid:
            # (1) either simply as grid of indices
            # input grid
            datagrid0      = (np.repeat([np.arange(designshape[0], dtype = np.float)], axis = 0, repeats = designshape[1])).reshape(fitshape)[:,0]
            datagrid1      = (np.repeat([np.arange(designshape[1], dtype = np.float)], axis = 0, repeats = designshape[0]).T).reshape(fitshape)[:,0]
            design         = np.array([[x,y] for x,y in zip(datagrid0,datagrid1)])

            # fitting grid
            fitgrid0       = (np.repeat([np.linspace(0, datagrid0[-1], num = outshape[0])], axis = 0, repeats = outshape[1])).reshape(predshape)[:,0]
            fitgrid1       = (np.repeat([np.linspace(0, datagrid1[-1], num = outshape[1])], axis = 0, repeats = outshape[0]).T).reshape(predshape)[:,0]
            fitgr          = np.array([[x,y] for x,y in zip(fitgrid0,fitgrid1)])


        else:
            # (2) or on the orignal inoculum data:
            # reformat input data into correct array sizes
            # input grid
            datagrid0      = self.__designdata[self.__designassignment[dataID]][0].reshape(fitshape)[:,0]
            datagrid1      = self.__designdata[self.__designassignment[dataID]][1].reshape(fitshape)[:,0]
            design         = np.array([[np.log(x),np.log(y)] for x,y in zip(datagrid0,datagrid1)])

            # fitting grid
            fitgrid0          = (np.repeat([np.linspace(np.log(datagrid0[0]), np.log(datagrid0[-1]), num = outshape[0])], axis = 0, repeats = outshape[1])).reshape(predshape)[:,0]
            fitgrid1          = (np.repeat([np.linspace(np.log(datagrid1[0]), np.log(datagrid1[-1]), num = outshape[1])], axis = 0, repeats = outshape[0]).T).reshape(predshape)[:,0]
            fitgrid           = np.array([[x,y] for x,y in zip(fitgrid0,fitgrid1)])

        # define kernels for Gaussian Process
        kernel = generate_kernel(kernellist)

        # initiate GPR
        gp = sklgp.GaussianProcessRegressor(kernel = kernel, n_restarts_optimizer = restarts_optimizer)

        # set 'training data' to measurements from plates
        platedata = np.reshape(self[dataID],fitshape)

        # estimate hyperparamters for kernels
        gp.fit(design,platedata)

        # use GPR to estimate values on (fine) grid
        platedata_prediction = gp.predict(fitgrid)

        return platedata_prediction.reshape(outshape, order = 'A')


    def compute_growth_nogrowth_transition_GPR(self, dataID, threshold = None, gridsize = 24, kernellist = ['white','matern'], SaveGPRSurfaceToFile = False, FitToIndexGrid = False, onlyLongest = True, weight = None):

        if self.__UseBinarizedData:
            threshold = 0.5
        if threshold is None:
            threshold  = self.EstimateGrowthThreshold(dataID = None, weight = weight)

        outgridsize    = (gridsize,gridsize)
        outgridsize_m1 = np.array(outgridsize) - np.ones(2)
        pdpred         = self.GaussianProcessRegression(dataID, outputgrid = outgridsize, kernellist = kernellist, FitToIndexGrid = FitToIndexGrid)
        contours       = measure.find_contours(pdpred, threshold)


        if SaveGPRSurfaceToFile:
            filename = self.titles[dataID].replace(' ','_') + '.gprsurface'
            self.WriteArrayWithDesign(filename, data = pdpred, designID = self.__designassignment[dataID])

        threshold_contour = list()
        contour_length    = np.array([len(c) for c in contours],dtype = np.int)
        for i,c in enumerate(contours):
            if i == contour_length.argmax() or not onlyLongest:
                for pos in c:
                    threshold_contour.append(self.RelativePosToInoculum(position = pos/outgridsize_m1, designID = self.__designassignment[dataID]))

        threshold_contour = np.vstack(threshold_contour)
        return threshold_contour



    ################
    # python stuff #
    ################

    def __getitem__(self,key):
        if self.__UseBinarizedData:
            return self.BinarizedData(dataID = key)
        else:
            return self.__data[key]


    def __int__(self):
        return len(self.__data)


    def __len__(self):
        return len(self.__data)


    def __iter__(self):
        for i in range(len(self)):
            yield self.__filenames[i], self.__sheetnames[i], self[i], self.__designassignment[i]


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
        elif key == "threshold":
            if self.__theshold is None:
                if self.__UseBinarizedData:
                    self.__theshold = 0.5
                else:
                    self.ComputeThreshold()
            return self.__threshold








