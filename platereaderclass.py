#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import itertools
import openpyxl
import cairo

class PlateReaderData(object):
    def __init__(self,**kwargs):

        self.__infilenames = kwargs.get("infiles",[])
        self.__ignoresheetsparameter = kwargs.get("ignoresheets",[])

        self.extract_figure_file_parameters(kwargs)

        self.__designdata = list()
        self.__designtitle = list()

        self.__data = list()
        self.__filenames = list()
        self.__sheetnames = list()
        self.__designassignment = list()
        
        self.__rescale = kwargs.get("DataRescale",False)
        self.__logscale = kwargs.get("DataLogscale",False)
        self.__logmin  = kwargs.get("DataLogscaleMin",-20)
        
        self.__read_coordinates = { 'x0':kwargs.get("xstart",2),
                                    'y0':kwargs.get("ystart",2),
                                    'xwidth':kwargs.get("xwidth",12),
                                    'yheight':kwargs.get("yheight",8) }
        
        self.__ignoresheets = ['Plate Design*']
        if len(self.__ignoresheetsparameter) > 0:
            self.__ignoresheets += self.__ignoresheetsparameter
        
        # load all data at __init__()
        for fn in self.__infilenames:
            try:
                data = openpyxl.load_workbook(fn)
            except:
                continue
            
            for designsheet in [s for s in data if self.ignoresheet(s.title)]:
                self.__designdata.append(self.read_initial_conditions(data,designsheet.title))
                self.__designtitle.append(designsheet.title)
            i = 0
            for sheet in [s for s in data if not self.ignoresheet(s.title)]:
                self.__data.append(self.read_sheet(sheet))
                self.__sheetnames.append(sheet.title)
                self.__filenames.append(fn)
                try:
                    design = kwargs.get("designassignment",[])[i]
                    if not 0 <= design < len(self.__designdata):
                        raise KeyError
                    self.__designassignment.append(design)
                except:
                    self.__designassignment.append(0)
                i+=1
            
        
    
    
    def column_string(self,n):
        div=n
        string=""
        temp=0
        while div>0:
            module=(div-1)%26
            string=chr(65+module)+string
            div=int((div-module)/26)
        return string
    
    
    def rgb(self,color):
        r = int(color[0:2],16)
        g = int(color[2:4],16)
        b = int(color[4:6],16)
        return np.array([r/255.,g/255.,b/255.])


    def avg(self,values,geom = True):
        if len(values) > 0:
            if geom:    return np.power(np.product(values),1./len(values))
            else:       return np.mean(values)
        else:
            raise ValueError
    
    
    def read_sheet(self,sheet,x0 = None, y0 = None, width = None, height = None):
        if x0 is None:      x0 = self.__read_coordinates['x0']
        if y0 is None:      y0 = self.__read_coordinates['y0']
        if width is None:   width = self.__read_coordinates['xwidth']
        if height is None:  height = self.__read_coordinates['yheight']
        
        return np.array([[sheet['{:s}{:d}'.format(self.column_string(i),j)].value for i in range(x0,x0+width)] for j in range(y0,y0+height)])


    def read_initial_conditions(self,data,sheetname,xab = 4, yab = 14, xcells = 4, ycells = 3, width = 12, height = 8):
        if sheetname in data.sheetnames:
            return self.read_sheet(data[sheetname],x0 = xab,y0 = yab),self.read_sheet(data[sheetname],x0 = xcells,y0 = ycells)
        else:
            raise KeyError
    

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
    
    
    def rescale(self,g):
        
        if self.__logscale:
            r = np.log(g)
            r[r<self.__logmin] = self.__logmin
        else:
            r = g[:,:]
            
        r = (r - np.min(r))/(np.max(r) - np.min(r))
        return r
    
    
    def extract_figure_file_parameters(self,kwargs):
        self.figureparameters = {   'colors':   ['3465a4','ffffff','2e3436','eeeeec'],
                                    'wellradius': 20,
                                    'wellsize':50,
                                    'linewidth':3}




    def write_PNG(self,dataid,outfilename = None):
        if 0 <= dataid < len(self.__data):
            if outfilename is None:
                outfilename = self.__sheetnames[dataid] + '.png'
                
            cFull   = self.rgb(self.figureparameters['colors'][0])
            cEmpty  = self.rgb(self.figureparameters['colors'][1])
            cBorder = self.rgb(self.figureparameters['colors'][2])
            cBack   = self.rgb(self.figureparameters['colors'][3])

            CairoImage = cairo.ImageSurface(cairo.FORMAT_ARGB32,self.__data[dataid].shape[1] * self.figureparameters['wellsize'],self.__data[dataid].shape[0] * self.figureparameters['wellsize'])
            context    = cairo.Context(CairoImage)

            context.rectangle(0,0,self.__data[dataid].shape[1] * self.figureparameters['wellsize'],self.__data[dataid].shape[0] * self.figureparameters['wellsize'])
            context.set_source_rgb(cBack[0],cBack[1],cBack[2])
            context.fill_preserve()
            context.new_path()

            context.set_line_width(self.figureparameters['linewidth'])
            context.translate(.5 * self.figureparameters['wellsize'],.5 * self.figureparameters['wellsize'])
            datamax   = np.amax(self.__data[dataid])
            datarange = np.amax(self.__data[dataid]) - np.amin(self.__data[dataid])
            for x in range(int(self.__data[dataid].shape[1])):
                for y in range(int(self.__data[dataid].shape[0])):
                    context.new_path()
                    context.arc(0,0,self.figureparameters['wellradius'],0,2*math.pi)
                    r = (datamax - self.__data[dataid][y,x])/datarange
                    c = cFull * (1 - r) + cEmpty * r
                    context.set_source_rgb(c[0],c[1],c[2])
                    context.fill_preserve()
                    context.set_source_rgb(cBorder[0],cBorder[1],cBorder[2])
                    context.stroke_preserve()
                    context.translate(0,self.figureparameters['wellsize'])
                context.translate(self.figureparameters['wellsize'],-self.__data[dataid].shape[0] * self.figureparameters['wellsize'])
            
            CairoImage.write_to_png(outfilename)

    def Export_All_PNGs(self):
        for i in range(self.count):
            self.write_PNG(i)

    def get_design(self,designid = 0, designname = None):
        if not designname is None:
            if designname in self.__designtitle:
                return self.__designdata[self.__designtitle.index(designname)]
        else:
            return self.__designdata[designid]
    
    
    def get_designassignment(self,dataid = None):
        if dataid is None:
            return self.__designassignment
        else:
            if 0 <= dataid < len(self.__designassignment):
                return self.__designassignment[dataid]
            else:
                raise KeyError
    
    
    def count_design(self):
        assert len(self.__designdata) == len(self.__designtitle)
        return len(self.__designdata)


    def all_values(self):
        return np.concatenate([x.flatten() for x in self.__data])
    

    def compute_growth_nogrowth_transition(self,dataid,threshold,design = 0,geom = True):
        r=list()
        for i,j in itertools.product(np.arange(np.shape(self.__data[dataid])[0]),np.arange(np.shape(self.__data[dataid])[1])):
            try:
                if (self.__data[dataid][i,j] - threshold) * (self.__data[dataid][i+1,j] - threshold) < 0:
                    r.append(np.array([self.avg([self.__designdata[design][0][i,j],self.__designdata[design][0][i+1,j]],geom),self.avg([self.__designdata[design][1][i,j],self.__designdata[design][1][i+1,j]],geom)]))
            except:
                pass
            try:
                if (self.__data[dataid][i,j] - threshold) * (self.__data[dataid][i,j+1] - threshold) < 0:
                    r.append(np.array([self.avg([self.__designdata[design][0][i,j],self.__designdata[design][0][i,j+1]],geom),self.avg([self.__designdata[design][1][i,j],self.__designdata[design][1][i,j+1]],geom)]))
            except:
                pass
        if len(r) > 0:
            return np.vstack(r)
        else:
            return None


    def transitions(self,threshold):
        return [(self.__filenames[i],self.__sheetnames[i],self.compute_growth_nogrowth_transition(i,threshold,self.__designassignment[i])) for i in range(len(self.__data))]


    def __iter__(self):
        for fn,title,platedata in zip(self.__filenames,self.__sheetnames,self.__data):
            yield fn,title,self.rescale(platedata)
    
    
    def __len__(self):
        return len(self.__data)
    
    def __getattr__(self,key):
        if key == "count":
            return len(self.__data)
        elif key == "count_design":
            return len(self.__designdata)
        else:
            super(PlateReaderData,self).__getitem__(key)
    
