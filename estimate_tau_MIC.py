#!/usr/bin/env python

import numpy as np
import argparse
import openpyxl
import cairo
import time
import dateutil
import sys,math



def column_string(n):
    div=n
    string=""
    temp=0
    while div>0:
        module=(div-1)%26
        string=chr(65+module)+string
        div=int((div-module)/26)
    return string


def read_sheet(sheet,x0 = 2, y0 = 2, width = 12, height = 8):
    return np.array([[sheet['{:s}{:d}'.format(column_string(i),j)].value for i in range(x0,x0+width)] for j in range(y0,y0+height)],dtype=np.float)


def rescale(g,logscale=False,logmin=-20):
    
    if logscale:
        r = np.log(g)
        r[r<logmin] = 0
    else:
        r = g[:,:]
        
    r = (r - np.min(r))/(np.max(r) - np.min(r))
    return r


def rgb(color):
    r = int(color[0:2],16)
    g = int(color[2:4],16)
    b = int(color[4:6],16)
    return np.array([r/255.,g/255.,b/255.])


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
    




def main():
    ignoresheets = ["Plate design"]

    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*")
    args = parser.parse_args()
    
    for fn in args.infiles:
        try:
            data = openpyxl.load_workbook(fn)
        except:
            raise IOError("could not open file")
    
        print fn
        for sheet in [s for s in data if not s.title in ignoresheets]:
            print sheet.title
            
            outfile = sheet.title.replace(' ','_') + '.png'
            
            growth  = read_sheet(sheet)
            growth  = rescale(growth)
            
            write_PNG(growth,outfile)
            


if __name__ == "__main__":
    main()

