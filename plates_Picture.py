#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import cairo

def rgb(color):
    r = int(color[0:2],16)
    g = int(color[2:4],16)
    b = int(color[4:6],16)
    return r/255.,g/255.,b/255.


parser = argparse.ArgumentParser()

parserIO = parser.add_argument_group(description = "=== I/O parameters ===")
parserIO.add_argument("-i","--infile")
parserIO.add_argument("-o","--outfile",default="out.png")

parserFigure = parser.add_argument_group(description = "Option to change output picture")
parserFigure.add_argument("-x","--wellx",default = 50,type=int)
parserFigure.add_argument("-y","--welly",default = 50,type=int)
parserFigure.add_argument("-r","--wellradius",default=18,type=int)
parserFigure.add_argument("-L","--linewidth",default=3,type=int)
parserFigure.add_argument("-c","--ColorEmpty",default="d3d7cf")
parserFigure.add_argument("-C","--ColorFull",default="cc0000")
parserFigure.add_argument("-B","--ColorBackground",default=None)
parserFigure.add_argument("-b","--ColorBorder",default="2e3436")
args = parser.parse_args()


try:
    data = np.genfromtxt(args.infile)
except:
    raise IOError

cE = rgb(args.ColorEmpty)
cF = rgb(args.ColorFull)
cB = rgb(args.ColorBorder)

CairoImage = cairo.ImageSurface(cairo.FORMAT_ARGB32,data.shape[1] * args.wellx,data.shape[0] * args.welly)
context    = cairo.Context(CairoImage)

if not args.ColorBackground is None:
    cBG = rgb(args.ColorBackground)
    context.rectangle(0,0,data.shape[1] * args.wellx,data.shape[0] * args.welly)
    context.set_source_rgb(cBG[0],cBG[1],cBG[2])
    context.fill_preserve()
    context.new_path()

context.set_line_width(args.linewidth)
context.translate(.5 * args.wellx,.5 * args.welly)
for x in range(data.shape[1]):
    for y in range(data.shape[0]):
        #print x,y
        context.new_path()
        context.arc(0,0,args.wellradius,0,2*math.pi)
        if data[y,x] <= 0:
            context.set_source_rgb(cE[0],cE[1],cE[2])
        else:
            context.set_source_rgb(cF[0],cF[1],cF[2])
        context.fill_preserve()
        context.set_source_rgb(cB[0],cB[1],cB[2])
        context.stroke_preserve()
        context.translate(0,args.wellx)
    context.translate(args.welly,-data.shape[0] * args.wellx)

CairoImage.write_to_png(args.outfile)

