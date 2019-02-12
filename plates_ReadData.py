#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

import platereaderclass as prc

def main():
    parser = argparse.ArgumentParser()
    parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
    parser_io.add_argument("-i","--infiles",          nargs   = "*")
    parser_io.add_argument("-D","--designassignment", default = None, nargs = "*")
    parser_io.add_argument("-o","--outfileprefix",    default = "", type = str)
    parser_io.add_argument("-P","--WritePNG",         default = False, action = "store_true")
    parser_io.add_argument("-S","--WriteSVG",         default = False, action = "store_true")
    parser_io.add_argument("-E","--ErrorEstimates",   default = False, action = "store_true")
    parser_io.add_argument("-v","--verbose",          default = False, action = "store_true")
    
    
    parser_Figure = parser.add_argument_group(description = "==== parameters for graphical output ====")
    parser_Figure.add_argument("-x","--FigureWellDistance",        default = 50, type = int)
    parser_Figure.add_argument("-r","--FigureWellRadius",          default = 18, type = int)
    parser_Figure.add_argument("-L","--FigureLinewidth",           default =  3, type = int)
    parser_Figure.add_argument("-c","--FigureColorEmpty",          default = "d3d7cf")
    parser_Figure.add_argument("-C","--FigureColorFull",           default = "3465a4")
    parser_Figure.add_argument("-B","--FigureColorBackground",     default = None)
    parser_Figure.add_argument("-b","--FigureColorBorder",         default = "2e3436")
    parser_Figure.add_argument("-N","--FigureColorBorderNoGrowth", default = "a40000")
    parser_Figure.add_argument("-T","--FigureEstimateThreshold",   default = False, action = "store_true")
    
    args = parser.parse_args()

    data = prc.PlateReaderData(**vars(args))

    for dataID in range(int(data)):
        title = data.titles[dataID]
        if args.verbose: print title
        if len(args.outfileprefix) > 0: fn = (args.outfileprefix + '_' + title).replace(' ','_')
        else:                           fn = title.replace(' ','_')
        fn += '.data'
        
        data.write_data_to_file(dataID = dataID,filename = fn,include_error_estimates = args.ErrorEstimates)

        threshold = None
        if args.FigureEstimateThreshold:
            threshold = data.EstimateGrowthThreshold(dataID,historange=(-7,1),bins = 20)
        
        if args.WritePNG:
            prc.PlateImage(data[dataID], outfilename = data.titles[dataID], growththreshold = threshold, outputformat = 'png', **vars(args))
            
        if args.WriteSVG:
            prc.PlateImage(data[dataID], outfilename = data.titles[dataID], growththreshold = threshold, outputformat = 'svg', **vars(args))
            
            
if __name__ == "__main__":
    main()
    
