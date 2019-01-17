#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

import platereaderclass as prc

def main():
    parser = argparse.ArgumentParser()
    parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
    parser_io.add_argument("-i","--infiles", nargs = "*")
    parser_io.add_argument("-o","--outfileprefix", default = "", type = str)
    parser_io.add_argument("-D","--designassignment", default = None, nargs = "*")
    parser_io.add_argument("-P","--WritePNG", default = False, action = "store_true")
    parser_io.add_argument("-E","--ErrorEstimates", default = False, action = "store_true")
    parser_io.add_argument("-v","--verbose", default = False, action = "store_true")
    args = parser.parse_args()

    data = prc.PlateReaderData(**vars(args))

    for dataID in range(int(data)):
        title = data.titles[dataID]
        print title
        if len(args.outfileprefix) > 0: fn = (args.outfileprefix + '_' + title).replace(' ','_')
        else:                           fn = title.replace(' ','_')
        
        data.write_data_to_file(dataID = dataID,filename = fn,include_error_estimates = args.ErrorEstimates)
        if args.WritePNG:
            threshold = data.EstimateGrowthThreshold(dataID,historange=(-7,1),bins = 20)
            data.write_PNG(dataID,growththreshold = threshold)
            
if __name__ == "__main__":
    main()
    
