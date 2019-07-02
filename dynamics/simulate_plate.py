#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math

sys.path.append(sys.path[0] + '/..')

import popdyn
import platereaderclass as prc


def computeThreshold(x):
    x   = list(x.flatten())
    sx  = np.sort(x)
    lx  = len(x)
    w   = np.arange(lx,dtype=np.float)/lx
    m   = np.array([np.mean(sx[:k]) * w[k] if k > 0 else 0 for k in range(lx)])
    mT  = np.mean(sx)
    sB  = np.array([(mT * w[k] - m[k])**2/(w[k]*(1.-w[k])) if w[k]*(1.-w[k]) > 0 else 0 for k in range(lx)])
    idx = np.argmax(sB)
    return sx[idx]



def main():

    if not '-B' in sys.argv:
        sys.argv.append('-B')
        sys.argv.append('0')
    B0index = sys.argv.index('-B') + 1
    
    if not '-N' in sys.argv:
        sys.argv.append('-N')
        sys.argv.append('1')
    N0index = sys.argv.index('-N') + 1
    
    if not '-S' in sys.argv:
        sys.argv.append('-S')
        sys.argv.append('1e7')

    outfilename = 'out'
    if '-o' in sys.argv:
        outfilename = sys.argv[sys.argv.index('-o') + 1]

    bsize = 12
    nsize =  8
    
    plate = np.zeros((bsize,nsize))
    
    for ndilution in np.arange(nsize):
        N0 = 2.4e1 * 4**ndilution
        sys.argv[N0index] = str(N0)

        for bdilution in np.arange(bsize):
            B0 = 2e-3 * 2**bdilution
            sys.argv[B0index] = str(B0)

            run = popdyn.main()
            popsize = np.array(run['N0'])
            print('{:.2e} {:.2e} {:.2e}'.format(N0,B0,popsize[-1]))
            
            plate[bdilution,nsize - ndilution - 1] = popsize[-1]

    threshold = computeThreshold(plate)
    prc.PlateImage(data = plate.T, growththreshold = threshold, outfilename = outfilename)



if __name__ == "__main__":
    main()
