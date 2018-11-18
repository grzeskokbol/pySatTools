#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""Usage:
    CoordConvert.py <date> <X> <Y> <Z>  <VX> <VY> <VZ>  <systemIn2systemOut> [-v] [-h]

Arguments:
    <date> Date in the ISO format: 2017-01-01T12:00:00
    <X> <Y> <Z> : object's coordinates [m]
    <VX> <VY> <VZ> : object's velocity [m/s]

Optional arguments:

    -v, --verbose           Verbose mode
    -h, --help              Shows help message and exit
"""

import dateutil.parser
import argparse
import logging

import numpy as np


import libpySat as pySat

if __name__ == '__main__':
    # Just some parsing to datetime formats, then calling the actual computation file.
    parser = argparse.ArgumentParser(description='Tranformation between coordinate systems.')
    parser.add_argument('date', type=str, help='Date in ISO format: 2017-01-01T12:00:00')
    parser.add_argument('X', type=float,      help='X coordinate')
    parser.add_argument('Y', type=float, help='Y coordinate')
    parser.add_argument('Z', type=float, help='Z coordinate')
    parser.add_argument('VX', type=float,      help='X velocity')
    parser.add_argument('VY', type=float, help='Y velocity')
    parser.add_argument('VZ', type=float, help='Z velocity')
    parser.add_argument('systemIn2systemOut', type=str, help="Define transformation")
    parser.add_argument("-v", "--verbose", required=False,action="store_true", help="Increase output verbosity" )
    args = parser.parse_args()
    epoch = dateutil.parser.parse(args.date)
    systemInOut=args.systemIn2systemOut

    # Get the state [X,V]
    posXYZ_in=np.arange((3),dtype=float)
    velXYZ_in = np.arange((3), dtype=float)

    posXYZ_in[0] = args.X
    posXYZ_in[1] = args.Y
    posXYZ_in[2] = args.Z

    velXYZ_in[0] = args.VX
    velXYZ_in[1] = args.VY
    velXYZ_in[2] = args.VZ

    stateXYZ_in = np.arange((6), dtype=float)
    stateXYZ_in[0:3]=posXYZ_in
    stateXYZ_in[3:6] =velXYZ_in

    if (args.verbose):
        logging.basicConfig(level=logging.INFO)
    logging.info('# Epoch: %s ' % epoch )
    logging.info('# Input State [X Y Z]: %12.4f %12.4f %12.4f' %
                 (posXYZ_in[0],posXYZ_in[1],posXYZ_in[2]))
    logging.info('# Input State [Vx Vy Vz]: %8.4f %8.4f %8.4f' %
                 (velXYZ_in[0],velXYZ_in[1],velXYZ_in[2]))
    logging.info('#Transformation: %s' % systemInOut)

    #pySatTr.pySatTransform() # getUT1(epoch) # = pySatTransform.pySatTransform()
    tr= pySat.Transform()
    if systemInOut == "ITRF2ICRF" :
        state_out = tr.TRF2CRF(epoch, stateXYZ_in)
    if systemInOut == "ICRF2ITRF":
        state_out = tr.CRF2TRF(epoch, stateXYZ_in)
    print('#Output coordinates at %s : %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f ' % (epoch, state_out[0],state_out[1],state_out[2],state_out[3],state_out[4],state_out[5]))

    exit(0)
