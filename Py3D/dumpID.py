#######################################################################
#                                                                     #
#                  Progs        :  dumpID.py                          #
#                  Aruthor      :  Colby Haggerty                     #
#                  Date         :  2016.01.31                         #
#                                                                     #
#                                                                     #
#######################################################################


import os
import sys
import pdb
import numpy as np
from _methods import load_param
from _methods import _num_to_ext

def location_to_dumpindex(x0=.5,
                          y0=.5,
                          z0=.5,
                          param_file=None)
    """
    Takes the real x,y,z postion and returns what dump file
    partilces coresponding to that position will be on, as well as
    the index position of the list of processeors on that dump file.
    """
    pass

def _location_to_proc(param, x0=0., y0=0., z0=0.):
    """ Returns the px,py,pz processeor for a given values of x, y, and z
    """

    lx = param['lx']
    ly = param['ly']
    lz = param['lz']
    
    err_msg = '{0} value {1} is outside of the simulation boundry [0.,{2}].'+\
              'Setting {0} = {3}'

    if x0 < 0.:
        print err_msg.format('X',x0,lx,0.)
        x0 = 0.
    if x0 > lx:
        print err_msg.format('X',x0,lx,lx)
        x0 = lx 

    if y0 < 0.:
        print err_msg.format('Y',y0,ly,0.)
        y0 = 0.
    if y0 > ly:
        print err_msg.format('Y',y0,ly,ly)
        y0 = ly

    if z0 < 0.:
        print err_msg.format('Z',z0,lz,0.)
        z0 = 0.
    if z0 > lz:
        print err_msg.format('Z',z0,lz,lz)
        z0 = lz

    px = int(np.floor(x0/param['lx']*param['pex'])) + 1
    py = int(np.floor(y0/param['ly']*param['pey'])) + 1
    pz = int(np.floor(z0/param['lz']*param['pez'])) + 1

    return px,py,pz


def _proc_to_dumpindex(param, px=1, py=1, pz=1):
    """ Returns the dump index (di), as well as the postion in the array 
        returned in _get_particles(dump_index=di)

    Big Note: There are two ways marc stores procs on dump files:
              an old way and a new way. We need a to distingush
              which way we are using.
    Old Way:
        Scan over Y, Scan over X then Scan over Z
    New Way:
        Scan over X, Scan over Y then Scan over Z
    """

    pex = param['pex']
    pey = param['pey']
    pez = param['pez']
    nch = param['nchannels']

    if pex*pey*pez%nch == 0:
        npes_per_dump = pex*pey*pez/nch
    else:
        raise NotImplementedError()
    
    # Code for new way
    pe = (pz - 1)*pex*pey + (py - 1)*pex + (px - 1)

    N = pe/npes_per_dump + 1
    R = pe%npes_per_dump

    # Code for the old way
    #if pz != 1:
    #    raise NotImplementedError()

    # I think Every Proc has all Y and Z, and only varys on x 
    N = (px - 1)%nch + 1
    R = (pz - 1)*(pex/nch)*(pey) + (pex/nch)*(py - 1) + (px - 1)/nch

    return num_to_ext(N),R
