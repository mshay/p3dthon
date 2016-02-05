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

class DumpID(object):

    def __init__(self, param_file=None):
        self.param = load_param(param_file)


    def get_part_in_box(self,r0,dx):
        pass


    def _get_procs_in_box(self, x0, dx, y0, dy, z0=None, dz=None):
        """
        Takes the real r postion and returns what dump file
        partilces coresponding to that position will be on, as well as
        the index position of the list of processeors on that dump file.

        """

        proc_dx = np.array([self.param['lx']/self.param['pex'],
                            self.param['ly']/self.param['pey'],
                            self.param['lz']/self.param['pez']])

        if z0 is None and dz is None:
            z0 = self.param['lz']/2.
            dz = self.param['lz']

        r0 = np.array([x0,y0,z0])
        dx = np.array([dx,dy,dz])

        procs_needed = [] # The corners of the cube 
        
        # find the lower left most proc
        procs_needed.append(self._r0_to_proc(*(r0 - dx/2.))
    
        for c in range(3):
            rng.append(np.arange(r0[c] - dx[c]/2., r0[c] + dx[c]/2., proc_dx[c]))
            if rng[c][-1] < r0[c] + dx[c]/2.:
                rng[c] = np.hstack((rng[c],r0[c] + dx[c]/2.))


        pass


    def _r0_to_proc(self, x0, y0, z0):
        """ Returns the px,py,pz processeor for a given values of x, y, and z
        """

        lx = self.param['lx']
        ly = self.param['ly']
        lz = self.param['lz']
        
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

        px = int(np.floor(x0/self.param['lx']*self.param['pex'])) + 1
        py = int(np.floor(y0/self.param['ly']*self.param['pey'])) + 1
        pz = int(np.floor(z0/self.param['lz']*self.param['pez'])) + 1

        return px,py,pz


    def _proc_to_dumplocation(self, px=1, py=1, pz=1):
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

        pex = self.param['pex']
        pey = self.param['pey']
        pez = self.param['pez']
        nch = self.param['nchannels']

        if pex*pey*pez%nch != 0:
            raise NotImplementedError()

        dump_IO_version = 'V1'

        if self.param.has_key('USE_IO_V2'):
            dump_IO_version = 'V2'
       
        if dump_IO_version == 'V1':
            N = (px - 1)%nch + 1
            R = (pz - 1)*(pex/nch)*(pey) + (pex/nch)*(py - 1) + (px - 1)/nch

        else: # dump_IO_version == 'V2'

            npes_per_dump = pex*pey*pez/nch

            pe = (pz - 1)*pex*pey + (py - 1)*pex + (px - 1)

            N = pe/npes_per_dump + 1
            R = pe%npes_per_dump

        return num_to_ext(N),R



