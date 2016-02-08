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
from dump import Dump
from _methods import load_param
from _methods import _num_to_ext

class DumpID(object):

    def __init__(self, 
                 num=None,
                 param_file=None,
                 path='./'):

        self.dump = Dump(num, param_file, path)
        self.param = self.dump.param


    def get_part_in_box(self,
                        r=[1.,1.],
                        dx=[.5,.5],
                        species=None):

        r0  = [1., 1., .5]
        dx0 = [.5, .5, 1.]

        for c,(r_i,dx_i) in enumerate(zip(r,dx)):
            r0[c]  = r_i
            dx0[c] = dx_i

        dump_and_index = self._get_procs_in_box(r0[0],dx0[0],
                                                r0[1],dx0[1],
                                                r0[2],dx0[2])

        if species is None:
            parts = {'i':[], 'e':[]}
        else:
            parts = {species:[]}

        for d in dump_and_index:
            data = self.dump.read_dump_file(d)
            try:
                data = data[1]
            except KeyError:
                pass
            for sp in parts:
                parts[sp] += [data[sp][g] for g in dump_and_index[d]]
            pdb.set_trace()
        
        for sp in parts:
            for c,p in enumerate(parts[sp]):
                parts[sp][c] = self._trim_parts(p, r0, dx0)
                

        return parts


    def _trim_parts(p, r0, dx0):
        pass


    def _get_procs_in_box(self, x0, dx, y0, dy, z0, dz):
        """
        Takes the real r postion and returns what dump file
        partilces coresponding to that position will be on, as well as
        the index position of the list of processeors on that dump file.

        """

        proc_dx = np.array([self.param['lx']/self.param['pex'],
                            self.param['ly']/self.param['pey'],
                            self.param['lz']/self.param['pez']])

        r0 = np.array([x0,y0,z0])
        dx = np.array([dx,dy,dz])

        procs_needed = [] # The corners of the cube 
        
        # find the lower left most proc
        procs_needed.append(self._r0_to_proc(*(r0 - dx/2.)))
         
        r0_rng = []
        for c in range(3):
            r0_rng.append(np.arange(r0[c] - dx[c]/2., 
                                    r0[c] + dx[c]/2.,
                                    proc_dx[c]))

            if r0_rng[c][-1] < r0[c] + dx[c]/2.:
                r0_rng[c] = np.hstack((r0_rng[c],r0[c] + dx[c]/2.))

        print r0_rng

        p0_rng = []
        for x in r0_rng[0]:
            for y in r0_rng[1]:
                for z in r0_rng[2]:
                    p0_rng.append(self._r0_to_proc(x,y,z))

        p0_rng = set(p0_rng) #This removes duplicates

        di_dict = {}
        for p in p0_rng:
            d = self._proc_to_dumplocation(*p)
            if d[0] in di_dict:
                di_dict[d[0]].append(d[1])
            else:
                di_dict[d[0]] = [d[1]]

        for k in di_dict:
            di_dict[k].sort()
            di_dict[k] = list(set(di_dict[k]))
            #print k, di_dict[k]

        return di_dict


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
            px = 1
        elif x0 >= lx:
            print err_msg.format('X',x0,lx,lx)
            px = self.param['pex']
        else:
            px = int(np.floor(x0/self.param['lx']*self.param['pex'])) + 1

        if y0 < 0.:
            print err_msg.format('Y',y0,ly,0.)
            py = 1
        elif y0 >= ly:
            print err_msg.format('Y',y0,ly,ly)
            py = self.param['pey']
        else:
            py = int(np.floor(y0/self.param['ly']*self.param['pey'])) + 1

        if z0 < 0.:
            print err_msg.format('Z',z0,lz,0.)
            pz = 1
        elif z0 >= lz:
            print err_msg.format('Z',z0,lz,lz)
            pz = self.param['pez']
        else:
            pz = int(np.floor(z0/self.param['lz']*self.param['pez'])) + 1

        return px,py,pz


    def _proc_to_dumplocation(self, px, py, pz):
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

        dump_IO_version = 'V2'

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

        return _num_to_ext(N),R



