
# this module codes the test particle run base

import os
import numpy as np
import scipy.ndimage as ndimage

# now for linking C to python
from numpy.ctypeslib import ndpointer
import ctypes

# for plotting routines
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class TPRun:
    """
    Test Particle run
    """

    #==========================================================
    #==========================================================
    def __init__(self,
                 npart=1,
                 charge=-1.,
                 mass=0.04,
                 run=000,
                 tstart=0,
                 tend=5.,
                 t0=.0,
                 r0=[130.,30.0],
                 dr0=[0.1,0.1],
                 v0=[0.,0.,0.],
                 dv0=[0.1,0.1],
                 dt=.001,
                 loading     = '',     #randu,randn,copy
                 fieldinterp = False):

        """ constructor of the testparticle run object

        @param npart       : # of particles
        @param charge      : charge of the particles
        @param mass        : mass of the particles
        @param run         : run from which the fiels will be read
        @param tstart      : starting time of the integration
        @param tend        : end time of the integration
        @param t0          : time at which positions/velocities are given
        @param r0          : initial position
        @param dr0         : initial spatial deviation
        @param v0          : initial velocity
        @param dv0         : initial velocity deviation
        @param dt          : time step
        @param loading     : either 'randu', 'randn', or 'copy'
        @param fieldinterp : delfault(False) inteprolate fields in time or not

        if loading == 'randu' particles will be loaded randomly in a rectangle
        of size dr0
        if 'randn' is chosen, they will be loaded in a gaussian of spatial std=dr0
        if 'copy' is chosen, r0,v0 have to eb arrays of size (3,npart), dr0 and dv0
        are then disregarded


        @return: a TPRun object

        Creation : 2013-05-01 11:17:42.987369

        """
        self._npart     = npart
        self._charge    = charge
        self._mass      = mass
        self._r0        = None
        self._v0        = None
        self._r0        = np.array([[r0[0]],[r0[1]],[0.]])
        self._v0        = np.array([[v0[0]],[v0[1]],[v0[2]]])

        # keep the user values for plotting etc.

        #TODO attention si loading==user
        # r0,v0 sont des tableaux et dr0,dv0 
        # ne doivent pas etre utilises !
        #- 2013-05-20 08:06:57.972466
        self._r0u       = r0
        self._v0u       = v0
        self._dr0u      = dr0
        self._dv0u      = dv0

        self.dt         = dt
        self.tstart     = tstart
        self.tend       = tend

        self.t0         = t0
        self._nt        = int((tend-tstart)/dt) + 1


        # checks the time interval is well defined 
        #----------------------------------------------------------------------
        if t0 < self.tstart or t0 > self.tend:
            print 'time (%5.3f) should be between tstart(%5.3f)\
                   and tend(%5.3f) (both included)' \
                   % (t0 ,self.tstart,self.tend)

            return None
        #----------------------------------------------------------------------
        # that is the selection time index.
        self._it0 = int((self.t0   - self.tstart)/self.dt) + 1


        #c# # loading method
        #c# self._loading = loading
        #c# if loading.lower() == 'randu':
        #c#     self.load_randu(r0,dr0,v0,dv0)

        #c# elif loading.lower() == 'randn':
        #c#     self.load_randn(r0,dr0,v0,dv0)

        #c# elif loading.lower() == 'user':
        #c#     self._r0 = np.zeros((3,self._npart))
        #c#     self._r0[0,:] = r0[0,:]
        #c#     self._r0[1,:] = r0[1,:]
        #c#     self._v0 = np.zeros((3,self._npart))
        #c#     self._v0 = v0

        #c# else:
        #c#     print 'Ptest : warning, no loading method specified'


        # position and velocity arrays

        self.r   = np.zeros((3, self._npart, self._nt),
                                    order='FORTRAN',
                                    dtype=np.float64)

        self.v   = np.zeros((3, self._npart, self._nt),
                                    order='FORTRAN',
                                    dtype=np.float64)


        # electric and magnetic field seen by each particle
        self.Ep   = np.zeros((3, self._npart, self._nt),
                             order = 'FORTRAN',
                             dtype=np.float64)

        self.Bp   = np.zeros((3, self._npart, self._nt),
                             order = 'FORTRAN',
                             dtype = np.float64)


        # do we interpolate fields in time ?
        self._fldinterp = fieldinterp

        # the PIC run from which he get the fields
        self._run       = run


        # if we interpolate fields in time
        # we actually need to read all the files between the time
        # interval

#cc# Right now we are just assuming that we have an idl file correctly loaded in
#cc# This would be a good point to add in stuf that would use the p3d run object
        if self._fldinterp == False:
            #c# self._E = np.array([CR['exav'],CR['eyav'],CR['ezav']])
            #c# self._B = np.array([CR['bxav'],CR['byav'],CR['bzav']])
            self._E = np.zeros((3,CR['bzav'].shape[1],CR['bzav'].shape[0]),"float32",order='FORTRAN')
            self._B = np.zeros((3,CR['bzav'].shape[1],CR['bzav'].shape[0]),"float32",order='FORTRAN')
            self._E[0,:,:] = CR['exav'].transpose()
            self._E[1,:,:] = CR['eyav'].transpose()
            self._E[2,:,:] = CR['ezav'].transpose()
            self._B[0,:,:] = CR['bxav'].transpose()
            self._B[1,:,:] = CR['byav'].transpose()
            self._B[2,:,:] = CR['bzav'].transpose()
            #c# self._E = np.array([transpose(CR['exav']),transpose(CR['eyav']),transpose(CR['ezav'])])
            #c# self._B = np.array([transpose(CR['bxav']),transpose(CR['byav']),transpose(CR['bzav'])])

            # smooth fields may help... B is fine, E is noisy !!
            for c in range(3):
                self._E[c,:,:] = ndimage.gaussian_filter(self._E[c,:,:],
                                                         sigma=6,
                                                         order=0)

            # checks which component is the out of plane
            #c# if self._run.outofplane == 1:
            #c#     Ey = self._E[1,:,:].copy()
            #c#     Ez = self._E[2,:,:].copy()
            #c#     self._E[1,:,:] = np.copy(Ez)
            #c#     self._E[2,:,:] = np.copy(-Ey)

            #c#     By = self._B[1,:,:].copy()
            #c#     Bz = self._B[2,:,:].copy()
            #c#     self._B[1,:,:] = np.copy(Bz)
            #c#     self._B[2,:,:] = np.copy(-By)



        # in case we want to debug, use analytic fields
        debug = False

        if debug == True:
            self._B[0,:,:]  = 0.#1. + np.random.randn(self._B.shape[1],self._B.shape[2])*1e-3
            self._B[1,:,:]  = 0.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-3
            self._B[2,:,:]  = 1.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-3
            self._E[0,:,:]  = 0.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-2
            self._E[1,:,:]  = 1.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-2
            self._E[2,:,:]  = 0.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-2


        # linking to the C library
        pathlibdir   = os.path.dirname(os.path.realpath(__file__))
        pathlib      = os.path.join(pathlibdir,'functions.so')
        self._lib = ctypes.cdll.LoadLibrary(pathlib)
    #==========================================================

    #==========================================================
    def move(self):
        """moves all the particles

        @return: @todo

        Exemple  : 

        Creation : 2013-05-01 14:28:08.724842

        """

        #c# xr = self._run.GetCoord(axis=0)
        #c# yr = self._run.GetCoord(axis=1)

        #c# xmin = xr[0]
        #c# ymin = yr[0]

        # setup the initial condition

        self.r[:,:,self._it0] = self._r0
        self.v[:,:,self._it0] = self._v0

        # call super-fast cython functions now haha

        # the following function will move all the particles
        # for all time steps
        print 'calling _moveall'
        self._moveall()
        print 'calling _pfields'
        self._pfields()


    #==========================================================




    #==========================================================
    #==========================================================
    def _moveall(self):
        """move all particles backward and forward as necessary

        Exemple  : 

        Creation : 2013-05-04 15:40:32.257893

        """
        func          = self._lib.moveall
        func.restype  = None
        func.argtypes = [ndpointer(ctypes.c_double),    # pos
                         ndpointer(ctypes.c_double),    # vel
                         ndpointer(ctypes.c_float),     # E
                         ndpointer(ctypes.c_float),     # B
                         ctypes.c_uint,                 # nx
                         ctypes.c_uint,                 # ny
                         ctypes.c_double,               # dx
                         ctypes.c_double,               # dy
                         ctypes.c_double,               # xmin
                         ctypes.c_double,               # ymin
                         ctypes.c_double,               # dt
                         ctypes.c_double,               # charge
                         ctypes.c_double,               # mass
                         ctypes.c_uint,                 # nt
                         ctypes.c_uint,                 # it0
                         ctypes.c_uint]                 # npart


        #c# xc = self._run.GetCoord(axis=0)
        #c# yc = self._run.GetCoord(axis=1)

        xc = CR['xx'][0]
        yc = CR['yy'][0]
        self._dx = CR['xx'][1] - CR['xx'][0]
        self._dy = CR['yy'][1] - CR['yy'][0]

        func(self.r,
             self.v,
             self._E,
             self._B,
             self._B[0,:,:].shape[0],
             self._B[0,:,:].shape[1],
             self._dx,
             self._dy,
             xc, yc,
             self.dt,
             self._charge,
             self._mass,
             self._nt,
             self._it0,
             self._npart)
    #==========================================================




    #==========================================================
    #==========================================================
    def _pfields(self):
        """interpolates the fields at the particle positions

        Exemple  :

        Creation : 2013-05-04 15:40:32.257893

        """
        func          = self._lib.pfields
        func.restype  = None
        func.argtypes = [ndpointer(ctypes.c_double),    # pos
                         ndpointer(ctypes.c_float),     # E
                         ndpointer(ctypes.c_float),     # B
                         ctypes.c_uint,                 # nx
                         ctypes.c_uint,                 # ny
                         ctypes.c_uint,                 # npart
                         ctypes.c_uint,                 # nt
                         ctypes.c_double,               # xmin
                         ctypes.c_double,               # ymin
                         ctypes.c_double,               # dx
                         ctypes.c_double,               # dy
                         ndpointer(ctypes.c_double),    # Ep
                         ndpointer(ctypes.c_double)]    # Bp


        self._dx = CR['xx'][1] - CR['xx'][0]
        self._dy = CR['yy'][1] - CR['yy'][0]
        xc = CR['xx'][0]-self._dx/2.
        yc = CR['yy'][0]-self._dy/2.

        print 

        func(self.r,
             self._E,
             self._B,
             self._E[0,:,:].shape[0],
             self._E[0,:,:].shape[1],
             self._npart,
             self._nt,
             xc,yc,
             self._dx,
             self._dy,
             self.Ep,
             self.Bp)

    #==========================================================

