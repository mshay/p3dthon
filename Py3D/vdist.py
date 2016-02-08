#######################################################################
#                                                                     #
#                  Python Progs :  vdist.py                           #
#                  Aruthor      :  Colby Haggerty                     #
#                  Date         :  2016.02.08                         #
#                                                                     #
#                                                                     #
#######################################################################
import os
import sys 
import datetime
import numpy as np
import struct
import glob
import pdb
from scipy.io.idl import readsav
from scipy.ndimage import gaussian_filter

sys.exit()

class VDist(object):
    """ velocity distribution fucntion calculator
    """

    def __init__(self, ): 
        """ Initilazition Routine for the p3d_run object
        """

        self.parts = parts

    def _vdist_2d(self, v1, v2, bins=None):
        """ Simple 2D velocity histogram
        """

        H, xedges, yedges = np.histogram2d(velo[species]['v1'],velo[species]['v2'],**kwargs)
# H needs to be rotated and flipped
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(['V_1 vs V_2',H,xedges,yedges])

                H, xedges, yedges = np.histogram2d(velo[species]['v1'],velo[species]['v3'],**kwargs)
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(['V_1 vs V_3',H,xedges,yedges])

                H, xedges, yedges = np.histogram2d(velo[species]['v2'],velo[species]['v3'],**kwargs)
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(['V_2 vs V_3',H,xedges,yedges])
        else:
            for species in self.species:
                H, xedges, yedges = np.histogram2d(self.particles[species]['vx'],self.particles[species]['vy'],**kwargs)
# H needs to be rotated and flipped
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(['V_X vs V_Y',H,xedges,yedges])

                H, xedges, yedges = np.histogram2d(self.particles[species]['vx'],self.particles[species]['vz'],**kwargs)
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(['V_X vs V_Z',H,xedges,yedges])

                H, xedges, yedges = np.histogram2d(self.particles[species]['vy'],self.particles[species]['vz'],**kwargs)
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(['V_Y vs V_Z',H,xedges,yedges])

# Mask zeros
         #Hmasked = np.ma.masked_where(H==0,H)

        return return_hist_dict

    def vdist_2d( **kwargs):
        """ Generates differnent 2-Dimensional histograms for particles 
            Please god colby add a doc string!!!
        """

        if not kwargs.has_key('bins'): kwargs['bins']=51
        if par or Bvec or pitch or pizza:
            #qc#print 'Reading in the Fields form the Dump File'
            self.dump_field_dict = self.read_dump_file(fields=True)
            if pitch:
                return_hist = self._vdist_pitch(par=par,**kwargs)
            elif pizza:
                return_hist = self._vdist_pizza(par,**kwargs)
            elif par and not pitch and not pizza:
                return_hist = self._vdist_2d_par(**kwargs)
            else:
                return_hist = self._vdist_2d(**kwargs)
                #qc#print 'Interpolating the Bfield at the given r0 value'
                return_hist['B'] = np.array([
                   self.interp_field(self.dump_field_dict['bx']),
                   self.interp_field(self.dump_field_dict['by']), 
                   self.interp_field(self.dump_field_dict['bz'])])
        else:
            return_hist = self._vdist_2d(**kwargs)

        self._box = [self._r0[0]-self._dx[0]/2.,
                     self._r0[0]+self._dx[0]/2.,
                     self._r0[1]-self._dx[1]/2.,
                     self._r0[1]+self._dx[1]/2.]

        return_hist['box'] = self._box

        return return_hist



    def _vdist_pitch(self,par=False,energy=False,**kwargs): 
        """ Pitch party!!!!
        """
# I have a lot of code that can return a 2D distrubution function for a givn pitch angle range
# but I don't think that this will ever be usefull... I should remove it
        pa  = kwargs.pop('pa',90.)
        dpa = kwargs.pop('dpa',5.)
        wax = kwargs.pop('wax',0)
        v0_frame = kwargs.pop('v0_frame',None)

        if energy: par = True; kwargs['normed'] = True
        self._wax = wax
##############################################

        if par:
            b_interp = np.array([ self.interp_field(self.dump_field_dict['bx'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                                  self.interp_field(self.dump_field_dict['by'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                                  self.interp_field(self.dump_field_dict['bz'],self.param_dict['lx'],self.param_dict['ly'],self._r0)])

            e_interp = np.array([ self.interp_field(self.dump_field_dict['ex'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                                  self.interp_field(self.dump_field_dict['ey'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                                  self.interp_field(self.dump_field_dict['ez'],self.param_dict['lx'],self.param_dict['ly'],self._r0)])

            exb = np.cross(e_interp,b_interp)/sum(b_interp**2)
            if abs(np.sqrt(sum(exb**2))) < np.spacing(10):
                exb = np.cross(np.array([0.,1.,1.]),b_interp)
                exb = exb/np.sqrt(sum(exb**2))
            else:
                exb = exb/np.sqrt(sum(exb**2))

            b_interp = b_interp/np.sqrt(sum(b_interp**2))

            bxexb = np.cross(b_interp,exb) 
            #print b_interp
            #print exb
            #print bxexb

        velo={}
        return_hist_dict = {}
        self.subpart = {}
        #print 'Mucking this up!!'
        #Npart = np.shape(self.particles['e']['vx'])
        #self.particles['e']['vx'] = 1.*np.random.normal(scale=2.5,size=Npart)
        #self.particles['e']['vy'] = 1.*np.random.normal(scale=2.5,size=Npart)
        #self.particles['e']['vz'] = 1.*np.random.normal(scale=2.5,size=Npart)
        #print 'Done Mucking this up!!!'
        for species in self.species:
       
            if par:
                v0 = (b_interp[0]*self.particles[species]['vx']+
                      b_interp[1]*self.particles[species]['vy']+
                      b_interp[2]*self.particles[species]['vz'])

                v1 = (exb[0]*self.particles[species]['vx']+
                      exb[1]*self.particles[species]['vy']+
                      exb[2]*self.particles[species]['vz'])

                v2 = (bxexb[0]*self.particles[species]['vx']+
                      bxexb[1]*self.particles[species]['vy']+
                      bxexb[2]*self.particles[species]['vz'])
# Muck ===============
#               print 'Mucking this up!!!'
#               Npart = 1000000
#               v0 = 1.*np.random.normal(scale=np.sqrt(50.),size=Npart)
#               v1 = 1.*np.random.normal(scale=np.sqrt(50.),size=Npart)
#               v2 = 1.*np.random.normal(scale=np.sqrt(50.),size=Npart)
# Muck ===============
            else:
                v0 = self.particles[species]['vx']
                v1 = self.particles[species]['vy']
                v2 = self.particles[species]['vz']

            if wax == 0: # wax is just which axis defines the plane we are looking at
                vp0 = v1
                vp1 = v2
                vax = v0
            elif wax == 1:
                vp0 = v0
                vp1 = v1
                vax = v2
            elif wax == 2:
                vp0 = v0
                vp1 = v1
                vax = v2
            else :
                print ''
                print 'The plane axis is out of bounds fo wax = ',wax
                print 'vdist_2d is crashing!!!!'
                print ''
            
            if v0_frame:
                if energy:
                    vax = vax - np.mean(vax)
                    vp0 = vp0 - np.mean(vp0)
                    vp1 = vp1 - np.mean(vp1)
                else:
                    print '~'*80
                    print 'Shifting vax by ',np.mean(vax)
                    vax = vax - np.mean(vax)

            #self.vax = vax
            #self.vp = np.sqrt(vp0**2+vp1**2)
            #self.exb = exb
            #self.bxexb = bxexb
            # This means field aligned is 0, perp is 90 and anti aligned is 180
            pitch_angle = -1.0*(np.arctan(vax/np.sqrt(vp0**2+vp1**2))/np.pi*180. - 90.)
            self.ptc = pitch_angle

# Not sure which of these two is right
            #pitch_angle = np.arccos(vpar/vmag)/np.pi*180.
            #pitch_angle = np.arctan(np.sqrt((vmag**2 - vpar**2)/2.)/vpar)/np.pi*180. + 90.

            subpartind = np.where(abs(pitch_angle - pa) < dpa/2.)

            #subpart = self.particles[species][subpartind]
            #self.subpart[species] = subpart
            #subpitch_angle = pitch_angle[subpartind] 

            #print 'Total %s in pitch angle range are: %i'%(species,len(subpart))
            return_hist_dict[species] = []

#This is silly just do the whole distro right here
            if energy:
                if species == 'e':
                    reltv = True
                    if reltv:
                        v2 = (self.particles[species]['vx']**2+
                              self.particles[species]['vy']**2+
                              self.particles[species]['vz']**2)

                        KE = self.param_dict['m_e']* self.param_dict['c_2']*\
                             (1./np.sqrt(1. - v2/self.param_dict['c_2']) - 1.)

                    else:
                        KE = self.param_dict['m_e']/2.0*\
                            (self.particles[species]['vx']**2+
                             self.particles[species]['vy']**2+
                             self.particles[species]['vz']**2)
# Muck ===============
#                   print 'More Mucking!!!'
#                   KE = self.param_dict['m_e']/2.0*(v0**2+
#                                                    v1**2+
#                                                    v2**2)
#                   
# Muck ===============
                else:
                    KE = 1.0/2.0*(self.particles[species]['vx']**2+
                                  self.particles[species]['vy']**2+
                                  self.particles[species]['vz']**2)
# Muck ===============
#                   print 'More Mucking!!!'
#                   KE = 1./2.0*(v0**2+
#                                v1**2+
#                                v2**2)
#                   
# Muck ===============
                #H,xedges = np.histogram(KE,**kwargs)
                #return_hist_dict[species].append(H)
                #return_hist_dict[species].append((xedges[:-1]+xedges[1:])/2.0)


                print 'pa_size = ',np.size(pitch_angle)
                print 'KE_size = ',np.size(KE)
                H,xedges,yedges = np.histogram2d(pitch_angle,KE,**kwargs)
                xx,yy = np.meshgrid(yedges,xedges)
                eng = (xx[1:,1:]+xx[:-1,:-1])/2.
                ynorm = (np.cos(yy[:-1,:-1]/180.*np.pi) - np.cos(yy[1:,1:]/180.*np.pi))
                ynorm = ynorm/sum(ynorm)
                if species == 'e':
                    self.ynorm=ynorm
                    self.pa = pitch_angle
                    self.KE = KE
                    self.eng = eng
                    self.H = H
                    self.yyy = yy

# We need to account for a geometric factor because
# E is basicly v^2, so we are binning circles?

                #H = H/np.sqrt(eng)
                #return_hist_dict[species].append(H/ynorm*eng**2.*np.size(pitch_angle))

                Ebar = eng/self.param_dict['m_e']/self.param_dict['c_2']
                rel_vel = np.sqrt((Ebar**2 + 2.*Ebar)/(Ebar**2 + 2.*Ebar + 1.))
                rel_vel = rel_vel*np.sqrt(self.param_dict['c_2'])

                return_hist_dict[species].append(H/ynorm*eng*rel_vel*np.size(pitch_angle))

                return_hist_dict[species].append(xedges)
                return_hist_dict[species].append(yedges)

                
            else:
                #Colby devide by bin size!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                H, xedges, yedges = np.histogram2d(vp0[subpartind],
                                                   vp1[subpartind],
                                                   **kwargs)
# Muck ===============
# Mucking!!!!!
#                H = gaussian_filter(H, sigma=1.)
# Muck ===============

                xx,yy = np.meshgrid(yedges,xedges)
                print pa,dpa

                #int_cone = 1./18.*(
                self.inc = 1./18.*(
                       (-xx[:-1,:-1]**3+6.*xx[:-1,:-1]*yy[:-1,:-1]*np.sqrt(xx[:-1,:-1]**2+yy[:-1,:-1]**2) + 
                        3.*yy[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[:-1,:-1]**2) + xx[:-1,:-1] + np.spacing(4)) + 
                        3.*xx[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[:-1,:-1]**2)+yy[:-1,:-1] + np.spacing(4))) -
                       (-xx[1:,1:]**3+6.*xx[1:,1:]*yy[:-1,:-1]*np.sqrt(xx[1:,1:]**2+yy[:-1,:-1]**2) + 
                        3.*yy[:-1,:-1]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[:-1,:-1]**2) + xx[1:,1:] + np.spacing(4)) + 
                        3.*xx[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[:-1,:-1]**2)+yy[:-1,:-1] + np.spacing(4))) -
                       (-xx[:-1,:-1]**3+6.*xx[:-1,:-1]*yy[1:,1:]*np.sqrt(xx[:-1,:-1]**2+yy[1:,1:]**2) + 
                        3.*yy[1:,1:]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[1:,1:]**2) + xx[:-1,:-1] + np.spacing(4)) + 
                        3.*xx[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[1:,1:]**2)+yy[1:,1:] + np.spacing(4))) +
                       (-xx[1:,1:]**3+6.*xx[1:,1:]*yy[1:,1:]*np.sqrt(xx[1:,1:]**2+yy[1:,1:]**2) + 
                        3.*yy[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[1:,1:]**2) + xx[1:,1:] + np.spacing(4)) + 
                        3.*xx[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[1:,1:]**2)+yy[1:,1:] + np.spacing(4))))

                self.inc = self._int_cone(pa,dpa,xedges,yedges)
                #H = H/self.inc
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(H)
                return_hist_dict[species].append(xedges)
                return_hist_dict[species].append(yedges)


        return return_hist_dict

    def _vdist_pizza(self,
                     par=False,
                     z0=0.,
                     dz=1.,
                     wax=0,
                     v0_frame=False,
                     **kwargs): 
        """ Pizza Party!!!
            
            Also I am sooo sory if you ever have to
            use this code in this version... :(

            A little help:
            wax = which axis do we integarte over?
            wax = (0,1,2) -> (x,y,z) 
                or if in field line coordinates
            wax = (0,1,2) -> (||, ExB, bxExB)

        """

        self._z0  = z0
        self._dz  = dz
        self._wax = wax

        if par:
            b_interp = np.array([ 
                       self.interp_field(self.dump_field_dict['bx'],
                                         self.param_dict['lx'],
                                         self.param_dict['ly'],self._r0),
                       self.interp_field(self.dump_field_dict['by'],
                                         self.param_dict['lx'],
                                         self.param_dict['ly'],self._r0),
                       self.interp_field(self.dump_field_dict['bz'],
                                         self.param_dict['lx'],
                                         self.param_dict['ly'],self._r0)])

            e_interp = np.array([ 
                       self.interp_field(self.dump_field_dict['ex'],
                                         self.param_dict['lx'],
                                         self.param_dict['ly'],self._r0),
                       self.interp_field(self.dump_field_dict['ey'],
                                         self.param_dict['lx'],
                                         self.param_dict['ly'],self._r0),
                       self.interp_field(self.dump_field_dict['ez'],
                                         self.param_dict['lx'],
                                         self.param_dict['ly'],self._r0)])
            # Mike wanted me to do it this way
            fld = self.box_avg()

            #print 'Frist the interp fields:'
            #print 'B = ',b_interp
            #print 'E = ',e_interp
            b_interp = np.array([fld['bx'],fld['by'],fld['bz']])
            e_interp = np.array([fld['ex'],fld['ey'],fld['ez']])

            #print 'Then the avged fields:'
            #print 'B = ',b_interp
            #print 'E = ',e_interp

            exb = np.cross(e_interp,b_interp)/sum(b_interp**2)
            exb = exb/np.sqrt(sum(exb**2))

            b_interp = b_interp/np.sqrt(sum(b_interp**2))

            bxexb = np.cross(b_interp,exb) 
            #print b_interp
            #print exb
            #print bxexb

        velo={}
        return_hist_dict = {}
        self.subpart = {}
        #print 'Mucking this up!!'
        #Npart_2 = int(np.shape(self.particles['i']['vx'])[0]/2)

        #self.particles['i']['vx'][:Npart_2] = 1.*np.random.normal(scale=1.0,size=Npart_2)
        #self.particles['i']['vy'][:Npart_2] = 1.*np.random.normal(scale=1.0,size=Npart_2)
        #self.particles['i']['vz'][:Npart_2] = 1.*np.random.normal(scale=1.0,size=Npart_2)

        #self.particles['i']['vx'][Npart_2:] = 1.*np.random.normal(loc=10.,scale=1.0,size=Npart_2)
        #self.particles['i']['vy'][Npart_2:] = 1.*np.random.normal(loc=4.,scale=1.0,size=Npart_2)
        #self.particles['i']['vz'][Npart_2:] = 1.*np.random.normal(loc=4.,scale=1.0,size=Npart_2)
        ##print 'Done Mucking this up!!!'
        for species in self.species:
            return_hist_dict[species] = []
       
            if par:
# MuckMuckMuckMuck!!!!!
#                vxxx = self.particles[species]['vx'] + 1.0
# MuckMuckMuckMuck!!!!!

                #v0 = (b_interp[0]*vxxx+
                v0 = (b_interp[0]*self.particles[species]['vx']+
                      b_interp[1]*self.particles[species]['vy']+
                      b_interp[2]*self.particles[species]['vz'])

                #v1 = (exb[0]*vxxx+
                v1 = (exb[0]*self.particles[species]['vx']+
                      exb[1]*self.particles[species]['vy']+
                      exb[2]*self.particles[species]['vz'])

                #v2 = (bxexb[0]*vxxx+
                v2 = (bxexb[0]*self.particles[species]['vx']+
                      bxexb[1]*self.particles[species]['vy']+
                      bxexb[2]*self.particles[species]['vz'])
#Light Debuging
                self.bbb = b_interp
                self.exb = exb
                self.beb = bxexb
                #print 'Mucking this up!!!'
                #Npart = 100000
                #v0 = 1.*np.random.normal(scale=3.0,size=Npart)
                #v1 = 1.*np.random.normal(scale=1.0,size=Npart)
                #v2 = 1.*np.random.normal(scale=1.0,size=Npart)
            else:
                v0 = self.particles[species]['vx']
                v1 = self.particles[species]['vy']
                v2 = self.particles[species]['vz']

            if wax == 0: # wax is just which axis defines the plane we are looking at
                vp0 = v1
                vp1 = v2
                vax = v0
            elif wax == 1: 
                vp0 = v2
                vp1 = v0
                vax = v1
            elif wax == 2:
                vp0 = v0
                vp1 = v1
                vax = v2
            else :
                print ''
                print 'The plane axis is out of bounds fo wax = ',wax
                print 'vdist_2d is crashing!!!!'
                print ''

            if v0_frame:
                print '#'*80
                print 'Shifting vax by ',np.mean(vax)
                print '#'*80
                vax = vax - np.mean(vax)

            subpartind = np.where(abs(vax - z0) < dz/2.)
            print ''
            print 'The reduced number of {0} is {1}'.format(species,np.size(subpartind))
            print ''

            H, xedges, yedges = np.histogram2d(vp0[subpartind],vp1[subpartind],**kwargs)

            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(H)
            return_hist_dict[species].append(xedges)
            return_hist_dict[species].append(yedges)


        return return_hist_dict




    def _vdist_2d_par(self,**kwargs):
        """
        I dont know what is better
            To rotate all particles to a particular b feild
            or rotate all particles to their own b field?

            right now each part has its own b field
        """
        
        b_interp = np.array([ self.interp_field(self.dump_field_dict['bx'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                           self.interp_field(self.dump_field_dict['by'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                           self.interp_field(self.dump_field_dict['bz'],self.param_dict['lx'],self.param_dict['ly'],self._r0)])
        b_interp = b_interp/np.sqrt(sum(b_interp**2))

        e_interp = np.array([ self.interp_field(self.dump_field_dict['ex'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                           self.interp_field(self.dump_field_dict['ey'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                           self.interp_field(self.dump_field_dict['ez'],self.param_dict['lx'],self.param_dict['ly'],self._r0)])

        print 'b_interp'
        print b_interp
        self.trash1 = b_interp
        print 'e_interp'
        print e_interp
        self.trash2 = e_interp

        exb = np.cross(e_interp,b_interp) 
        exb = exb/np.sqrt(sum(exb**2))

        bxexb = np.cross(b_interp,exb) 

        if abs(sum(bxexb**2) - 1.0) > .001:
            print '###################### WARNING ######################'
            print '      your b and exb seem to not be perpandicular!!!'
            bxexb = bxexb/np.sqrt(sum(bxexb**2))

        #qc#print 'Rotating Velocties'
        velo={}
        for species in self.species:
            velo[species] = {}

            rotate_all_parts = False
            if rotate_all_parts == True:

                delx = self.param_dict['lx']*1.0/self.param_dict['pex']/self.param_dict['nx']
                dely = self.param_dict['ly']*1.0/(self.param_dict['pey']*self.param_dict['ny'])

                xind = (np.floor((self.particles[species]['x']-delx/2.0)/delx)).astype(int)
                yind = (np.floor((self.particles[species]['y']-dely/2.0)/dely)).astype(int)

                wx = (self.particles[species]['x']-delx/2.0)%delx
                wy = (self.particles[species]['y']-dely/2.0)%dely

                partbx = wx     *wy     *self.dump_field_dict['bx'][xind.tolist()    ,yind.tolist()] + \
                         (1.-wx)*wy     *self.dump_field_dict['bx'][(xind+1).tolist(),yind.tolist()] + \
                         wx     *(1.-wy)*self.dump_field_dict['bx'][xind.tolist()    ,(yind+1).tolist()] + \
                         (1.-wx)*(1.-wy)*self.dump_field_dict['bx'][(xind+1).tolist(),(yind+1).tolist()] 

                partby = wx     *wy     *self.dump_field_dict['by'][xind.tolist()    ,yind.tolist()] + \
                         (1.-wx)*wy     *self.dump_field_dict['by'][(xind+1).tolist(),yind.tolist()] + \
                         wx     *(1.-wy)*self.dump_field_dict['by'][xind.tolist()    ,(yind+1).tolist()] + \
                         (1.-wx)*(1.-wy)*self.dump_field_dict['by'][(xind+1).tolist(),(yind+1).tolist()] 

                partbz = wx     *wy     *self.dump_field_dict['bz'][xind.tolist()    ,yind.tolist()] + \
                         (1.-wx)*wy     *self.dump_field_dict['bz'][(xind+1).tolist(),yind.tolist()] + \
                         wx     *(1.-wy)*self.dump_field_dict['bz'][xind.tolist()    ,(yind+1).tolist()] + \
                         (1.-wx)*(1.-wy)*self.dump_field_dict['bz'][(xind+1).tolist(),(yind+1).tolist()] 

                (partbx,partby,partbz) = (partbx/np.sqrt(partbx**2+partby**2+partbz**2),
                                          partby/np.sqrt(partbx**2+partby**2+partbz**2),
                                          partbz/np.sqrt(partbx**2+partby**2+partbz**2))

                partpb1 =  0.0*partbx
                partpb2 = -1.*np.sign(by_interp)*partbz/(partbz**2 + partby**2)**(.5)
                partpb3 = 1.*np.sign(by_interp)*partby/(partbx**2 + partby**2)**(.5)

                (partpb1,partpb2,partpb3) = (partpb1/np.sqrt(partpb1**2+partpb2**2+partpb3**2),
                                             partpb2/np.sqrt(partpb1**2+partpb2**2+partpb3**2),
                                             partpb3/np.sqrt(partpb1**2+partpb2**2+partpb3**2))

                partpp1 =  (partby*partpb3 - partbz*partpb2)
                partpp2 =  (partbz*partpb1 - partbx*partpb3)
                partpp3 =  (partbx*partpb2 - partby*partpb1)

                velo[species]['par']   = (partbx*self.particles[species]['vx']+partby*self.particles[species]['vy']+partbz*self.particles[species]['vz'])
                velo[species]['perp1'] = (partpb1*self.particles[species]['vx']+partpb2*self.particles[species]['vy']+partpb3*self.particles[species]['vz'])
                velo[species]['perp2'] = (partpp1*self.particles[species]['vx']+partpp2*self.particles[species]['vy']+partpp3*self.particles[species]['vz'])

            else:
                velo[species]['par']   = (b_interp[0]*self.particles[species]['vx']+
                                          b_interp[1]*self.particles[species]['vy']+
                                          b_interp[2]*self.particles[species]['vz'])

                velo[species]['perp1'] = (exb[0]*self.particles[species]['vx']+
                                          exb[1]*self.particles[species]['vy']+
                                          exb[2]*self.particles[species]['vz'])

                velo[species]['perp2'] = (bxexb[0]*self.particles[species]['vx']+
                                          bxexb[1]*self.particles[species]['vy']+
                                          bxexb[2]*self.particles[species]['vz'])


        #velo['e']['par']   = (bx_interp*self.particles['e']['vx']+by_interp*self.particles['e']['vy']+bz_interp*self.particles['e']['vz'])/bmag_interp
        #velo['e']['perp1'] = self.particles['e']['vx']*b_perp1x+self.particles['e']['vy']*b_perp1y+self.particles['e']['vz']*b_perp1z
        #velo['e']['perp2'] = self.particles['e']['vx']*b_perp2x+self.particles['e']['vy']*b_perp2y+self.particles['e']['vz']*b_perp2z

        #qc#print 'Generating Histograms'
>>>>>>> f1535e5f06cfecfdc6dafa3bb35419013221891d

        return_hist_dict = {}
        for species in self.species:
            return_hist_dict[species] = []

        for species in velo.keys():
            H, xedges, yedges = np.histogram2d(velo[species]['par'],velo[species]['perp1'],**kwargs)
# H needs to be rotated and flipped
            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(['V_b vs V_exb',H,xedges,yedges])

            H, xedges, yedges = np.histogram2d(velo[species]['par'],velo[species]['perp2'],**kwargs)
            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(['V_b vs V_bxexb 2',H,xedges,yedges])

            H, xedges, yedges = np.histogram2d(velo[species]['perp1'],velo[species]['perp2'],**kwargs)
            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(['V_exb 1 vs V_bxexb',H,xedges,yedges])

# Mask zeros
         #Hmasked = np.ma.masked_where(H==0,H)

        return return_hist_dict



    #def get_part_in_box([location, width]):
    def _part_in_box(self,r0=[0.5,0.5],dx=[1.,1.]):
        """
        #--------------------------------------------------------------
        #   Method      : get_part_in_box
        #
        #   Discription : This method accepts a point and a width in 
        #                 simulation units (c/wpi) to define a box.
        #               : In that box we bin all of the particles to 
                          form the effective distrobution function
        #
        #   Args     r0 : location [x,y] ( where you want the center 
                          of you box to be located at)
        #            dx : width [x,y] (the width of the box to bin 
                           particles 
        #               : dump_num (This spesifies the particular 
                          runs dump file 
        #
        #   Comments    : It would be pretty easy and potential 
                          usefull to allow this to wrap around 
                          the edges
        #               : so we can pick zero as a boundry and the 
                          code will know what todo.
                          #--------------------------------------------------------------
        """
        x0 = r0[0]
        y0 = r0[1]
        if isinstance(dx,float) or isinstance(dx,float):
            dx = dx*1.0
            dy = dx       # Square box is assumed
        else:
            dy = dx[1]
            dx = dx[0]

        #qc#print 'r0 = [%f,%f] and dx = [%f,%f]'%(x0,y0,dx,dy)
# Figure out which set of processors we are on
        xlb = x0 - dx/2.
        xub = x0 + dx/2.
        ylb = y0 - dy/2.
        yub = y0 + dy/2.

# BIGNOTE: This seems iffy you should come back and double check this colby. It doesnt make sence that we should have to add 1
        #xproc_lb = (int(np.floor(1.0*self.param_dict['pex']*xlb/self.param_dict['lx']))+1)%int(self.param_dict['nchannels'])
        #xproc_lb = (int(np.floor(1.0*self.param_dict['pex']*xlb/self.param_dict['lx']))+1)%int(self.param_dict['nchannels'])
        if type(self.param_dict['nchannels']) == str:
# Man you should fix this colby!!!
            #self.param_dict['nchannels'] =  self.param_dict['pex'] # a lot of times we just set nchannels as pex
            self.param_dict['nchannels'] =  self.param_dict[self.param_dict['nchannels']]
            
        
        #xproc_lb = (int(np.floor(1.0*self.param_dict['pex']*xlb/self.param_dict['lx']))+1)%self.param_dict['nchannels']
# We are chaning this but it may not work for all stuff so you know 
        xproc_lb = (int(np.floor(1.0*self.param_dict['pex']*xlb/self.param_dict['lx'])))%self.param_dict['nchannels'] +1
# Colby this needs to be cooded better!!!
# if you have a diffent number of channels than pex you run into some shit
        yproc_lb = int(np.floor(1.0*self.param_dict['pey']*ylb/self.param_dict['ly']))
        #yproc_lb = int(np.floor(1.0*self.param_dict['pey']*ylb/self.param_dict['ly']))*2 #Uncomment for yishin
        #xproc_ub = (int(np.floor(1.0*self.param_dict['pex']*xub/self.param_dict['lx']))+1)%int(self.param_dict['nchannels'])


        if abs(xub - self.param_dict['lx']) < abs(np.spacing(2)): 
            xub = self.param_dict['lx'] - np.spacing(2)
        xproc_ub = (int(np.floor(1.0*self.param_dict['pex']*xub/self.param_dict['lx'])))%self.param_dict['nchannels'] +1
        yproc_ub = int(np.floor(1.0*self.param_dict['pey']*yub/self.param_dict['ly'])) 
        #yproc_ub = int(np.floor(1.0*self.param_dict['pey']*yub/self.param_dict['ly']))*2 #Uncomment for yishin

        if xproc_lb > xproc_ub:
            print 'Lower Bound greater than upper bound! That is obviously an issue!'
            return -1

        if xproc_lb < 1: xproc_lb = 1
        if xproc_ub > self.param_dict['nchannels']: xproc_ub = self.param_dict['nchannels'] 

#Colby! you can code this smarter but presently you are not doing that!
# To code smarter, insted of making this simple upper bound you should code
# to allow for nchanles to be a non multiple of pex
        max_yproc =  int(round(self.param_dict['pex']*self.param_dict['pey']/self.param_dict['nchannels']))
        if yproc_lb < 0: yproc_lb = 0
        if yproc_ub > max_yproc:
            yproc_ub = max_yproc


        xprocs_2load = range(xproc_lb,xproc_ub+1)
        yprocs_2load = range(yproc_lb,yproc_ub+1)
        
        print 'The x processors we will be loading (i.e. the dump files) are: {}'.format(xprocs_2load)
        print 'The y processors we will be loading (i.e. the sub arrays) are: {}'.format(yprocs_2load)
        #c# print 'That means we have {} processors to load and you can expect an aprox. {:.2f} min wait time'.format(
        #c#      len(xprocs_2load)*len(yprocs_2load),1.41*len(xprocs_2load)*len(yprocs_2load)*4./60.)

# Load in the appropriate Processors
        temp_dump_pruned = {}# List to hold all the particles
        for species in self.species:
            temp_dump_pruned[species] = []
        #first for loop over px
        for xprocs_index in xprocs_2load:
            dump_dat_dict = {}
            dump_index = self._num_to_ext(xprocs_index)
            #qc#print 'Loading in xproc number '+dump_index
            temp_dump_dat = self.read_dump_file(dump_index)

            #c# new_tdd = {}
# We need to#c#  throw out the extra data
            #c# if (int(np.floor(1.0*self.param_dict['pex']*xlb/self.param_dict['lx']))+1)/self.param_dict['nchannels'] < 1:
            #c#     new_tdd['i'] = temp_dump_dat['i'][:len(temp_dump_dat['i'])/2-1]
            #c#     new_tdd['e'] = temp_dump_dat['e'][:len(temp_dump_dat['e'])/2-1]
            #c# else:
            #c#     new_tdd['i'] = temp_dump_dat['i'][len(temp_dump_dat['i'])/2:]
            #c#     new_tdd['e'] = temp_dump_dat['e'][len(temp_dump_dat['e'])/2:]

            #c# temp_dump_dat = new_tdd



# We need to looop over Ions and Electrons
            for species in self.species:
                #if species == 'i':
                #    print '\tSelecting Ions'
                #else:
                #    print '\tSelecting Electrons'
# second for loop over the py
# also just doing electron for right now
                for yprocs_index in yprocs_2load:
                    self._debug = temp_dump_dat 
                    temp_dump_yproc = temp_dump_dat[species][yprocs_index]
### Lets try somthing new, that might be faster.  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    temp_index = np.where((temp_dump_dat[species][yprocs_index]['y'] - yub)**2 +
                    (temp_dump_dat[species][yprocs_index]['y'] - ylb)**2 < (yub-ylb)**2 )
                    temp_dump_xproc = temp_dump_dat[species][yprocs_index][temp_index]

                    temp_index = np.where((temp_dump_xproc['x'] - xub)**2 +
                    (temp_dump_xproc['x'] - xlb)**2 < (xub-xlb)**2 )

                    temp_dump_pruned[species].append(temp_dump_xproc[temp_index])


#test ### This is how we were sorting this But I think it is MUCH SLOWER
#test # You only need to sort if you are on the edge processors
#test                     if yprocs_index == yprocs_2load[0] or yprocs_index == yprocs_2load[-1]: 
#test                         #qc#print '\t\tSorting yproc number '+str(yprocs_index)
#test                         sorted_index = temp_dump_dat[species][yprocs_index].argsort(order='y')
#test                         temp_dump_yproc = temp_dump_dat[species][yprocs_index][sorted_index]
#test # Here we need kind of a complecated if structure to get all the poible cases since
#test # we are scaning over muliple processors.
#test # If you are on your first y processor then you need to find a lower boundry
#test                     if yprocs_index == yprocs_2load[0]: 
#test                         #qc#print '\t\t\tFinding lower yboundry index '
#test                         lower_yboundry_index = np.searchsorted(temp_dump_yproc['y'],ylb)
#test                     else:
#test                         lower_yboundry_index = 0#np.searchsorted(temp_dump_yproc['y'],ylb)
#test # If you are on your last y processor then you need to find a upper boundry
#test                     if yprocs_index == yprocs_2load[-1]: 
#test                         #qc#print '\t\t\tFinding upper yboundry index '
#test                         upper_yboundry_index = np.searchsorted(temp_dump_yproc['y'],yub)
#test                     else:
#test                         upper_yboundry_index = -1#np.searchsorted(temp_dump_yproc['y'],yub)
#test                     # You only need to sort if you are on the edge processors
#test                     temp_dump_xproc = temp_dump_yproc[lower_yboundry_index:upper_yboundry_index]
#test                     if xprocs_index == xprocs_2load[0] or xprocs_index == xprocs_2load[-1]: 
#test                         #qc#print '\t\tNow sorting x values for remaing data'
#test                         sorted_index = temp_dump_xproc.argsort(order='x')
#test                         temp_dump_xproc = temp_dump_xproc[sorted_index] 
#test # If you are on your first x processor then you need to find a lower boundry
#test                     if xprocs_index == xprocs_2load[0]: 
#test                         #qc#print '\t\t\tFinding lower xboundry index '
#test                         lower_xboundry_index = np.searchsorted(temp_dump_xproc['x'],xlb)
#test                     else:
#test                         lower_xboundry_index = 0#np.searchsorted(temp_dump_xproc['x'],xlb)
#test # If you are on your last x processor then you need to find a upper boundry
#test                     if xprocs_index == xprocs_2load[-1]: 
#test                         #qc#print '\t\t\tFinding upper xboundry index '
#test                         upper_xboundry_index = np.searchsorted(temp_dump_xproc['x'],xub)
#test                     else: 
#test                         upper_xboundry_index = -1#np.searchsorted(temp_dump_xproc['x'],xub)
#test                     temp_dump_pruned[species].append(temp_dump_xproc[lower_xboundry_index:upper_xboundry_index])
        for species in self.species:
            temp_dump_pruned[species] = np.concatenate(temp_dump_pruned[species])
            print 'Total %s in box are: %i'%(species,len(temp_dump_pruned[species]))
        return temp_dump_pruned



    def _num_to_ext(self,num):
        if type(num) is str: 
            return num
        else:
            if int(np.floor(num / 10)) > 0: 
                if int(np.floor(num / 100)) > 0: 
                    return str(num)
                else:
                    return '0'+str(num)
            else:
                return '00'+str(num)

    def _int_cone(self,pa,dpa,xedges,yedges):
        """
#---------------------------------------------------------------------------------------------------------------
#   Method      : _int_cone
#
#   Discription : This integrates a cone over a differental grid
#
#   Args        : xedges the x edges of the histogram
#               : yedges the y edges of the histogram
#
#   Comments    : I think this is working ok? 
#---------------------------------------------------------------------------------------------------------------
        """
        

        xx,yy = np.meshgrid(yedges,xedges)
        intcone = 1./18.*(
               (-xx[:-1,:-1]**3+6.*xx[:-1,:-1]*yy[:-1,:-1]*np.sqrt(xx[:-1,:-1]**2+yy[:-1,:-1]**2) + 
                3.*yy[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[:-1,:-1]**2) + xx[:-1,:-1] + np.spacing(4)) + 
                3.*xx[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[:-1,:-1]**2)+yy[:-1,:-1] + np.spacing(4))) -
               (-xx[1:,1:]**3+6.*xx[1:,1:]*yy[:-1,:-1]*np.sqrt(xx[1:,1:]**2+yy[:-1,:-1]**2) + 
                3.*yy[:-1,:-1]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[:-1,:-1]**2) + xx[1:,1:] + np.spacing(4)) + 
                3.*xx[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[:-1,:-1]**2)+yy[:-1,:-1] + np.spacing(4))) -
               (-xx[:-1,:-1]**3+6.*xx[:-1,:-1]*yy[1:,1:]*np.sqrt(xx[:-1,:-1]**2+yy[1:,1:]**2) + 
                3.*yy[1:,1:]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[1:,1:]**2) + xx[:-1,:-1] + np.spacing(4)) + 
                3.*xx[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[1:,1:]**2)+yy[1:,1:] + np.spacing(4))) +
               (-xx[1:,1:]**3+6.*xx[1:,1:]*yy[1:,1:]*np.sqrt(xx[1:,1:]**2+yy[1:,1:]**2) + 
                3.*yy[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[1:,1:]**2) + xx[1:,1:] + np.spacing(4)) + 
                3.*xx[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[1:,1:]**2)+yy[1:,1:] + np.spacing(4))))

        norm = 1./18.*(
               (-xx[0,0]**3+6.*xx[0,0]*yy[0,0]*np.sqrt(xx[0,0]**2+yy[0,0]**2) + 
                3.*yy[0,0]**3*np.log(np.sqrt(xx[0,0]**2+yy[0,0]**2) + xx[0,0] + np.spacing(4)) + 
                3.*xx[0,0]**3*np.log(np.sqrt(xx[0,0]**2+yy[0,0]**2)+yy[0,0] + np.spacing(4))) -
               (-xx[-1,-1]**3+6.*xx[-1,-1]*yy[0,0]*np.sqrt(xx[-1,-1]**2+yy[0,0]**2) + 
                3.*yy[0,0]**3*np.log(np.sqrt(xx[-1,-1]**2+yy[0,0]**2) + xx[-1,-1] + np.spacing(4)) + 
                3.*xx[-1,-1]**3*np.log(np.sqrt(xx[-1,-1]**2+yy[0,0]**2)+yy[0,0] + np.spacing(4))) -
               (-xx[0,0]**3+6.*xx[0,0]*yy[-1,-1]*np.sqrt(xx[0,0]**2+yy[-1,-1]**2) + 
                3.*yy[-1,-1]**3*np.log(np.sqrt(xx[0,0]**2+yy[-1,-1]**2) + xx[0,0] + np.spacing(4)) + 
                3.*xx[0,0]**3*np.log(np.sqrt(xx[0,0]**2+yy[-1,-1]**2)+yy[-1,-1] + np.spacing(4))) +
               (-xx[-1,-1]**3+6.*xx[-1,-1]*yy[-1,-1]*np.sqrt(xx[-1,-1]**2+yy[-1,-1]**2) + 
                3.*yy[-1,-1]**3*np.log(np.sqrt(xx[-1,-1]**2+yy[-1,-1]**2) + xx[-1,-1] + np.spacing(4)) + 
                3.*xx[-1,-1]**3*np.log(np.sqrt(xx[-1,-1]**2+yy[-1,-1]**2)+yy[-1,-1] + np.spacing(4))))
        
        self.normie=norm

        print 'norm = ',norm
        print 'otha = ',abs(1./np.tan((pa-dpa/2.)/180.*np.pi) - 1./np.tan((pa+dpa/2.)/180.*np.pi))

        #return intcone/norm#*abs(1./np.tan((pa-dpa/2.)/180.*np.pi) - 1./np.tan((pa+dpa/2.)/180.*np.pi))
        return intcone/intcone.min()#*abs(1./np.tan((pa-dpa/2.)/180.*np.pi) - 1./np.tan((pa+dpa/2.)/180.*np.pi))

# Ripped this code from else where nee to add it
    def box_avg(self):
        
        if not hasattr(self, 'dump_field_dict'):
            print 'Field data not loaded somthings wrong!!!!'
            return None
        
        r0 = self._r0
        dx = self._dx

        x0 = r0[0] - dx[0]/2.
        x1 = r0[0] + dx[0]/2.
        y0 = r0[1] - dx[1]/2.
        y1 = r0[1] + dx[1]/2.

        dx0 = 1.0*self.param_dict['lx']/(self.param_dict['nx']*self.param_dict['pex'])

        ip0 = int(np.ceil((x0 -  dx0/2.)/dx0))
        ip1 = int(np.floor((x1 - dx0/2.)/dx0))
        jp0 = int(np.ceil((y0 -  dx0/2.)/dx0))
        jp1 = int(np.floor((y1 - dx0/2.)/dx0))

        if ip1 < ip0 or jp1 < jp0:
            print 'No grid points found in box!!!!'
            return None
        #else:
        #    print 'DEBUG: ',ip0,ip1,jp0,jp1
       
        iprg, jprg = [],[]
        for c in range(ip0,ip1+1):
            for d in range(jp0,jp1+1):
                iprg.append(c)
                jprg.append(d)

        avg_flds = {}
        for k,fld in self.dump_field_dict.iteritems():
            avg_flds[k] = np.mean(fld[jprg,iprg])

        return avg_flds


    def _shift_frame(self):
        """ Shifts all the particles velocties into their fluid frame.
            This finds the mean ux uy and uz for e and i and shifts everying
            by thoes values
        """

        for s in self.species:
            for c in ['vx','vy','vz']:
                uu = np.mean(self.particles[s][c])
                print '{0}_{1} = {2}'.format(c,s,uu)
                self.particles[s][c] = self.particles[s][c] - uu


    def interp_field(self,field,lx=None,ly=None,r0=None):
        """
#----------------------------------------------------------------------
#        
#   Method      : interp_field
#
#   Discription : This method takes a field and a floating point, and returns the linear fit value 
#               : between the grid points
#
#   Args        : field  The field you are interpolating
#               : r0[0] The xpoint to interpolate at
#               : r0[1] The ypoint to interpolate at
#
#   Comments    : I think this is working ok? It would be smart to make
#                 this an object method that just reads
#               : the internally saved field. so CODE IN THE FUTURE
#----------------------------------------------------------------------
        """
        if lx is None:lx=self.param_dict['lx']
        if ly is None:ly=self.param_dict['ly']
        if r0 is None:r0=self._r0

        nx = len(field[0,:])
        ny = len(field[:,0])
        ip = int(np.floor(1.0*r0[0]/lx*nx))
        jp = int(np.floor(1.0*r0[1]/ly*ny))

        if ip + 1 > nx-1: ipp1 = 0
        else: ipp1 = ip+1

        if jp + 1 > ny-1: jpp1 = 0
        else: jpp1 = jp+1

       # print ip,jp,nx,ny

        wx = 1.0*r0[0]/lx*nx - np.floor(1.0*r0[0]/lx*nx)
        wy = 1.0*r0[1]/ly*ny - np.floor(1.0*r0[1]/ly*ny)

        return (1.-wx)*(1.-wy)*field[jp  , ip]+\
               (wx)   *(1.-wy)*field[jpp1, ip]+\
               (1.-wx)*   (wy)*field[jp  ,ipp1]+\
               (wx)   *   (wy)*field[jpp1,ipp1]
    


#Orphaned method that I use in load param but not quite sure how to fit it in
def convert(val):
    constructors = [int, float, str]
    for c in constructors:
        try:
            return c(val)
        except ValueError:
            pass
            

#c# This was a set of code I used to rotate each particle indviduly
#c# Mike insists this is silly, so we will not be using this
#c#            rotate_all_parts = False
#c#            if rotate_all_parts == True:
#c#
#c#                xind = (np.floor((self.particles[species]['x']-delx/2.0)/delx)).astype(int)
#c#                yind = (np.floor((self.particles[species]['y']-dely/2.0)/dely)).astype(int)
#c#
#c#                wx = (self.particles[species]['x']-delx/2.0)%delx
#c#                wy = (self.particles[species]['y']-dely/2.0)%dely
#c#
#c#                (1.-wx)*(1.-wy)
#c#                (1.-wx)*wy     
#c#                wx     *(1.-wy)
#c#                wx     *wy     
#c#
#c#                partbx = (1.-wx)*(1.-wy)*self.dump_field_dict['bx'][xind.tolist()    ,yind.tolist()] + \
#c#                         (1.-wx)*wy     *self.dump_field_dict['bx'][(xind+1).tolist(),yind.tolist()] + \
#c#                         wx     *(1.-wy)*self.dump_field_dict['bx'][xind.tolist()    ,(yind+1).tolist()] + \
#c#                         wx     *wy     *self.dump_field_dict['bx'][(xind+1).tolist(),(yind+1).tolist()] 
#c#
#c#                partby = (1.-wx)*(1.-wy)*self.dump_field_dict['by'][xind.tolist()    ,yind.tolist()] + \
#c#                         (1.-wx)*wy     *self.dump_field_dict['by'][(xind+1).tolist(),yind.tolist()] + \
#c#                         wx     *(1.-wy)*self.dump_field_dict['by'][xind.tolist()    ,(yind+1).tolist()] + \
#c#                         wx     *wy     *self.dump_field_dict['by'][(xind+1).tolist(),(yind+1).tolist()] 
#c#
#c#                partbz = (1.-wx)*(1.-wy)*self.dump_field_dict['bz'][xind.tolist()    ,yind.tolist()] + \
#c#                         (1.-wx)*wy     *self.dump_field_dict['bz'][(xind+1).tolist(),yind.tolist()] + \
#c#                         wx     *(1.-wy)*self.dump_field_dict['bz'][xind.tolist()    ,(yind+1).tolist()] + \
#c#                         wx     *wy     *self.dump_field_dict['bz'][(xind+1).tolist(),(yind+1).tolist()] 
#c#
#c#                vmag = np.sqrt(self.particles[species]['vx']**2+self.particles[species]['vy']**2+self.particles[species]['vz']**2)
#c#                vpar = (self.particles[species]['vx']*partbx+self.particles[species]['vy']*partby+self.particles[species]['vz']*partbz)/np.sqrt(partbx**2+partby**2+partbz**2)



        

