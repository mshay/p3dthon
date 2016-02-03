import os
import pdb
import glob
import numpy as np
from ._methods import load_param
from ._methods import _num_to_ext

class Movie(object):
    """Class to load p3d movie data"""

    def __init__(self,
                 num=None,
                 param=None,
                 path='./'):
        """ Initlize a movie object
        """

        self.path       = self._get_movie_path(path)
        self.param      = load_param(param)
        self.num        = self._get_movie_num(num)
        self.movie_vars = self._get_movie_vars()
        self.log        = self._load_log()
        self.ntimes     = len(self.log[self.movie_vars[0]])

    def get_fields(self):
        pass

#    def load_movie(self,var='NOT_SET',time=None,local=False):
#        """ A method to load a particular value for a given time
#
#        Load movie files for a given set of varibles.
#        You can pass as a list, or a single string, or a keyword all
#
#        @return: @todo
#
#        Exemple  :
#
#        Creation
#        :
#        2014-06-16
#        """
#
#        if type(var) is not list:
#            if var.lower() == 'all': var_arr = self.movie_arr
#            else: var_arr = [var]
#        else: var_arr = var
#        return_dict = {}
#        for cosa in var_arr:
#            if (cosa not in self.movie_arr):
#                print 'Varable %s not found in movie_arr. Nothing was loaded!'%cosa
#                cosa = raw_input('Please Enter a Varible: ')
#            if (cosa not in self.movie_arr):
#                print 'Varable %s not found in movie_arr. Nothing was loaded!'%cosa
#                print 'You dont get a second try!'
#                return -1
#
#            if time is None:
#                time = raw_input('Time %s out of range [0 - %i]\n'% \
#                                 (time,self.num_of_times-1) + \
#                                 'Please Enter a time: ')
#
#            if time == 'all':
#                time = range(self.num_of_times)
#            elif type(time) == str:
#                time = [int(x) for x in time.split()]
#            elif type(time) is int:
#                time = [time]
#            else:
#                time = list(time)
#
#            for chose in time:
#                if -1 < chose < self.num_of_times:
#                    pass
#                else:
#                    print 'Time %i is out of time range [%i - %i]'\
#                          %(chose,0,self.num_of_times-1)
#                    return None
#
#            fname = self.movie_path+'/movie.'+cosa+'.'+self.movie_num_str
#            fname = os.path.abspath(fname)

    def _read_movie(self,var, time):
        
        # Insert Comment about werid movie shape
        pt = self.ntimes
        pz = self.param['pez']*self.param['nz'] 
        py = self.param['pey']*self.param['ny'] 
        px = self.param['pex']*self.param['nx'] 
        
        movie_shape = (pt,pz,py,px)

        fname = self.path+'/movie.'+var+'.'+self.num
        print "Loading {0}".format(fname)

        # It seems that Marc Swisdak hates us and wants to be unhappy because 
        # the byte data is unsigned and the doulbe byte is signed so that is 
        # why one has a uint and the other is just int
        if 'double_byte' in self.param:
            dat_type = np.dtype('int16')
            norm_cst = 256**2-1
            shft_cst = 1.0*256**2/2
        else: #single byte precision
            dat_type = np.dtype('uint8')
            norm_cst = 256-1
            shft_cst = 0.0

        t = 0 # This keeps track of where we are
        cmin = self.log[var][:,0][time]
        cmax = self.log[var][:,1][time]

        fp = np.memmap(fname, dtype=dat_type, mode='r', shape=movie_shape)

        fp = fp[time].copy()
        
#        pdb.set_trace()

        return cmin,cmax,fp

            
#            with open(fname,'r') as F:
#                for c,t in enumerate(time):
#
#                    F.seek((chose - _)*grid_pts*dat_type.itemsize, os.SEEK_SET)
#                    byte_arr[i,:,:] = np.fromfile(f,
#                                                  dtype=dat_type,
#                                                  count=grid_pts
#                                                  ).reshape(ney,nex)
#
#                    byte_arr[i,:,:] = (1.0*byte_arr[i,:,:] + shft_cst)* \
#                                      (lmax[chose]-lmin[chose]) \
#                                      /(1.0*norm_cst) + lmin[chose]
#
#                    _ = chose + 1
#
#            return_dict[cosa] = np.squeeze(byte_arr)
#
#        return_dict.update(zip(('xx','yy'), self.get_domain_arrays()))
#        return return_dict



    def _load_log(self):
        """ Loads the log file for a given set of moives files
            It creates a dictoary 
        """

        fname = self.path+'/movie.log.'+self.num

        print "Loading {0}".format(fname)
        clims = np.loadtxt(fname)
        
        if len(clims)%len(self.movie_vars) != 0:
            raise Exception('Param/Moive Incompatibility')

        log = {}
        for c,k in enumerate(self.movie_vars):
            log[k] = clims[c::30,:]
    
        return log
        # usefull use later
        #print "movie.log '%s' has %i time slices"%(fname,self.num_of_times)

# STRUCTURE OF movie_log_dict{}
#   movie_log_dict is a dictionary of all the of the varibles that could be read in a movie file
#   you pass the standered name of the varible as a string and you get back an array.
#   in the array each element coresponds to a diffrent time slice
#   so      movie.movie_log_dict['bz'] = [

    def _get_movie_path(self,path):

        attempt_tol = 5
        path = os.path.abspath(path)
        choices = glob.glob(path+'/movie.log.*')

        c = 0
        while not choices and c < attempt_tol:
            print '='*20 + ' No movie files found ' + '='*20
            path = os.path.abspath(raw_input('Please Enter Path: '))
            c =+ 1

        assert choices, 'No movie log files found!' 

        return path


    def _get_movie_num(self,num):

        choices = glob.glob(self.path+'/movie.log.*')
        choices = [k[-3:] for k in choices]

        num = _num_to_ext(num)

        if num not in choices:
            _ =  'Select from the following possible moive numbers:'\
                 '\n{0}'.format(choices)
            num = int(raw_input(_))
 
        return _num_to_ext(num)


    def _get_movie_vars(self):
        #NOTE: The movie_vars are in an order, please do not switch around 
        #      unless you want incidious bugs

        #Check the moive header type
        if self.param['movie_header'] == '"movie2dC.h"':
            return ['rho',
                    'jx','jy','jz',
                    'bx','by','bz',
                    'ex','ey','ez',
                    'ne',
                    'jex','jey','jez',
                    'pexx','peyy','pezz','pexy','peyz','pexz',
                    'ni',
                    'pixx','piyy','pizz','pixy','piyz','pixz']

        elif self.param['movie_header'] == '"movie4b.h"':
            return ['rho',
                    'jx','jy','jz',
                    'bx','by','bz',
                    'ex','ey','ez',
                    'ne','jex','jey','jez',
                    'pexx','peyy','pezz','pexz','peyz','pexy',
                    'ni','jix','jiy','jiz',
                    'pixx','piyy','pizz' 'pixz','piyz','pixy']

        elif self.param['movie_header'] == '"movie2dD.h"':
            return ['rho',
                    'jx','jy','jz',
                    'bx','by','bz',
                    'ex','ey','ez',
                    'ne',
                    'jex','jey','jez',
                    'pexx','peyy','pezz','pexy','peyz','pexz',
                    'ni',
                    'jix','jiy','jiz',
                    'pixx','piyy','pizz','pixy','piyz','pixz']

        elif self.param['movie_header'] == '"movie_pic3.0.h"':
            return ['rho',
                    'jx','jy','jz',
                    'bx','by','bz',
                    'ex','ey','ez',
                    'ne',
                    'jex','jey','jez',
                    'pexx','peyy','pezz','pexy','peyz','pexz',
                    'ni',
                    'jix','jiy','jiz',
                    'pixx','piyy','pizz','pixy','piyz','pixz']
        else:
            err_msg = '='*80 + \
                      '\t This particular moive headder has not been coded!\n'\
                      '\tTalk to Colby to fit it, or fix it yourself.\n'\
                      '\tI dont care, Im a computer not a cop'\
                      '='*80

            print err_msg
            raise NotImplementedError()

################################################################################
############   Any thing below here has not been added yet  ####################
################################################################################

#c##Colby you need to figure out a way to make sure that
#c## the path is ok this should likly be done one level up.
#c#
#c#
#c#
#c#
#c##---------------------------------------------------------------------------------------------------------------
#c##   Method: load_movie
#c##   Args  : movie_num (to identify which of the posible several movie files to read from)
#c##         : movie_var (to identify which varible you want to read)
#c##         : movie_time (to identify which slice of time should be read)
#c##       This accepts the run name idetifies the coresponding
#c##       information in the run_list.dat file.
#c##---------------------------------------------------------------------------------------------------------------
#c#
#c#
#c#
#c#
#c#    def get_domain_arrays(self):
#c#        lx = self.param_dict['lx']
#c#        ly = self.param_dict['ly']
#c#        nx = self.param_dict['pex']*self.param_dict['nx']
#c#        ny = self.param_dict['pey']*self.param_dict['ny']
#c#        dx = lx/nx
#c#        dy = ly/ny
#c#        return (np.arange(dx/2.,lx,dx),np.arange(dy/2.,ly,dy))
