from scipy.io.idl import readsav
init = False
if init:
    filename = '/glade/p/work/colbyh/2014.Winter.electron.heating.paper/dat_files/reconn621_lower_timeave.dat'
    filename = '/glade/p/work/colbyh/2014.Winter.electron.heating.paper/dat_files/tave100r305upperhalf-calcT.dat'
    filename = '/glade/p/work/colbyh/2014.Winter.electron.heating.paper/dat_files/reconn674_upper_timeave.dat'
    filename = '/glade/p/work/colbyh/2014.Winter.electron.heating.paper/dat_files/reconn691_upper_timeave.dat'
    filename = '/glade/p/work/colbyh/2014.Winter.electron.heating.paper/dat_files/tave100r305upperhalf-calcT.dat'
    CR = readsav(filename)
    filename = '/glade/p/work/colbyh/2014.Winter.electron.heating.paper/dat_files/reconn691_AVG_lower_timeave.dat'
    #CR1 = readsav(filename)
    filename = '/glade/p/work/colbyh/2014.Winter.electron.heating.paper/dat_files/reconn603_upper_timeave.dat'
    #CR2 = readsav(filename)
    extent = [CR['xx'][0],CR['xx'][-1], CR['yy'][0], CR['yy'][-1] ]
    #extent1 = [CR1['xx'][0],CR1['xx'][-1], CR1['yy'][0], CR1['yy'][-1] ]
    #extent2 = [CR2['xx'][0],CR2['xx'][-1], CR2['yy'][0], CR2['yy'][-1] ]

extent = [CR['xx'][0],CR['xx'][-1], CR['yy'][0], CR['yy'][-1] ]
# Pick location
r0 = array([250.,295.0])
dx0 = CR['xx'][2] - CR['xx'][1]

#r1 = array([126.,34.5])
#dx1 = CR1['xx'][2] - CR1['xx'][1]

#r2 = array([65., 76.8 + (34.5 - 25.6)])
#dx2 = CR2['xx'][2] - CR2['xx'][1]

aE0 = array([CR['exav'][int((r0[1]-CR['yy'][0])/dx0),int(r0[0]/dx0)], \
             CR['eyav'][int((r0[1]-CR['yy'][0])/dx0),int(r0[0]/dx0)], \
             CR['ezav'][int((r0[1]-CR['yy'][0])/dx0),int(r0[0]/dx0)] ])

aB0 = array([CR['bxav'][int((r0[1]-CR['yy'][0])/dx0),int(r0[0]/dx0)], \
             CR['byav'][int((r0[1]-CR['yy'][0])/dx0),int(r0[0]/dx0)], \
             CR['bzav'][int((r0[1]-CR['yy'][0])/dx0),int(r0[0]/dx0)] ])

v0 = cross(aE0,aB0)/dot(aB0,aB0)

#aE1 = array([CR1['exav'][int((r1[1]-CR1['yy'][0])/dx1),int(r1[0]/dx1)], \
#             CR1['eyav'][int((r1[1]-CR1['yy'][0])/dx1),int(r1[0]/dx1)], \
#             CR1['ezav'][int((r1[1]-CR1['yy'][0])/dx1),int(r1[0]/dx1)] ])

#aB1 = array([CR1['bxav'][int((r1[1]-CR1['yy'][0])/dx1),int(r1[0]/dx1)], \
#             CR1['byav'][int((r1[1]-CR1['yy'][0])/dx1),int(r1[0]/dx1)], \
#             CR1['bzav'][int((r1[1]-CR1['yy'][0])/dx1),int(r1[0]/dx1)] ])
#
#v1 = cross(aE1,aB1)/dot(aB1,aB1)
#
#aE2 = array([CR2['exav'][int((r2[1]-CR2['yy'][0])/dx2),int(r2[0]/dx2)], \
#             CR2['eyav'][int((r2[1]-CR2['yy'][0])/dx2),int(r2[0]/dx2)], \
#             CR2['ezav'][int((r2[1]-CR2['yy'][0])/dx2),int(r2[0]/dx2)] ])
#
#aB2 = array([CR2['bxav'][int((r2[1]-CR2['yy'][0])/dx2),int(r2[0]/dx2)], \
#             CR2['byav'][int((r2[1]-CR2['yy'][0])/dx2),int(r2[0]/dx2)], \
#             CR2['bzav'][int((r2[1]-CR2['yy'][0])/dx2),int(r2[0]/dx2)] ])
#
#v2 = cross(aE2,aB2)/dot(aB2,aB2)

TR = TPRun(mass=1.,charge=1.,npart=5,t0=5.,tend=150.,r0=r0,dr0=[2.*dx0,2.*dx0],v0=v0,dv0=[.05,.05,.05])
TR.move()

#CR = CR1
#TR_691 = TPRun(mass=1.,charge=1.,t0=5.,tend=50.,r0=r1,dr0=[2.*dx1,2.*dx1],npart=5,v0=v1,dv0=[.05,.05,.05])
#TR_691.move()

#CR = CR2
#TR_603 = TPRun(mass=1.,charge=1.,t0=5.,tend=50.,r0=r2,dr0=[2.*dx2,2.*dx2],npart=5,v0=v2,dv0=[.05,.05,.05])
#TR_603.move()

#TR = TR_691
#CR = CR1
#
#TR = TR_603
#CR = CR2
#
