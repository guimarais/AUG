#!/usr/bin/env python
#import dd
#import dd_20180216 as dd
import dd
import kk_abock
import ww
#import ww_20180130 as ww #19/06/2019: "import ww" is not working anymore
from scipy.interpolate import interp1d
import numpy as np
import math as math
import sys

# Class to handle equilibrium
class EQUILLIBRIUM(object):
  
  def __init__(self,shotnumber):
    self.eq = kk_abock.kk()
    self.eq.Open(shotnumber, diag=equidiag)
  
  def get_kkBrzt(self,time,rs,zs):
    brzt = self.eq.get_B(time,rs,zs)

    bfield = np.sqrt(np.power(brzt['Br'],2)+np.power(brzt['Bt'],2)+np.power(brzt['Bz'],2))
    
    return rs,zs,bfield;

  def get_kkRZ2Rhop(self, time, r, z):
    rhop = self.eq.Rz_to_rhopol(time, r, z);

    return rhop

#################################################
def eq_usage(shotnumber, time, r, z):
    equi = EQUILLIBRIUM(shotnumber)
#  tlen = 10
#  npoints = 50
#  time = np.ones(tlen)
#  r = np.ones(tlen,npoints)
#  z = np.ones(tlen,npoints)
    rhop = np.empty(r.shape)
    for i in range(len(time)):
        rhop[i,:] = equi.get_kkRZ2Rhop(time[i], r[i,:], z[i,:])

    return rhop

############################################

shotnr = int(sys.argv[1])
equidiag = 'EQH'

##################################################
###############     Option Parsing
##################################################
##################################################
# Default Settings
#dens = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]

#Old campaigns with only Q
#dens = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25]
dens = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]

# Low Density Settings
#dens  =np.linspace(0.4, 2.0, 12)
#dens = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25]
#dems = [0.4, 0.6, 0.7, 0.8, 1.0, 1.15, 1.3, 1.4, 1.6, 1.7, 1.85, 2.0]
#dens = [0.5,0.7,0.9,1.05,1.1,1.15,1.25,1.3,1.35,1.5,1.7,1.9]
#dens = [0.402, 0.6, 0.839 ,1.0 ,1.2, 1.4, 1.608, 1.8, 2.0, 2.4, 2.8, 3.354]
#dens = [0.402, 0.6, 0.8, 0.973, 1.2, 1.4, 1.608, 1.8, 2.0, 2.4, 2.8, 3.354]
#dens = [0.402, 0.6, 0.839, 1.0, 1.2, 1.4, 1.608, 1.8, 2.0, 2.4, 2.8, 3.354]

nlayers = len(dens)

# Open the shotfile
rps = dd.shotfile('RPS', shotnr, experiment='AUGD')

edition = rps.edition

# Get LFS and HFS Radii and densitites/group delays

# Get the correct item and re-scale them:
nl = rps('neb_LFS')
nh = rps('neb_HFS')

rl = rps('RB_LFS')
rh = rps('RB_HFS')

#if shotnr > 33900:#Some prblem with writing RPS
#  tl = rl
#  th = rh#np.zeros_like(rh)
#else:

#2018/9/27: try again the group delay writing
tl = rps('GD_LFS')
th = rps('GD_HFS')

nl.data = nl.data * 1e-19
nh.data = nh.data * 1e-19

# Get the times for the REF
time = rps('TIME')

# Close the shotfile
rps.close()

poslfs    = []
poslfs_gd = []
poshfs    = []
poshfs_gd = []


for n, d in enumerate(dens):
    newposlfs = []
    newposhfs = []
    newposlfs_gd = []
    newposhfs_gd = []
    for i in range(time.size):
        #Interpolate Radii
        finterpl = interp1d(nl.data[i,:], rl.data[i,:], bounds_error=False) 
        finterph = interp1d(nh.data[i,:], rh.data[i,:], bounds_error=False) 
        newposlfs.append(finterpl(d))
        newposhfs.append(finterph(d))
        #Interpolate Group Delays
        ginterpl = interp1d(nl.data[i,:], tl.data[i,:], bounds_error=False) #
        ginterph = interp1d(nh.data[i,:], th.data[i,:], bounds_error=False) 
        newposlfs_gd.append(ginterpl(d))
        newposhfs_gd.append(ginterph(d))
    poslfs.append(newposlfs)
    poshfs.append(newposhfs)
    poslfs_gd.append(newposlfs_gd)
    poshfs_gd.append(newposhfs_gd)

#poslfs = np.zeros([len(dens), time.size]).T
#poslfs_gd = np.zeros([len(dens), time.size]).T
#poshfs = np.zeros([len(dens), time.size]).T
#poshfs_gd = np.zeros([len(dens), time.size]).T

#Convert to numpy arrays for easeness of handling, should be transposed
poslfs    = np.array(poslfs)
poshfs    = np.array(poshfs)
poslfs_gd = np.array(poslfs_gd)
poshfs_gd = np.array(poshfs_gd)

poslfs    = poslfs.T
poshfs    = poshfs.T
poslfs_gd = poslfs_gd.T
poshfs_gd = poshfs_gd.T

### Open the shotfile
diagnostic = 'RDL'
experiment = 'AUGD'
wwsf = ww.shotfile()
#if not wwsf.Open(experiment=experiment, diagnostic=diagnostic, shotnumber=shotnr, mode='new'):
if not wwsf.Open(experiment, diagnostic, shotnr):
    raise Exception("Error creating the new shotfile")
    exit()  

##Write the used densities
paramsetname = 'Aux'
signalname = 'ne'
nevalarray = np.array( dens ) * 1e19 
payload = np.array(nevalarray, dtype='float32')
wwsf.SetParameter(paramsetname,  signalname, payload)

##Write the RPS edition used for this shotfile
paramsetname = 'Aux'
signalname = 'RPS_ed'
edition = np.array( edition )
payload = np.array(edition, dtype='int')
wwsf.SetParameter(paramsetname, signalname, payload)

##Write the Magnetic Diagnostic
paramsetname = 'Aux'
signalname = 'Eq'
wwsf.SetParameter(paramsetname, signalname, equidiag)

###WRITE TIMEBASE
timename = 'TIME'
tdata = np.array(time, dtype='float32')
#wwsf.SetTimebase(timename, tdata)
wwsf.SetSignal(timename, tdata)

### LFS Radial positions and rho poloidal
siggrname = 'LFSR'
data = np.array(poslfs, dtype='float32')
wwsf.SetSignalGroup(siggrname, data)

zlfs = np.ones(poslfs.shape) * 0.14

rhol = eq_usage(shotnr, time, poslfs, zlfs)
siggrname = 'LFSRHO'
data = np.array(rhol, dtype='float32')
wwsf.SetSignalGroup(siggrname, data)

###
siggrname = 'HFSR'
data = np.array(poshfs, dtype='float32')
wwsf.SetSignalGroup(siggrname, data)

zhfs = np.ones(poshfs.shape) * 0.07

rhoh = eq_usage(shotnr, time, poshfs, zhfs)
siggrname = 'HFSRHO'
data = np.array(rhoh, dtype='float32')
wwsf.SetSignalGroup(siggrname, data)

###
siggrname = 'LFSGD'
data = np.array(poslfs_gd, dtype='float32')
wwsf.SetSignalGroup(siggrname, data)

###
siggrname = 'HFSGD'
data = np.array(poshfs_gd, dtype='float32')
wwsf.SetSignalGroup(siggrname, data)

###Close the shotfile
wwsf.Close()
