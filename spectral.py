from __future__ import print_function
from __future__ import division
from builtins import str
from past.utils import old_div
from numpy import *
import numpy as np
import time
import os
from os.path import *
from os import fsync
# import harm_script as h
import glob
# import sys
# from reader import *
import reader as re
import scipy.interpolate as si
from matplotlib.cbook import flatten
# import matplotlib.pyplot as plt

# makes a Fourier transform of a patch of the simulation domain
def snapshot(prefix, x1=10., x2=12., y1=-5., y2=5.):
    re.rg("gdump")
    re.rd(prefix)

    # need to set up the coordinate network
    x=re.r*np.sin(re.h) ; y=re.r*cos(re.h)
    wbox=((((x-x1)*(x-x2))<0.)&(((y-y1)*(y-y2))<0.))
#    print("wbox shape = "+str(shape(wbox)))
    qty=re.rho
    # we need a way to switch between different variables
#    print("wbox shape = "+str(shape(re.rho[wbox])))
    run=unique(re.r[wbox]) ;   hun=unique(re.h[wbox])
    nx=np.size(run) ; ny=np.size(hun) ; nphi=re.nz

    rmean=re.r[wbox].mean() ; hmean=re.h[wbox].mean()
    
    xreg=(x2-x1)*np.arange(nx)/np.double(nx)+x1
    yreg=(y2-y1)*np.arange(ny)/np.double(ny)+y1
    dx=(x2-x1)/np.double(nx) ; dy=(y2-y1)/np.double(ny) ; dphi=2.*pi/np.double(nphi)

#    print("nx = "+str(nx)+", ny = "+str(ny))
    
    zint=np.zeros([nx,ny,nphi], dtype=double)    
    
    for kphi in arange(nphi):
        #        wbox=np.where((((re.r[:,:,kphi]-r1)*(re.r[:,:,kphi]-r2))<0.)&(((re.h[:,:,kphi]-th1)*(re.h[:,:,kphi]-th2))<0.))[0]
        # rr=(re.r[:,:,kphi])[wbox] ; hh=(re.h[:,:,kphi])[wbox] ;
        xx=x[wbox] ; yy=y[wbox] ; qq=np.reshape(qty[:,:,kphi],[re.nx,re.ny,1])
        qq=qq[wbox]
#        print("shape "+str(np.shape(xx))+", "+str(np.shape(yy)))
#        print("shape "+str(np.shape(qq)))
#        print("to "+str(np.shape(xreg))+", "+str(np.shape(yreg)))
        fint=si.SmoothBivariateSpline(xx,yy,qq,w=qq, kx=5, ky=5)
        #si.interp2d(yy, xx, qty)
        #        print(fint.ev(x1,y1))
        zint[:,:,kphi]=fint.__call__(xreg, yreg)
        # now we have a regularly-spaced 3Darray

    #    np.fft.rfftn
    if(nphi>1):
        zint_sp=np.fft.rfftn(zint)
    else:
        zint_sp=np.fft.rfft2(zint[:,:,0])
        zs=np.shape(zint_sp)
        zint_sp=np.reshape(zint_sp, [zs[0], zs[1], 1])
    zs=np.shape(zint_sp)
    kx=np.fft.fftfreq(zs[0], d=dx) ; ky=np.fft.fftfreq(zs[1], d=dy) ; kphi=np.fft.fftfreq(zs[2], d=dphi)/rmean/np.sin(hmean)
    return kx, ky, kphi, zint_sp
    
def powerstack(n1, n2, x1=10., x2=20., y1=-5., y2=5.):
    
    for k in np.arange(n2-n1)+n1:
        print("powerstack: reading "+re.dumpname(k))
        kx, ky, kphi, zint_sp = snapshot(re.dumpname(k), x1=x1, x2=x2, y1=y1,y2=-y1)
        if(k==n1):
            pds3d=np.abs(zint_sp)**2 ; dpds3d=np.abs(zint_sp)**4
            nkx=np.size(kx) ; nky=np.size(ky) ; nkphi=np.size(kphi)
        else:
            pds3d+=np.abs(zint_sp)**2 ; dpds3d+=np.abs(zint_sp)**4

    pds3d/=np.double(n2-n1)
    dpds3d=np.sqrt((dpds3d/np.double(n2-n1)-pds3d**2)/np.double(n2-n1-1)) # RMS

    # ascii output:
    fpds=open('pdsout.dat', 'w')
    for kkx in np.arange(nkx):
        for kky in np.arange(nky):
            for kkphi in np.arange(nkphi):
                fpds.write(str(kx[kkx])+" "+str(ky[kky])+" "+str(kphi[kkphi])+" "+str(pds3d[kkx,kky,kkphi])+" "+str(dpds3d[kkx,kky,kkphi])+"\n")
    fpds.close()
