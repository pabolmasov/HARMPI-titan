from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import sys
import glob
import os

from numpy import sin, cos, tan, pi
# import numpy as np
from numpy import *

def readascmeshes(prefix):
    # reads stored r and theta meshes (in ASCII)    
    rfile=prefix+'_r.dat'
    hfile=prefix+'_h.dat'
    
    # radial mesh:
    fr=open(rfile, 'r')
    s=str.split(str.strip(fr.readline()))
    r=[]
    while(s):
        r.append(s[0])
        s=str.split(str.strip(fr.readline()))
    fr.close()
    r=asarray(r, dtype=double)
    nr=size(r)
    # polar angle mesh:
    fh=open(hfile, 'r')
    s=str.split(str.strip(fh.readline()))
    th=[]
    while(s):
        th.append(s[0])
        s=str.split(str.strip(fh.readline()))
    fh.close()
    th=asarray(th, dtype=double)
    nh=size(th)

    # 2d-grid (order??)
    h2,r2=meshgrid(th,r)
    print(shape(r2))
    print(nr, nh)
    return r2,h2

# reading a four-column text file:
def uread(infile, dims):

    lines=loadtxt('titania/merge_uu.dat')
    u0=reshape(asarray(lines[:,0], dtype=double), dims)
    u1=reshape(asarray(lines[:,1], dtype=double), dims)
    u2=reshape(asarray(lines[:,2], dtype=double), dims)
    u3=reshape(asarray(lines[:,3], dtype=double), dims)
    return u0, u1, u2, u3

def readasc_ori(orifile, newshape):
    fori=open(orifile, 'r')
    s=str.split(str.strip(fori.readline()))
    orr=[] ; orth=[] ; orphi=[]
    while(s):
        orr.append(s[0]) ;  orth.append(s[1]) ;  orphi.append(s[2])
        s=str.split(str.strip(fori.readline()))
    fori.close()
    orr=reshape(asarray(orr, dtype=double), newshape) 
    orth=reshape(asarray(orth, dtype=double), newshape) 
    orphi=reshape(asarray(orphi, dtype=double), newshape) 
    return orr, orth, orphi
    
def dinforead(prefix):

    # mesh data:
    dfile=prefix+'_dinfo.dat'
    fd=open(dfile, 'r')
    s=str.split(str.strip(fd.readline()))
    nr=int(s[0]) ; nh=int(s[1]) ; nphi=int(s[2])
    s=str.split(str.strip(fd.readline()))
    a=double(s[0])
    s=str.split(str.strip(fd.readline()))
    if(s):
        t=double(s[0])
    else:
        t=0.
    print("Kerr a = "+str(a))
    print(nr, nh)
    fd.close()
    return nr, nh, nphi, a, t

def readascone(fname, newshape):
    '''
    reading a stored quantity from ASCII file fname and reshaping it into newshape
    '''
    f=open(fname, 'r')
    s=str.split(str.strip(f.readline()))
    q=[]
    while(s):
        q.append(s[0])
        s=str.split(str.strip(f.readline()))
    f.close()
    q=asarray(q, dtype=double)
    q=reshape(q, newshape)
    return q
