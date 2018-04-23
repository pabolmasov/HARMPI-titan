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

"""
this file was created specifically for remote processing of extensively large dumps outputs on a cluster. 
the challenge here is to get without plotting and ipython libraries; some of the procedures are cloned from harm_script
"""

def get_sorted_file_list(prefix="dump"):
    flist0 = np.sort(glob.glob( os.path.join("dumps/", "%s[0-9][0-9][0-9]_0000"%prefix) ) )
    flist1 = np.sort(glob.glob( os.path.join("dumps/", "%s[0-9][0-9][0-9][0-9]_0000"%prefix) ) )
    flist2 = np.sort(glob.glob( os.path.join("dumps/", "%s[0-9][0-9][0-9][0-9][0-9]_0000"%prefix) ) )
    flist0.sort()
    flist1.sort()
    flist2.sort()
    flist = np.concatenate((flist0,flist1,flist2))
    nlist=size(flist)
    for k in arange(nlist):
        flist[k]=flist[k][:-5] # cutting off the tails
    return flist

def faraday():
    #    global re.omegaf1, re.omegaf2
    if 're.omegaf1' in globals():
        del re.omegaf1
    if 're.omemaf2' in globals():
        del re.omegaf2
    omegaf1=old_div(re.fFdd(0,1),re.fFdd(1,3))
    omegaf2=old_div(re.fFdd(0,2),re.fFdd(2,3))
    return omegaf1, omegaf2

def Tcalcud():
    global Tud, TudEM, TudMA
    global mu, sigma
    global enth
    global unb, isunbound
    pg = (re.gam-1)*re.ug
    w=re.rho+re.ug+pg
    eta=w+re.bsq
    Tud = np.zeros((4,4,re.nx,re.ny,re.nz),dtype=np.float32,order='F')
    TudMA = np.zeros((4,4,re.nx,re.ny,re.nz),dtype=np.float32,order='F')
    TudEM = np.zeros((4,4,re.nx,re.ny,re.nz),dtype=np.float32,order='F')
    for kapa in np.arange(4):
        for nu in np.arange(4):
            if(kapa==nu): delta = 1
            else: delta = 0
            TudEM[kapa,nu] = re.bsq*re.uu[kapa]*re.ud[nu] + 0.5*re.bsq*delta - re.bu[kapa]*re.bd[nu]
            TudMA[kapa,nu] = w*re.uu[kapa]*re.ud[nu]+pg*delta
            #Tud[kapa,nu] = eta*uu[kapa]*ud[nu]+(pg+0.5*bsq)*delta-bu[kapa]*bd[nu]
            re.Tud[kapa,nu] = re.TudEM[kapa,nu] + re.TudMA[kapa,nu]
    mu = old_div(-Tud[1,0],(re.rho*re.uu[1]))
    sigma = old_div(re.TudEM[1,0],re.TudMA[1,0])
    enth=1+re.ug*re.gam/re.rho
    unb=enth*re.ud[0]
    isunbound=(-unb>1.0)

#############################################################
# tetrad covariant components... and all the grid parameters should be set as well
# on my laptop, the description is in ~/Arts/illum/frominside.tex
def tetrad_t(uumean, udmean):
    return old_div(-udmean[0],drdx[0,0]), old_div(-udmean[1],drdx[1,1]), 0.*udmean[0]/drdx[2,2], old_div(-udmean[3],drdx[3,3])

def tetrad_r(uumean, udmean):
    guu=re.guu ; drdx=re.drdx
    gg=guu[0,3]**2-guu[0,0]*guu[3,3] #1.+uumean[3]*udmean[3]
    hh=guu[3,3]+uumean[3]*uumean[3]
    s=sqrt(fabs(gg*hh*guu[1,1]))
    return guu[3,3]*uumean[1]/s/drdx[0,0], gg*udmean[0]/s/drdx[1,1], 0.*uumean[0]/drdx[2,2], -guu[0,3]*uumean[1]/s/drdx[3,3]

def tetrad_h(uumean, udmean):
    gdd=re.gdd ; drdx=re.drdx
    return 0.*uumean[0]/drdx[0,0], 0.*uumean[0]/drdx[1,1], old_div(sqrt(gdd[2,2]),drdx[2,2]), 0.*uumean[0]/drdx[3,3]

def tetrad_p(uumean, udmean):
    guu=re.guu ; drdx=re.drdx
    hh=sqrt(guu[3,3]+uumean[3]*uumean[3])
    return uumean[3]*udmean[0]/hh/drdx[0,0], uumean[3]*udmean[1]/hh/drdx[1,1], 0.*uumean[0]/drdx[2,2], (1.+uumean[3]*udmean[3])/hh/drdx[3,3]

# outputting all the basic information about the dump to a file
def dumpinfo(prefix):
    
    re.rg("gdump")
    re.rd(prefix)
    fout=open("dumps/"+prefix+"_dinfo.dat", "w")

    fout.write(str(re.nx)+" "+str(re.ny)+" "+str(re.nz)+"\n")
    print("nr x ntheta x nphi = "+str(re.nx)+" "+str(re.ny)+" "+str(re.nz))
    
    fout.write(str(re.a)+"\n")
    print("Kerr parameter a = "+str(re.a))

    fout.write(str(re.t)+"\n")
    print("time "+str(re.t))
    fout.close()

# calculated and writes out evolution of some global parameters; rref is the radius at which the mass and momentum flows are calculated
def glevol(nmax, rref):
#    global rho, uu
    re.rg("gdump")
    gdet=re.gdet ; gcov=re.gcov ; _dx2=re._dx2 ; _dx3=re._dx3
    # rref is reference radius
    run=unique(re.r)
    nr=run[where(run<rref)].argmax()

#    maccre=zeros(nmax, dtype=double)

    fmdot=open("merge_mevol"+str(rref)+".dat", "w")

    for k in arange(nmax):
        fname=re.dumpname(k)
        print("reading "+str(fname))
        re.rd(fname)
        print("rho is "+str(shape(re.rho)))
        #        Tcalcud()
        # accretion rate at rref
        rhom=(re.rho).mean(axis=2) ;  uum=(re.uu).mean(axis=3)
        rhom_south(re.rho*cos(re.h)).mean(axis=2)
        rhom_east(re.rho*sin(re.h)).mean(axis=2)
        # do we need to multiply this by drdx??
        #        uum[0]=(uu[0]*drdx[0,0]).mean(axis=3) ;    uum[1]=(uu[1]*drdx[1,1]).mean(axis=3)
        #        uum[2]=(uu[2]*drdx[2,2]).mean(axis=3) ;    uum[2]=(uu[2]*drdx[2,2]).mean(axis=3)
        gm=sqrt(old_div(gdet,gcov[1,1]))[:,0,:]
        maccre=-trapz((rhom*uum[1]*(uum[1]<0.)*gm)[nr,:], x=re.h[nr,:,0])*_dx2*_dx3
        maccre_south=-trapz((rhom_south*uum[1]*(uum[1]<0.)*gm)[nr,:], x=re.h[nr,:,0])*_dx2*_dx3
        maccre_east=-trapz((rhom_east*uum[1]*(uum[1]<0.)*gm)[nr,:], x=re.h[nr,:,0])*_dx2*_dx3
        mwind=-trapz((rhom*uum[1]*(uum[1]>0.)*gm)[nr,:], x=re.h[nr,:,0])*_dx2*_dx3
        laccre=-trapz((rho*fabs(old_div(ud[3],ud[0]))*uu[1]*(uu[1]<0.)*sqrt(old_div(gdet,gcov[1,1]))).mean(axis=2)[nr,:], x=re.h[nr,:,0])*_dx2*_dx3
        lwind=-trapz((rho*fabs(old_div(ud[3],ud[0]))*uu[1]*(uu[1]>0.)*sqrt(old_div(gdet,gcov[1,1]))).mean(axis=2)[nr,:], x=re.h[nr,:,0])*_dx2*_dx3
        # maccre=-((rho*uu[1]*sqrt(gdet/gcov[1,1]))[nr,:,:]).sum()*_dx2*_dx3
        fmdot.write(str(t)+" "+str(maccre)+" "+str(mwind)+" "+str(old_div(laccre,maccre))+" "+str(old_div(lwind,mwind))+" "+str(maccre_south)+" "+str(maccre_east)+"\n")
    
    fmdot.close()


def mint(rref):
    print("mint:")
    print("shape(rho) = "+str(shape(re.rho)))
    print("shape(uu1) = "+str(shape(re.uu[1])))
    print("shape gdet = "+str(shape(re.gdet)))
    hun=unique(re.h) ;     run=unique(re.r)
    gdet=re.gdet ; gcov=re.gcov ; _dx2=re._dx2 ; _dx3=re._dx3 ; drdx=re.drdx
    uu=re.uu ; ud=re.ud ; rho=re.rho ; _dx2=re._dx2 ; _dx3=re._dx3
    nr=run[where(run<rref)].argmax()
    indd=(rho*uu[1]*(uu[1]<0.)*sqrt(old_div(gdet,gcov[1,1])))[nr,:,:]
    maccre=-trapz((indd).mean(axis=-1), x=hun)*_dx2*_dx3
    indd=(rho*uu[1]*(uu[1]>0.)*sqrt(old_div(gdet,gcov[1,1])))[nr,:,:]
    mwind=trapz((indd).mean(axis=-1), x=hun)*_dx2*_dx3
    indd=(rho*uu[1]*(uu[1]<0.)*fabs(old_div(ud[3],drdx[3,3]))*sqrt(old_div(gdet,gcov[1,1])))[nr,:,:]
    laccre=-trapz((indd).mean(axis=-1), x=hun)*_dx2*_dx3
    indd=(rho*uu[1]*(uu[1]>0.)*fabs(old_div(ud[3],drdx[3,3]))*sqrt(old_div(gdet,gcov[1,1])))[nr,:,:]
    lwind=-trapz((indd).mean(axis=-1), x=hun)*_dx2*_dx3

    return maccre, mwind, old_div(laccre,maccre), old_div(lwind,mwind)

def readndump(n1, n2, rref=5.0):

#    run=unique(r)
#    nr=run[where(run<rref)].argmax()

    re.rg("gdump")
    nx=re.nx ; ny=re.ny ; nz=re.nz
    gdet=re.gdet ; gcov=re.gcov ; _dx2=re._dx2 ; _dx3=re._dx3 ; drdx=re.drdx
    r=re.r ; h=re.h ; phi=re.phi # importing coordinate mesh
    if (n2<n1):
        print("readndump: invalid file number range")
        exit()

    fmdot=open("merge_mevol"+str(rref)+".dat", "w")

    nframes=n2-n1+1
    n=n1+arange(nframes)

    for k in n:
        fname=re.dumpname(k)
        re.rd(fname)
        Tcalcud()
        p=(re.gam-1.)*re.ug
        magp=old_div(re.bsq,2.)
        rho=re.rho ; uu=re.uu ; ud=re.ud 
        if(k==n1):
            rhomean=rho
            rhosqmean=rho**2
            # velocity components:
            uu0=uu[0]*rho ; ud0=ud[0]*rho
            uur=uu[1]*rho ; udr=ud[1]*rho
            uuh=uu[2]*rho ; udh=ud[2]*rho
            uup=uu[3]*rho ; udp=ud[3]*rho
            puu0=uu[0]*p ; pud0=ud[0]*p
            puur=uu[1]*p ; pudr=ud[1]*p
            puuh=uu[2]*p ; pudh=ud[2]*p
            puup=uu[3]*p ; pudp=ud[3]*p
            mpuu0=uu[0]*magp ; mpud0=ud[0]*magp
            mpuur=uu[1]*magp ; mpudr=ud[1]*magp
            mpuuh=uu[2]*magp ; mpudh=ud[2]*magp
            mpuup=uu[3]*magp ; mpudp=ud[3]*magp

            tudem=re.TudEM ; tudma=re.TudMA
            pmean=p
	    # unorm=uaver
            magp_mean=magp
            aphi=re.psicalc()
        else:
            rhomean+=rho
            rhosqmean+=rho**2
            # velocity components:
            uu0+=uu[0]*rho ; ud0+=ud[0]*rho
            uur+=uu[1]*rho ; udr+=ud[1]*rho
            uuh+=uu[2]*rho ; udh+=ud[2]*rho
            uup+=uu[3]*rho ; udp+=ud[3]*rho
            puu0+=uu[0]*p ; pud0+=ud[0]*p
            puur+=uu[1]*p ; pudr+=ud[1]*p
            puuh+=uu[2]*p ; pudh+=ud[2]*p
            puup+=uu[3]*p ; pudp+=ud[3]*p
            mpuu0+=uu[0]*magp ; mpud0+=ud[0]*magp
            mpuur+=uu[1]*magp ; mpudr+=ud[1]*magp
            mpuuh+=uu[2]*magp ; mpudh+=ud[2]*magp
            mpuup+=uu[3]*magp ; mpudp+=ud[3]*magp

            pmean+=p
            magp_mean+=magp
            aphi+=re.psicalc()
            #	    unorm+=uaver
            tudem+=re.TudEM ; tudma+=re.TudMA
        maccre, mwind, laccre, lwind = mint(rref)
        fmdot.write(str(re.t)+" "+str(maccre)+" "+str(mwind)+" "+str(old_div(laccre,maccre))+" "+str(old_div(lwind,mwind))+"\n")
    fmdot.close()
    # velocity normalization:
    uu0*=old_div(drdx[0,0],rhomean) ; uur*=old_div(drdx[1,1],rhomean) ; uuh*=old_div(drdx[2,2],rhomean) ; uup*=old_div(drdx[3,3],rhomean)
    ud0/=drdx[0,0]*rhomean ; udr/=drdx[1,1]*rhomean ; udh/=drdx[2,2]*rhomean ; udp/=drdx[3,3]*rhomean
    puu0*=old_div(drdx[0,0],pmean) ; puur*=old_div(drdx[1,1],pmean) ; puuh*=old_div(drdx[2,2],pmean) ; puup*=old_div(drdx[3,3],pmean)
    pud0/=drdx[0,0]*pmean ; pudr/=drdx[1,1]*pmean ; pudh/=drdx[2,2]*pmean ; pudp/=drdx[3,3]*pmean
    mpuu0*=old_div(drdx[0,0],magp_mean) ; mpuur*=old_div(drdx[1,1],magp_mean) ; mpuuh*=old_div(drdx[2,2],magp_mean) ; mpuup*=old_div(drdx[3,3],magp_mean)
    mpud0/=drdx[0,0]*magp_mean ; mpudr/=drdx[1,1]*magp_mean ; mpudh/=drdx[2,2]*magp_mean ; mpudp/=drdx[3,3]*magp_mean

#    uu0/=unorm ; uur/=unorm ; uuh/=unorm ; uup/=unorm
#    ud0/=unorm ; udr/=unorm ; udh/=unorm ; udp/=unorm
    # averaging the density:
    rhomean=old_div(rhomean,double(nframes))
    rhodisp=old_div(rhosqmean,double(nframes))-rhomean**2
    pmean=old_div(pmean,double(nframes))
    magp_mean=old_div(magp_mean,double(nframes))
    tudem/=double(nframes) ; tudma/=double(nframes)
    aphi/=double(nframes)
    # physical stress-energy tensor components:
    for k in arange(4):
        for j in arange(4):
            tudem[k,j]*=drdx[k,k]/drdx[j,j]*sqrt(guu[k,k]*gdd[j,j])
            tudma[k,j]*=drdx[k,k]/drdx[j,j]*sqrt(guu[k,k]*gdd[j,j])

	#   ss=shape(rhomean) 
 	#   nx=ss[0]
  	#  ny=ss[1]
    # we need some 3D data
    rho3d=rhomean ;  p3d=pmean ; ur3d=uur ; uh3d=uuh ; up3d=uup
    fout=open('merge_rho3d.dat', 'w')
    for kx in arange(nx):
        for ky in arange(ny):
            for kz in arange(nz):
                fout.write(str(rho3d[kx,ky, kz])+'\n')
    fout.close()
    fout=open('merge_p3d.dat', 'w')
    for kx in arange(nx):
        for ky in arange(ny):
            for kz in arange(nz):
                fout.write(str(p3d[kx,ky, kz])+'\n')
    fout.close()
    fout=open('merge_u3d.dat', 'w')
    # uu on the mesh
    for kx in arange(nx):
        for ky in arange(ny):
            for kz in arange(nz):
            	fout.write(str(uu0[kx,ky,kz])+' '+str(uur[kx,ky,kz])+' '+str(uuh[kx,ky,kz])+' '+str(uup[kx,ky,kz])+'\n')
    fout.close()
    fout=open('merge_pu3d.dat', 'w')
    # uu on the mesh
    for kx in arange(nx):
        for ky in arange(ny):
            for kz in arange(nz):
                fout.write(str(puu0[kx,ky,kz])+' '+str(puur[kx,ky,kz])+' '+str(puuh[kx,ky,kz])+' '+str(puup[kx,ky,kz])+'\n')
    fout.close()
    fout=open('merge_mpu3d.dat', 'w')
    # uu on the mesh
    for kx in arange(nx):
        for ky in arange(ny):
            for kz in arange(nz):
                fout.write(str(mpuu0[kx,ky,kz])+' '+str(mpuur[kx,ky,kz])+' '+str(mpuuh[kx,ky,kz])+' '+str(mpuup[kx,ky,kz])+'\n')
    fout.close()

    # 2D or 3D? averaging over phi
    rhomean=rhomean.mean(axis=2)
    rhodisp=rhodisp.mean(axis=2)
    pmean=pmean.mean(axis=2) ; magp_mean=magp_mean.mean(axis=2)
    uu0=uu0.mean(axis=2) ; uur=uur.mean(axis=2) ; uuh=uuh.mean(axis=2) ; uup=uup.mean(axis=2)
    ud0=ud0.mean(axis=2) ; udr=udr.mean(axis=2) ; udh=udh.mean(axis=2) ; udp=udp.mean(axis=2)
    puu0=puu0.mean(axis=2) ; puur=puur.mean(axis=2) ; puuh=puuh.mean(axis=2) ; puup=puup.mean(axis=2)
    pud0=pud0.mean(axis=2) ; pudr=pudr.mean(axis=2) ; pudh=pudh.mean(axis=2) ; pudp=pudp.mean(axis=2)
    mpuu0=mpuu0.mean(axis=2) ; mpuur=mpuur.mean(axis=2) ; mpuuh=mpuuh.mean(axis=2) ; mpuup=mpuup.mean(axis=2)
    mpud0=mpud0.mean(axis=2) ; mpudr=mpudr.mean(axis=2) ; mpudh=mpudh.mean(axis=2) ; mpudp=mpudp.mean(axis=2)
    #  drp=drp.mean(axis=2)  ;  dhp=dhp.mean(axis=2)  ;  drh=drh.mean(axis=2)
    tudem=tudem.mean(axis=-1)
    tudma=tudma.mean(axis=-1)
    # aphi=aphi.mean(axis=-1)	
    # r mesh:
    fout=open('merge_r.dat', 'w')
    for kx in arange(nx):
        fout.write(str(r[kx,0,0])+'\n')
    fout.close()
    fout=open('merge_h.dat', 'w')
    # theta mesh
    for ky in arange(ny):
        fout.write(str(h[0,ky,0])+'\n')
    fout.close()
    fout=open('merge_phi.dat', 'w')
    # phi mesh
    for kz in arange(nz):
        fout.write(str(phi[0,0,kz])+'\n')
    fout.close()   
    fout=open('merge_rho.dat', 'w')
    # rho on the mesh
    for kx in arange(nx):
        for ky in arange(ny):
            fout.write(str(rhomean[kx,ky])+'\n')
    fout.close()
    # p on the mesh
    fout=open('merge_p.dat', 'w')
    for kx in arange(nx):
        for ky in arange(ny):
            fout.write(str(pmean[kx,ky])+'\n')
    fout.close()
    # magnetic pressure on the mesh
    fout=open('merge_mp.dat', 'w')
    for kx in arange(nx):
        for ky in arange(ny):
            fout.write(str(magp_mean[kx,ky])+'\n')
    fout.close()
    fout=open('merge_uu.dat', 'w')
    # uu on the mesh
    for kx in arange(nx):
        for ky in arange(ny):
            fout.write(str(uu0[kx,ky])+' '+str(uur[kx,ky])+' '+str(uuh[kx,ky])+' '+str(uup[kx,ky])+'\n')
    fout.close()
    fout=open('merge_ud.dat', 'w')
    # ud on the mesh
    for kx in arange(nx):
        for ky in arange(ny):
            fout.write(str(ud0[kx,ky])+' '+str(udr[kx,ky])+' '+str(udh[kx,ky])+' '+str(udp[kx,ky])+'\n')
    fout.close()
    fout=open('merge_puu.dat', 'w')
    # uu on the mesh
    for kx in arange(nx):
        for ky in arange(ny):
            fout.write(str(puu0[kx,ky])+' '+str(puur[kx,ky])+' '+str(puuh[kx,ky])+' '+str(puup[kx,ky])+'\n')
    fout.close()
    fout=open('merge_pud.dat', 'w')
    # ud on the mesh
    for kx in arange(nx):
        for ky in arange(ny):
            fout.write(str(pud0[kx,ky])+' '+str(pudr[kx,ky])+' '+str(pudh[kx,ky])+' '+str(pudp[kx,ky])+'\n')
    fout.close()
    fout=open('merge_mpuu.dat', 'w')
    # uu on the mesh
    for kx in arange(nx):
        for ky in arange(ny):
            fout.write(str(mpuu0[kx,ky])+' '+str(mpuur[kx,ky])+' '+str(mpuuh[kx,ky])+' '+str(mpuup[kx,ky])+'\n')
    fout.close()
    fout=open('merge_mpud.dat', 'w')
    # ud on the mesh
    for kx in arange(nx):
        for ky in arange(ny):
            fout.write(str(mpud0[kx,ky])+' '+str(mpudr[kx,ky])+' '+str(mpudh[kx,ky])+' '+str(mpudp[kx,ky])+'\n')
    fout.close()

    fout=open('merge_tudem.dat', 'w')
    # EM energy-stress tensor
    for kx in arange(nx):
        for ky in arange(ny):
            fout.write(str(tudem[0,0][kx,ky])+' '+str(tudem[0,1][kx,ky])+' '+str(tudem[0,2][kx,ky])+' '+str(tudem[0,3][kx,ky])+' ')
            fout.write(str(tudem[1,0][kx,ky])+' '+str(tudem[1,1][kx,ky])+' '+str(tudem[1,2][kx,ky])+' '+str(tudem[1,3][kx,ky])+' ')
            fout.write(str(tudem[2,0][kx,ky])+' '+str(tudem[2,1][kx,ky])+' '+str(tudem[2,2][kx,ky])+' '+str(tudem[2,3][kx,ky])+' ')
            fout.write(str(tudem[3,0][kx,ky])+' '+str(tudem[3,1][kx,ky])+' '+str(tudem[3,2][kx,ky])+' '+str(tudem[3,3][kx,ky])+'\n')
    fout.close()
    fout=open('merge_tudma.dat', 'w')
    # matter energy-stress tensor
    for kx in arange(nx):
        for ky in arange(ny):
            fout.write(str(tudma[0,0][kx,ky])+' '+str(tudma[0,1][kx,ky])+' '+str(tudma[0,2][kx,ky])+' '+str(tudma[0,3][kx,ky])+' ')
            fout.write(str(tudma[1,0][kx,ky])+' '+str(tudma[1,1][kx,ky])+' '+str(tudma[1,2][kx,ky])+' '+str(tudma[1,3][kx,ky])+' ')
            fout.write(str(tudma[2,0][kx,ky])+' '+str(tudma[2,1][kx,ky])+' '+str(tudma[2,2][kx,ky])+' '+str(tudma[2,3][kx,ky])+' ')
            fout.write(str(tudma[3,0][kx,ky])+' '+str(tudma[3,1][kx,ky])+' '+str(tudma[3,2][kx,ky])+' '+str(tudma[3,3][kx,ky])+'\n')
    fout.close()
    fout=open('merge_aphi.dat', 'w')
    # poloidal magnetic field lines	
    for kx in arange(nx):
        for ky in arange(ny):
            fout.write(str(aphi[kx,ky])+'\n')
    fout.close()

# reading dumps again to estimate the ellipsoid of turbulent velocities
def corvee(n1,n2):
    
    re.rg("gdump")
    nx=re.nx ; ny=re.ny ; nz=re.nz
    gdet=re.gdet ; gcov=re.gcov ; _dx2=re._dx2 ; _dx3=re._dx3 ; drdx=re.drdx
    r=re.r ; h=re.h ; phi=re.phi # importing coordinate mesh

    # velocities:
    uufile='merge_uu.dat'
    udfile='merge_ud.dat'
    fuu=open(uufile, 'r')
    fud=open(udfile, 'r')
    #    s=str.split(str.strip(fuu.readline()))
    # mean velocity field
    uumean=zeros([4,nx,ny,nz], dtype=double)
    udmean=zeros([4,nx,ny,nz], dtype=double)

    for kx in arange(nx):
        for ky in arange(ny):
            s=str.split(str.strip(fuu.readline()))
            uumean[0,kx,ky,:]=double(s[0])
            uumean[1,kx,ky,:]=double(s[1])
            uumean[2,kx,ky,:]=double(s[2])
            uumean[3,kx,ky,:]=double(s[3])
            s=str.split(str.strip(fud.readline()))
            udmean[0,kx,ky,:]=double(s[0])
            udmean[1,kx,ky,:]=double(s[1])
            udmean[2,kx,ky,:]=double(s[2])
            udmean[3,kx,ky,:]=double(s[3])
    fuu.close()
    fud.close()

    # tetrad components:
    t0=tetrad_t(uumean, udmean)
    tr=tetrad_r(uumean, udmean)
    th=tetrad_h(uumean, udmean)
    tp=tetrad_p(uumean, udmean)
#    print shape(tr)

    nframes=n2-n1+1
    n=n1+arange(nframes)

    for k in n:
        fname=re.dumpname(k)
        re.rd(fname)
        uu=re.uu ; ud=re.ud ; rho=re.rho
        if(k==n1):
            rhomean=rho
            # velocity components:
            uu0=uu[0]*drdx[0,0]-uumean[0] ; ud0=old_div(ud[0],drdx[0,0])-udmean[0]
            uur=uu[1]*drdx[1,1]-uumean[1] ; udr=old_div(ud[1],drdx[1,1])-udmean[1]
            uuh=uu[2]*drdx[2,2]-uumean[2] ; udh=old_div(ud[2],drdx[2,2])-udmean[2]
            uup=uu[3]*drdx[3,3]-uumean[3] ; udp=old_div(ud[3],drdx[3,3])-udmean[3]
            tuur=(uu0*tr[0]+uur*tr[1]+uuh*tr[2]+uup*tr[3])  # co-moving velocity components
            tuuh=(uu0*th[0]+uur*th[1]+uuh*th[2]+uup*th[3])
            tuup=(uu0*tp[0]+uur*tp[1]+uuh*tp[2]+uup*tp[3])
            drh=rho*tuur*tuuh ; drp=rho*tuur*tuup ; dhp=rho*tuuh*tuup
            drr=rho*tuur*tuur ; dpp=rho*tuup*tuup ; dhh=rho*tuuh*tuuh
#            print(shape(drh))
#            print(shape(rho))
#            print(shape(tuur))
        else:
            rhomean+=rho
            # velocity components:
            uu0=uu[0]-uumean[0] ; ud0=ud[0]-udmean[0]
            uur=uu[1]-uumean[1] ; udr=ud[1]-udmean[1]
            uuh=uu[2]-uumean[2] ; udh=ud[2]-udmean[2]
            uup=uu[3]-uumean[3] ; udp=ud[3]-udmean[3]
            tuur=(uu0*tr[0]+uur*tr[1]+uuh*tr[2]+uup*tr[3])
            tuuh=(uu0*th[0]+uur*th[1]+uuh*th[2]+uup*th[3])
            tuup=(uu0*tp[0]+uur*tp[1]+uuh*tp[2]+uup*tp[3])
            drh+=rho*tuur*tuuh ; drp+=rho*tuur*tuup ; dhp+=rho*tuuh*tuup
            drr+=rho*tuur*tuur ; dpp+=rho*tuup*tuup ; dhh+=rho*tuuh*tuuh

    drh/=rhomean ;    drp/=rhomean ;    dhp/=rhomean 
    drr/=rhomean ;    dpp/=rhomean ;    dhh/=rhomean 
  
    drh=drh.mean(axis=2) ;  drp=drp.mean(axis=2) ;  dhp=dhp.mean(axis=2)
    drr=drr.mean(axis=2) ;  dpp=dpp.mean(axis=2) ;  dhh=dhh.mean(axis=2)
 
    fout=open('merge_corv.dat', 'w')
    for kx in arange(nx):
        for ky in arange(ny):
            # RR HH PP RH HP RP
            fout.write(str(drr[kx,ky])+' '+str(dhh[kx,ky])+' '+str(dpp[kx,ky])+' '+str(drh[kx,ky])+' '+str(dhp[kx,ky])+' '+str(drp[kx,ky])+'\n')
    fout.close()

# making a horizontal slice
def fromabove(fname, htor=0.1, alifactor=1):

    # alifactor = alias factor, every alifactorth point in r and phi will be outputted, phi averaged over
#    rg("gdump")
    re.rd(fname)

    run=unique(re.r)
    hun=unique(re.h)
    phun=unique(re.ph)

    indisk=double(abs(cos(re.h))<htor)
    innorm=indisk.mean(axis=1)

    ugmean=old_div((re.ug*indisk).mean(axis=1),innorm)
    drdx=re.drdx ; uu=re.uu ; ud=re.ud ; rho=re.rho ; ug=re.ug # can we just import all the data?
    uu[0]*=drdx[0,0]  ; uu[1]*=drdx[1,1]  ; uu[2]*=drdx[2,2]  ; uu[3]*=drdx[3,3]  # to physical
    ud[0]/=drdx[0,0]  ; ud[1]/=drdx[1,1]  ; ud[2]/=drdx[2,2]  ; ud[3]/=drdx[3,3]  # to physical

    uumean0=old_div((uu[0]*indisk*rho).mean(axis=1),(rho*indisk).mean(axis=1))
    uumeanr=old_div((uu[1]*indisk*rho).mean(axis=1),(rho*indisk).mean(axis=1))
    uumeanh=old_div((uu[2]*indisk*rho).mean(axis=1),(rho*indisk).mean(axis=1))
    uumeanp=old_div((uu[3]*indisk*rho).mean(axis=1),(rho*indisk).mean(axis=1))
    uumean0_p=old_div((uu[0]*indisk*ug).mean(axis=1),(ug*indisk).mean(axis=1))
    uumeanr_p=old_div((uu[1]*indisk*ug).mean(axis=1),(ug*indisk).mean(axis=1))
    uumeanh_p=old_div((uu[2]*indisk*ug).mean(axis=1),(ug*indisk).mean(axis=1))
    uumeanp_p=old_div((uu[3]*indisk*ug).mean(axis=1),(ug*indisk).mean(axis=1))
    udmean0=old_div((ud[0]*indisk*rho).mean(axis=1),(rho*indisk).mean(axis=1))
    udmeanr=old_div((ud[1]*indisk*rho).mean(axis=1),(rho*indisk).mean(axis=1))
    udmeanh=old_div((ud[2]*indisk*rho).mean(axis=1),(rho*indisk).mean(axis=1))
    udmeanp=old_div((ud[3]*indisk*rho).mean(axis=1),(rho*indisk).mean(axis=1))
    udmean0_p=old_div((ud[0]*indisk*ug).mean(axis=1),(ug*indisk).mean(axis=1))
    udmeanr_p=old_div((ud[1]*indisk*ug).mean(axis=1),(ug*indisk).mean(axis=1))
    udmeanh_p=old_div((ud[2]*indisk*ug).mean(axis=1),(ug*indisk).mean(axis=1))
    udmeanp_p=old_div((ud[3]*indisk*ug).mean(axis=1),(ug*indisk).mean(axis=1))

    rhomean=old_div((rho*indisk).mean(axis=1),innorm)
    pmag=(re.bsq*indisk).mean(axis=1)/innorm/2.

    # r mesh:
    fout=open("dumps/"+fname+'_eq_r.dat', 'w')
    for kx in arange(re.nx):
        if(kx%alifactor==0):
            fout.write(str(run[kx])+'\n')
    fout.close()
    fout=open("dumps/"+fname+'_eq_phi.dat', 'w')
    # phi mesh
    for kz in arange(re.nz):
	#        if(kz%alifactor==0):
        fout.write(str(phun[kz])+'\n')
    fout.close()

    foutrho=open("dumps/"+fname+'_eq_rho.dat', 'w') # density
    foutp=open("dumps/"+fname+'_eq_p.dat', 'w') # gas pressure
    foutpm=open("dumps/"+fname+'_eq_pm.dat', 'w') # magnetic pressure
    foutu=open("dumps/"+fname+'_eq_uu.dat', 'w') # contravariant velocities
    foutd=open("dumps/"+fname+'_eq_ud.dat', 'w') # covariant velocities
    foutup=open("dumps/"+fname+'_eq_puu.dat', 'w') # pressure-averaged contravariant velocities
    foutdp=open("dumps/"+fname+'_eq_pud.dat', 'w') # pressure-averaged covariant velocities

    # rho on the mesh
    for kx in arange(re.nx):
        if(kx%alifactor==0):
            for kz in arange(re.nz):
                foutrho.write(str(rhomean[kx,kz])+'\n')
                foutp.write(str(ugmean[kx,kz]*(re.gam-1.))+'\n')
                foutpm.write(str(pmag[kx,kz])+'\n')
                foutu.write(str(uumean0[kx,kz])+' '+str(uumeanr[kx,kz])+' '+str(uumeanh[kx,kz])+' '+str(uumeanp[kx,kz])+'\n')
                foutd.write(str(udmean0[kx,kz])+' '+str(udmeanr[kx,kz])+' '+str(udmeanh[kx,kz])+' '+str(udmeanp[kx,kz])+'\n')
                foutup.write(str(uumean0_p[kx,kz])+' '+str(uumeanr_p[kx,kz])+' '+str(uumeanh_p[kx,kz])+' '+str(uumeanp_p[kx,kz])+'\n')
                foutdp.write(str(udmean0_p[kx,kz])+' '+str(udmeanr_p[kx,kz])+' '+str(udmeanh_p[kx,kz])+' '+str(udmeanp_p[kx,kz])+'\n')
    foutrho.close()
    foutp.close()
    foutpm.close()
    foutu.close()
    foutd.close()
    foutup.close()
    foutdp.close()
    os.system('tar -cf dumps/'+fname+'_eq.tar dumps/'+fname+'_dinfo.dat dumps/'+fname+'_eq_*.dat')

########################################################################################################################    
# reading one frame and making it lighter and 2D (by averaging over phi) and saving as an ascii
def framerip(fname, alifactor=1):
#    global rho, ug, uu, ud, gam, B, aphi
    # alifactor = alias factor, every alifactorth point in r and th will be outputted, phi averaged over
#    rg("gdump")

    re.rd("dumps/"+fname)
	#    print(rho)
    run=unique(re.r) ;    hun=unique(re.h)
    drdx=re.drdx ; uu=re.uu ; ud=re.ud ; rho=re.rho ; ug=re.ug ; B=re.B; bsq=re.bsq # can we just import all the data?
    # origin variables:
    orr=(re.origin_r*rho).mean(axis=2)/rho.mean(axis=2)
    orth=(re.origin_th*rho).mean(axis=2)/rho.mean(axis=2)
    orphi=(re.origin_phi*rho).mean(axis=2)/rho.mean(axis=2)
    
    ug=ug.mean(axis=2)
#    uu[0]*=drdx[0,0]  ; uu[1]*=drdx[1,1]  ; uu[2]*=drdx[2,2]  ; uu[3]*=drdx[3,3]  # to physical
    uu0=old_div((uu[0]*rho*drdx[0,0]).mean(axis=2),rho.mean(axis=2))
    uur=old_div((uu[1]*rho*drdx[1,1]).mean(axis=2),rho.mean(axis=2))
    uuh=old_div((uu[2]*rho*drdx[2,2]).mean(axis=2),rho.mean(axis=2))
    uup=old_div((uu[3]*rho*drdx[3,3]).mean(axis=2),rho.mean(axis=2))

#   uu[0]*=drdx[0,0]  ; uu[1]*=drdx[1,1]  ; uu[2]*=drdx[2,2]  ; uu[3]*=drdx[3,3]  # to physical
#    ud[0]/=drdx[0,0]  ; ud[1]/=drdx[1,1]  ; ud[2]/=drdx[2,2]  ; ud[3]/=drdx[3,3]  # to physical

#    ud=(ud*rho).mean(axis=3)/rho.mean(axis=2)
#    ud[0]/=drdx[0,0]  ; ud[1]/=drdx[1,1]  ; ud[2]/=drdx[2,2]  ; ud[3]/=drdx[3,3]  # to physical
    ud0=old_div((ud[0]*rho/drdx[0,0]).mean(axis=2),rho.mean(axis=2))
    udr=old_div((ud[1]*rho/drdx[1,1]).mean(axis=2),rho.mean(axis=2))
    udh=old_div((ud[2]*rho/drdx[2,2]).mean(axis=2),rho.mean(axis=2))
    udp=old_div((ud[3]*rho/drdx[3,3]).mean(axis=2),rho.mean(axis=2))

    rhom=rho.mean(axis=2)
    pmag=old_div((bsq).mean(axis=2),2.)
    aphi=re.psicalc()
    B[1]=B[1]*drdx[1,1] ; B[2]=B[2]*drdx[2,2] ; B[3]=B[3]*drdx[3,3]
    B=B.mean(axis=3)
#    aphi=aphi.mean(axis=2)

    # r mesh:
    fout=open("dumps/"+fname+'_r.dat', 'w')
    for kx in arange(re.nx):
        if(kx%alifactor==0):
            fout.write(str(run[kx])+'\n')
    fout.close()
    fout=open("dumps/"+fname+'_h.dat', 'w')
    # theta mesh
    for ky in arange(re.ny):
        if(ky%alifactor==0):
            fout.write(str(hun[ky])+'\n')
    fout.close()
    
    foutrho=open("dumps/"+fname+'_rho.dat', 'w')
    foutp=open("dumps/"+fname+'_p.dat', 'w')
    foutpm=open("dumps/"+fname+'_pm.dat', 'w')
    foutu=open("dumps/"+fname+'_uu.dat', 'w')
    foutd=open("dumps/"+fname+'_ud.dat', 'w')
    foutb=open("dumps/"+fname+'_b.dat', 'w')
    foutori=open("dumps/"+fname+'_ori.dat', 'w')
    # rho on the mesh
    for kx in arange(re.nx):
        if(kx%alifactor==0):
            for ky in arange(re.ny):
                if(ky%alifactor==0):
                    foutrho.write(str(rhom[kx,ky])+'\n')
                    foutp.write(str(ug[kx,ky]*(re.gam-1.))+'\n')
                    foutpm.write(str(old_div(pmag[kx,ky],2.))+'\n')
                    foutu.write(str(uu0[kx,ky])+' '+str(uur[kx,ky])+' '+str(uuh[kx,ky])+' '+str(uup[kx,ky])+'\n')
                    foutd.write(str(ud0[kx,ky])+' '+str(udr[kx,ky])+' '+str(udh[kx,ky])+' '+str(udp[kx,ky])+'\n')
                    foutb.write(str(B[1, kx,ky])+' '+str(B[2, kx,ky])+' '+str(B[3, kx,ky])+' '+str(aphi[kx,ky])+'\n')
                    foutori.write(str(orr[kx,ky])+' '+str(orth[kx,ky])+' '+str(orphi[kx,ky])+'\n')
    foutb.close()
    foutori.close()
    foutrho.close()
    foutp.close()
    foutpm.close()
    foutu.close()
    foutd.close()
    os.system('tar -cf dumps/'+fname+'.tar dumps/'+fname+'_*.dat')

# makes ascii files with degraded resolution from all the dumps
def defaultrun():

    dire=''
    rref=10.
 #   os.chdir(dire)
#    dumpinfo()
    re.rg(dire+"gdump")
    run=unique(re.r)
    nr=run[where(run<rref)].argmax() # where the radius is closest to rref (from inside)

    fout=open(dire+"dumps_mevol.dat", "w")

    flist=get_sorted_file_list()
    nlist=size(flist)
    print(str(nlist)+" files")
    print(flist)
    for k in arange(nlist):
        dumpinfo("../"+flist[k])
        framerip("../"+flist[k], alifactor=3)
        maccre, mwind, laccre, lwind = mint(rref)
        fromabove("../"+flist[k], alifactor=3)
        print("merger_remote defaultrun: reducing "+str(flist[k]))
        fout.write(str(re.t)+" "+str(maccre)+" "+str(mwind)+" "+str(laccre)+" "+str(lwind)+"\n")
    fout.close()

# produces time-averaged frames:
def corveerun(nfirst=None, nlast=None):
    re.rg("gdump")
    readndump(nfirst,nlast)
    corvee(nfirst, nlast)
    os.system('tar -cf mergereads.tar merge_*.dat')

defaultrun()
# corveerun()

