from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import zip
from past.utils import old_div
import matplotlib
from matplotlib import rc
from matplotlib import axes
from numpy import *
from pylab import *
from scipy.integrate import *
from scipy.interpolate import *
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Ellipse

#Uncomment the following if you want to use LaTeX in figures 
rc('font',**{'family':'serif','serif':['Times']})
rc('mathtext',fontset='cm')
rc('mathtext',rm='stix')
rc('text', usetex=True)
# #add amsmath to the preamble
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amssymb,amsmath}"] 

import numpy as np
import matplotlib.pyplot as plt
import numpy.random
import time

import harm_script as h
import reader as rk

# last stable orbit
def arms(a):
    z1=1.+(1.-a*a)**(old_div(1.,3.))*((1.+a)**(old_div(1.,3.))+(1.-a)**(old_div(1.,3.)))
    z2=sqrt(3.*a*a+z1*z1)
#    if (a>0.):
#        r=3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2)) 
#    else: 
    r=3.+z2-sign(a)*sqrt((3.-z1)*(3.+z1+2.*z2))
    return r

# black hole horizon (Sasha's clone):
def bhole(rhor):

    ax = gca()
    el = Ellipse((0,0), 2*rhor, 2*rhor, facecolor='k', alpha=1)
    art=ax.add_artist(el)
    art.set_zorder(20)
    draw()

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
    fd.close()
    return nr, nh, nphi, a, t

def tedplotter(dire):
    nr, nh, nphi, a, t=dinforead(dire+'/merge')
    rfile=dire+'/merge_r.dat'
    fr=open(rfile, 'r')
    s=str.split(str.strip(fr.readline()))
    r=[]
    while(s):
        r.append(s[0])
        s=str.split(str.strip(fr.readline()))
    fr.close()
    r=asarray(r, dtype=double)
#    nr=size(r)
    # polar angle mesh:
    hfile='/home/pasha/harm/harmpi/'+dire+'/merge_h.dat'
    fh=open(hfile, 'r')
    s=str.split(str.strip(fh.readline()))
    th=[]
    while(s):
        th.append(s[0])
        s=str.split(str.strip(fh.readline()))
    fh.close()
    th=asarray(th, dtype=double)
#    nh=size(th)
    
    # 2d-grid (order??)
    h2,r2=meshgrid(th,r)
    print(shape(r2))
    print(nr, nh)

    # pressure:
    pfile=dire+'/merge_p.dat'
    fp=open(pfile, 'r')
    s=str.split(str.strip(fp.readline()))
    p=[]
    while(s):
        p.append(s[0])
        s=str.split(str.strip(fp.readline()))
    fp.close()
    p=asarray(p, dtype=double) ;    p=reshape(p, [nr, nh])

    # TudMA, TudEM
    trr=[] ; thh=[] ; tpp=[] ; trp=[] ; thp=[]
    tmafile=dire+'/merge_tudma.dat'
    ftma=open(tmafile, 'r')
    s=str.split(str.strip(ftma.readline()))
    rho=[]
    while(s):
        trr.append(s[5])
        thh.append(s[10])
        tpp.append(s[15])
        trp.append(s[7])
        thp.append(s[11])
        s=str.split(str.strip(ftma.readline()))
    ftma.close()
    trr=reshape(asarray(trr, dtype=double), [nr, nh]) 
    thh=reshape(asarray(thh, dtype=double), [nr, nh]) 
    tpp=reshape(asarray(tpp, dtype=double), [nr, nh]) 
    trp=reshape(asarray(trp, dtype=double), [nr, nh]) 
    thp=reshape(asarray(thp, dtype=double), [nr, nh]) 
    emtrr=[] ; emthh=[] ; emtpp=[] ; emtrp=[] ; emthp=[]
    temfile=dire+'/merge_tudem.dat'
    ftem=open(temfile, 'r')
    s=str.split(str.strip(ftem.readline()))
    rho=[]
    while(s):
        emtrr.append(s[5])
        emthh.append(s[10])
        emtpp.append(s[15])
        emtrp.append(s[7])
        emthp.append(s[11])
        s=str.split(str.strip(ftem.readline()))
    ftem.close()
    emtrr=reshape(asarray(emtrr, dtype=double), [nr, nh]) 
    emthh=reshape(asarray(emthh, dtype=double), [nr, nh]) 
    emtpp=reshape(asarray(emtpp, dtype=double), [nr, nh]) 
    emtrp=reshape(asarray(emtrp, dtype=double), [nr, nh]) 
    emthp=reshape(asarray(emthp, dtype=double), [nr, nh]) 
    
    alevs1=1e-3*0.5
    alevs2=1.0*0.5
    na=30
    alevs=(old_div(alevs2,alevs1))**(old_div(arange(na),double(na-1)))*alevs1
    alevs=around(alevs, 3)
    alevs[0]=0.
    alevs=unique(alevs)
    cmap = plt.get_cmap('jet')
    cmap.set_bad('white',1.)
    norm = BoundaryNorm(alevs, ncolors=cmap.N, clip=False)
    rmax=15.
    rhor = 1.+(1.-a**2)**0.5

    clf()
    fig=figure()
    subplot(121)
    contourf(r2*sin(h2), r2*cos(h2),fabs(old_div((trp+emtrp),p)), levels=alevs, norm=norm, cmap=cmap)
    colorbar()
    contour(r2*sin(h2), r2*cos(h2),(old_div((trp+emtrp),p)), colors='w', levels=[0.])
    contour(r2*sin(h2), r2*cos(h2),p, colors='w', linestyles='dotted')
    xlim(0., rmax)
    ylim(old_div(-rmax,2.), old_div(rmax,2.))
    xlabel(r'$\varpi$')
    ylabel(r'$z$')
    bhole(rhor)
    title(r'$\alpha_{r\varphi}$')
    subplot(122)
    contourf(r2*sin(h2), r2*cos(h2),fabs(old_div((thp+emthp),p)), levels=alevs, norm=norm, cmap=cmap)
    colorbar()
    contour(r2*sin(h2), r2*cos(h2),(old_div((thp+emthp),p)), colors='w', levels=[0.])
    contour(r2*sin(h2), r2*cos(h2),p, colors='w', linestyles='dotted')
    xlim(0., rmax)
    ylim(old_div(-rmax,2.), old_div(rmax,2.))
    xlabel(r'$\varpi$')
    ylabel(r'$z$')
    bhole(rhor)
    title(r'$\alpha_{z\varphi}$')
    fig.set_size_inches(15, 5)
    fig.tight_layout(pad=0., h_pad=-2.)
    savefig(dire+'/alphas.eps')
    close()

# reading all the averaged maps of the variables
def velread(prefix='merge_', nope=False, ifaphi=True):
    # prefix = 'titan2/merge_'
    rfile=prefix+'_r.dat'
    fr=open(rfile, 'r')
    s=str.split(str.strip(fr.readline()))
    r=[]
    while(s):
        r.append(s[0])
        s=str.split(str.strip(fr.readline()))
    fr.close()
    r=asarray(r, dtype=double)
    hfile=prefix+'_h.dat'
    fh=open(hfile, 'r')
    s=str.split(str.strip(fh.readline()))
    h=[]
    while(s):
        h.append(s[0])
        s=str.split(str.strip(fh.readline()))
    fh.close()
    h=asarray(h, dtype=double)
    h2,r2=meshgrid(h,r)
    nr,nh=shape(h2)

    # density:
    rhofile=prefix+'_rho.dat'
    frho=open(rhofile, 'r')
    s=str.split(str.strip(frho.readline()))
    rho=[]
    while(s):
        rho.append(s[0])
        s=str.split(str.strip(frho.readline()))
    frho.close()
    rho=asarray(rho, dtype=double)
    rho=reshape(rho, [nr, nh])

    # pressure:
    pfile=prefix+'_p.dat'
    fp=open(pfile, 'r')
    s=str.split(str.strip(fp.readline()))
    p=[]
    while(s):
        p.append(s[0])
        s=str.split(str.strip(fp.readline()))
    fp.close()
    p=asarray(p, dtype=double)
    p=reshape(p, [nr, nh])
    # magnetic pressure:
    pfile=prefix+'_mp.dat'
    fmp=open(pfile, 'r')
    s=str.split(str.strip(fmp.readline()))
    mp=[]
    while(s):
        mp.append(s[0])
        s=str.split(str.strip(fmp.readline()))
    fmp.close()
    mp=asarray(mp, dtype=double)
    mp=reshape(mp, [nr, nh])
    
    # velocities:
    uufile=prefix+'_uu.dat'
    udfile=prefix+'_ud.dat'
    puufile=prefix+'_puu.dat'
    pudfile=prefix+'_pud.dat'
    mpuufile=prefix+'_mpuu.dat'
    mpudfile=prefix+'_mpud.dat'
    uu0,uur, uuh, uup = rk.uread(uufile,[nr,nh]) # density-averaged velocity 
    ud0,udr, udh, udp = rk.uread(udfile,[nr,nh]) # density-averaged velocity 
    u0=sqrt(fabs(uu0*ud0))  ;   ur=sqrt(fabs(uur*udr))*sign(uur)  ;   uh=sqrt(fabs(uuh*udh))*sign(uuh)
    up=sqrt(fabs(uup*udp))*sign(uup)
    if(nope):
        pu0,pur, puh, pup = u0,ur, uh, up 
        mpu0,mpur, mpuh, mpup = u0,ur, uh, up 
    else:
        puu0,puur, puuh, puup = rk.uread(puufile,[nr,nh]) # pressure-averaged velocity 
        pud0,pudr, pudh, pudp = rk.uread(pudfile,[nr,nh]) # pressure-averaged velocity 
        mpuu0,mpuur, mpuuh, mpuup = rk.uread(mpuufile,[nr,nh]) # pressure-averaged velocity 
        mpud0,mpudr, mpudh, mpudp = rk.uread(mpudfile,[nr,nh]) # pressure-averaged velocity 
        pu0=sqrt(fabs(puu0*pud0))  ;   pur=sqrt(fabs(puur*pudr))*sign(puur)  ;   puh=sqrt(fabs(puuh*pudh))*sign(puuh) ; pup=sqrt(fabs(puup*pudp))*sign(puup)
        mpu0=sqrt(fabs(mpuu0*mpud0))  ;   mpur=sqrt(fabs(mpuur*mpudr))*sign(mpuur)  ;   mpuh=sqrt(fabs(mpuuh*mpudh))*sign(mpuuh) ; mpup=sqrt(fabs(mpuup*mpudp))*sign(mpuup)

    if(ifaphi):
        # vector potential A_phi:
        pfile=prefix+'_aphi.dat'
        faphi=open(pfile, 'r')
        s=str.split(str.strip(faphi.readline()))
        aphi=[]
        while(s):
            aphi.append(s[0])
            s=str.split(str.strip(faphi.readline()))
        faphi.close()
        aphi=asarray(aphi, dtype=double)
        aphi=reshape(aphi, [nr, nh])
    else:
        aphi=rho
            
    return r2, h2, rho, p, mp, u0, ur, uh, up, pu0, pur, puh, pup, mpu0, mpur, mpuh, mpup, aphi

# plots mean (phi- and t-averaged) maps of densities and velocities
def mplotter(dire,nope=False):

    dmatrix=True

    nr, nh, nphi, a, t=dinforead(dire+'/merge')
    r2, h2, rho, p, pm, u0, ur, uh, up, pu0, pur, puh, pup, mpu0, mpur, mpuh, mpup, aphi = velread(dire+'/merge')

    # 
    # velocity correlation matrix:
    if(dmatrix):
        dfile='/home/pasha/harm/harmpi/'+dire+'/merge_corv.dat'
        fd=open(dfile, 'r')
        #    s=str.split(str.strip(fd.readline()))
        dxy=zeros([3,3,nr,nh], dtype=double)
        #    vtrace1=zeros([nr,nh], dtype=double)
        for kx in arange(nr):
            for ky in arange(nh):
                s=str.split(str.strip(fd.readline()))
                dxy[0,0,kx,ky]=double(s[0])
                dxy[1,1,kx,ky]=double(s[1])
                dxy[2,2,kx,ky]=double(s[2])
                dxy[0,1,kx,ky]=double(s[3])
                dxy[1,0,kx,ky]=double(s[3])
                dxy[1,2,kx,ky]=double(s[4])
                dxy[2,1,kx,ky]=double(s[4])
                dxy[0,2,kx,ky]=double(s[5])
                dxy[2,0,kx,ky]=double(s[5])
                
        #
        vtrace=trace(dxy, axis1=0, axis2=1)
        #    print "vtrace = "+str(vtrace1.min())+" to "+str(vtrace1.max())
        #    print "vtrace = "+str(vtrace1.min())+" to "+str(vtrace1.max())
           # vertical slice:
        rrangemin=10. ; rrangemax=12.
        rrange=double((r2>rrangemin)*(r2<rrangemax))
        # averaging over radial velocity
        vtracemean=old_div((vtrace*rho*rrange).mean(axis=0),(rho*rrange).mean(axis=0))
        urmean=old_div((ur*rho*rrange).mean(axis=0),(rho*rrange).mean(axis=0))
        upmean=old_div((up*rho*rrange).mean(axis=0),(rho*rrange).mean(axis=0))
        uhmean=old_div((uh*rho*rrange).mean(axis=0),(rho*rrange).mean(axis=0))
        th=unique(h2)
        fig=figure()
        clf()
        plot(cos(th), sqrt(vtracemean), label='velocity RMS', color='b')
        plot(cos(th), urmean, label='radial velocity', color='r')
        plot(cos(th), -urmean, color='r', linestyle='dotted')
        plot(cos(th), upmean, label='rotation velocity', color='k')
        #        plot(cos(th), th*0.+rrangemin/(rrangemin**1.5+a), color='k', linestyle='dotted')
        #        plot(cos(th), th*0.+rrangemax/(rrangemax**1.5+a), color='k', linestyle='dotted')
        plot(cos(th), uhmean, label='latitudinal velocity', color='g')
        plot(cos(th), -uhmean, color='g', linestyle='dotted')
        yscale('log')
        legend(loc='best')
        xlabel(r'$\cos\theta$')
        ylabel('$v/c$')
        fig.set_size_inches(12, 6)
        savefig(dire+'/velcompare.eps')
        close()
        
    ono=20 # number of angular frequency levels
    rmin=old_div(h.Risco(a),2.)
    rhor = 1.+(1.-a**2)**0.5
        
    rmax=20.
    rlevs=(rmax/rmin*1.5)**(old_div(arange(ono),double(ono)))*rmin
    
    olevs=old_div(1.,(rlevs**1.5+a))
    olevs=olevs[::-1]
    olevs[ono-1]=olevs.max()*10.
    cmap = plt.get_cmap('jet')
    cmap.set_bad('white',1.)
    
#    grr=1./(1.-2./r+a**2/r**2)
        
    norm = BoundaryNorm(olevs, ncolors=cmap.N, clip=True)
    # density plot:
    clf()
    contourf(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), cmap=cmap, nlevels=30)
    contour(r2*sin(h2), r2*cos(h2), aphi, colors='k')
    xlim(0., rmax)
    ylim(old_div(-rmax,2.), old_div(rmax,2.))
    xlabel(r'$\varpi$')
    ylabel(r'$z$')
    bhole(rhor)
    savefig(dire+'/rho.eps')
    # beta magnetization plot:
    beta1=0.1 ; beta2=100. ; nbeta=30
    betalevs=log10((old_div(beta2,beta1))**(old_div(arange(nbeta),double(nbeta-1)))*beta1)
    clf()
    contourf(r2*sin(h2), r2*cos(h2), log10(old_div(p,pm)),levels=betalevs)
    colorbar()
    contour(r2*sin(h2), r2*cos(h2), log10(old_div(p,pm)),levels=[0.], colors='w', linewidths=2.)
    xlim(0., rmax)
    ylim(old_div(-rmax,2.), old_div(rmax,2.))
    xlabel(r'$\varpi$')
    ylabel(r'$z$')
    bhole(rhor)
    savefig(dire+'/beta.eps')
    # radial velocity
    clf()
    fig=figure()
    contourf(r2*sin(h2), r2*cos(h2), up,levels=olevs,norm=norm)
    colorbar()
    contour(r2*sin(h2), r2*cos(h2), old_div(1.,((r2*sin(h2))**1.5+a)), colors='k',levels=olevs,linewidths=1)
    contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), colors='w')
    xlim(0., 20.)
    ylim(-10., 10.)
    xlabel(r'$\varpi$')
    ylabel(r'$z$')
    bhole(rhor)
    fig.set_size_inches(8, 8)
    savefig(dire+'/omega.eps')
    vlevs=(arange(ono)/double(ono)*2.-1.)*0.01
    vlevs[0]=ur.min()*1.1
    vlevs[ono-1]=ur.max()*1.1
    norm = BoundaryNorm(vlevs, ncolors=cmap.N, clip=True)
    hdisk=0.25
    wdisk=double(fabs(cos(h2))<hdisk)
    wwind=double(fabs(cos(h2))>hdisk)
    urmean=old_div(((rho*ur)*wdisk).mean(axis=1),(rho*wdisk).mean(axis=1))
    urmeanp=old_div(((p*pur)*wdisk).mean(axis=1),(p*wdisk).mean(axis=1))

    # flow lines:
    nx=20 ; ny=21
    xmin=0. ; xmax=30.
    ymin=-10. ; ymax=10.

    cs=1. #p/rho/(4./3.)
    xflow=(xmax-xmin)*(arange(nx)+0.5)/double(nx)+xmin ;   yflow=(ymax-ymin)*(arange(ny)+0.5)/double(ny)+ymin
    x2,y2=meshgrid(xflow, yflow)
    #    vxfun=interp2d(r2*sin(h2), r2*cos(h2), ur*sin(h2)+uh*cos(h2),kind='linear')
    #    vyfun=interp2d(r2*sin(h2), r2*cos(h2), -uh*sin(h2)+ur*cos(h2),kind='linear')
    #    vx=vxfun(xflow, yflow)  ;  vy=vyfun(xflow, yflow)
    xgrid=(r2*sin(h2)).flatten()  ;  ygrid=(r2*cos(h2)).flatten()
    vxflow=(ur/cs*sin(h2)+uh/cs*cos(h2)).flatten()
    vyflow=(-uh/cs*sin(h2)+ur/cs*cos(h2)).flatten()
    pvxflow=(pur/cs*sin(h2)+puh/cs*cos(h2)).flatten()
    pvyflow=(-puh/cs*sin(h2)+pur/cs*cos(h2)).flatten()
#    vxflow=xgrid ; vyflow=ygrid
    vx=griddata(list(zip(xgrid, ygrid)),vxflow, (x2,y2),method='nearest')
    vy=griddata(list(zip(xgrid, ygrid)),vyflow, (x2,y2),method='nearest')
    pvx=griddata(list(zip(xgrid, ygrid)),pvxflow, (x2,y2),method='nearest')
    pvy=griddata(list(zip(xgrid, ygrid)),pvyflow, (x2,y2),method='nearest')

    vmin=1e-8 # sqrt((vx**2+vy**2)).min()*9.
    vmax=0.1 # sqrt((vx**2+vy**2)).max()*1.1
    vlevs=log10((old_div(vmax,vmin))**(old_div(arange(20),double(19)))*vmin)
    vlevs[0]=-30.
    norm = BoundaryNorm(vlevs, ncolors=cmap.N, clip=True)
    #    vmax=0.01
    #    vmin=-0.01
    #    vlevs=(vmax-vmin)*(arange(20)/double(19))+vmin
    #    norm = BoundaryNorm(vlevs, ncolors=cmap.N, clip=True)
    
    clf()
    fig=figure()
    contourf(r2*sin(h2), r2*cos(h2), log10(old_div(sqrt(ur**2+uh**2),cs)),levels=vlevs,norm=norm)
#    contourf(xflow, yflow, sqrt(vx**2+vy**2),levels=vlevs,norm=norm)
    colorbar()
    streamplot(xflow, yflow, pvx, pvy,color='k')
    streamplot(xflow, yflow, vx, vy,color='w')
    xlim(xmin, xmax)
    ylim(ymin, ymax)
    xlabel(r'$\varpi$')
    ylabel(r'$z$')
    bhole(rhor)
    fig.set_size_inches(15, 8)
    savefig(dire+'/stream.eps')
    close()
    # near eqplane:
    xscale=10.
    nx=7
    ny=5
    xflow=xscale*(arange(nx))/double(nx-1) ;   yflow=xscale*hdisk*((arange(ny))/double(ny-1)*2.-1.)
    x2,y2=meshgrid(xflow, yflow)    
    vx=griddata(list(zip(xgrid, ygrid)),vxflow, (x2,y2),method='nearest')
    vy=griddata(list(zip(xgrid, ygrid)),vyflow, (x2,y2),method='nearest')
    pvx=griddata(list(zip(xgrid, ygrid)),pvxflow, (x2,y2),method='nearest')
    pvy=griddata(list(zip(xgrid, ygrid)),pvyflow, (x2,y2),method='nearest')
    vratmin=0.6 # 0.2
    vratmax=1.1 # 1.
    nv=10
    vratlevs=(arange(nv+1))/double(nv)*(vratmax-vratmin)+vratmin
#    vratlevs[9]=1.3
    clf()
    fig=figure()
    # (sqrt(pur**2+puh**2))/(sqrt(ur**2+uh**2))
    contourf(r2*sin(h2), r2*cos(h2), old_div(pur,ur), levels=vratlevs, cmap='jet')
    colorbar()
    contour(r2*sin(h2), r2*cos(h2), old_div(pur,ur), levels=[1.], colors='w')
    plot([0.,xscale], [0.,0.], color='k', linestyle='dotted')
#    streamplot(xflow, yflow, pvx, pvy,color='k')
    streamplot(xflow, yflow, vx, vy,color='k')
    xlim(0.5, xscale)
    ylim(-xscale*hdisk, xscale*hdisk)
    xlabel(r'$\varpi$')
    ylabel(r'$z$')
    bhole(rhor)
    fig.set_size_inches(5*2+1, 5*2*hdisk+1.5)
    fig.tight_layout(pad=0.5)
    savefig(dire+'/streamband.eps')
    savefig(dire+'/streamband.jpg')
    close()
    vratmin=0.5 # 0.2
    vratmax=2.5 # 1.
    nv=10
    vratlevs=(arange(nv+1))/double(nv)*(vratmax-vratmin)+vratmin
    clf()
    fig=figure()
    # (sqrt(pur**2+puh**2))/(sqrt(ur**2+uh**2))
    contourf(r2*sin(h2), r2*cos(h2), old_div(mpur,pur), levels=vratlevs, cmap='jet')
    colorbar()
    contour(r2*sin(h2), r2*cos(h2), old_div(mpur,pur), levels=[1.], colors='w')
    plot([0.,xscale], [0.,0.], color='k', linestyle='dotted')
#    streamplot(xflow, yflow, pvx, pvy,color='k')
    streamplot(xflow, yflow, vx, vy,color='k')
    xlim(0.5, xscale)
    ylim(-xscale*hdisk, xscale*hdisk)
    xlabel(r'$\varpi$')
    ylabel(r'$z$')
    bhole(rhor)
    fig.set_size_inches(5*2+1, 5*2*hdisk+1.5)
    fig.tight_layout(pad=0.5)
    savefig(dire+'/streamband_mag.eps')
    savefig(dire+'/streamband_mag.jpg')
    close()

    # vertical slice:
    rrange=double((r2>5.)*(r2<10.))
    urhmean=old_div((ur*rho*rrange).mean(axis=0),(rho*rrange).mean(axis=0))
    uhhmean=old_div((uh*rho*rrange).mean(axis=0),(rho*rrange).mean(axis=0))
    urhmeanp=old_div((pur*p*rrange).mean(axis=0),(p*rrange).mean(axis=0))
    uhhmeanp=old_div((puh*p*rrange).mean(axis=0),(p*rrange).mean(axis=0))
    #    print shape(urhmean)
    #    print shape(h)
    th=unique(h2)
    clf()
    fig=figure()
    subplot(211)
    plot(cos(th), urhmean, color='k')
    plot(cos(th), urhmeanp, color='r')
#    plot(cos(th), uhhmean, color='k', linestyle='dotted')
#    plot(cos(th), uhhmeanp, color='r', linestyle='dotted')
    xlabel(r'$\cos\theta$')
    ylabel(r'$u^r$')
    ylim(-0.035,0.005)
    xlim(-hdisk, hdisk)
    subplot(212)
    plot(cos(th), old_div(urhmeanp,urhmean), color='k')
    xlabel(r'$\cos\theta$')
    ylabel(r'$\langle u^r\rangle_p / \langle u^r\rangle_\rho$')
    xlim(-hdisk, hdisk)
    ylim(0.,1.)
    fig.set_size_inches(8, 6)
    fig.tight_layout(pad=1.0,h_pad=0.5, w_pad=0.5)
    savefig('/home/pasha/harm/harmpi/'+dire+'/vverts.eps')
    close()
    clf()
    contourf(r2*sin(h2), r2*cos(h2), uh,levels=vlevs,norm=norm)
    colorbar()
#    contour(r2*sin(h2), r2*cos(h2), 1./((r2*sin(h2))**1.5+a), colors='k',levels=olevs,linewidths=1)
    contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), colors='w')
    xlim(0., 20.)
    ylim(-10., 10.)
    xlabel(r'$\varpi$')
    ylabel(r'$z$')
    bhole(rhor)
    savefig(dire+'/uh.eps')
    # turbulent velocity parameters:
    vmin=1e-8
    vmax=10.
    ono=100
    vlevs=(old_div(vmax,vmin))**(old_div(arange(ono),double(ono-1)))*vmin
    norm = BoundaryNorm(vlevs, ncolors=cmap.N, clip=True)
    if(dmatrix&False):
        clf()
        fig=figure()
        #    subplot(331)
        contourf(r2*sin(h2), r2*cos(h2), vtrace,levels=vlevs,norm=norm)
        colorbar()
        #    contour(r2*sin(h2), r2*cos(h2), 1./((r2*sin(h2))**1.5+a), colors='k',levels=olevs,linewidths=1)
        contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), colors='w',linestyles='dotted')
        #    contour(r2*sin(h2), r2*cos(h2), cos(h2), colors='y',linestyles='dashed',levels=[-hdisk, hdisk])
        xlim(0., 20.)
        ylim(-10., 10.)
        xlabel(r'$\varpi$')
        ylabel(r'$z$')
        bhole(rhor)
        savefig(dire+'/vturb.eps')
        close()
        vlevs=(2.*(old_div(arange(ono),double(ono-1)))-1.)*0.5
        vlevs[0]=-1.
        vlevs[ono-1]=1.
        norm = BoundaryNorm(vlevs, ncolors=cmap.N, clip=True)
        clf()
        fig=figure()
        subplot(331)
        contourf(r2*sin(h2), r2*cos(h2), old_div(dxy[0,0],vtrace),levels=vlevs,norm=norm)
        contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), colors='w',linestyles='dotted')
        title(r'$\Delta_{rr}$')
        xlim(0., 20.)
        ylim(-10., 10.)
        xlabel(r'$\varpi$')
        ylabel(r'$z$')
        bhole(rhor)
        subplot(332)
        contourf(r2*sin(h2), r2*cos(h2), old_div(dxy[0,1],vtrace),levels=vlevs,norm=norm)
        contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), colors='w',linestyles='dotted')
        title(r'$\Delta_{r\theta}$')
        xlim(0., 20.)
        ylim(-10., 10.)
        xlabel(r'$\varpi$')
        ylabel(r'$z$')
        bhole(rhor)
        subplot(333)
        contourf(r2*sin(h2), r2*cos(h2), old_div(dxy[0,2],vtrace),levels=vlevs,norm=norm)
        contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), colors='w',linestyles='dotted')
        title(r'$\Delta_{r\varphi}$')
        xlim(0., 20.)
        ylim(-10., 10.)
        xlabel(r'$\varpi$')
        ylabel(r'$z$')
        bhole(rhor)
        subplot(334)
        contourf(r2*sin(h2), r2*cos(h2), old_div(dxy[1,0],vtrace),levels=vlevs,norm=norm)
        contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), colors='w',linestyles='dotted')
        title(r'$\Delta_{\theta r}$')
        xlim(0., 20.)
        ylim(-10., 10.)
        xlabel(r'$\varpi$')
        ylabel(r'$z$')
        bhole(rhor)
        subplot(335)
        contourf(r2*sin(h2), r2*cos(h2), old_div(dxy[1,1],vtrace),levels=vlevs,norm=norm)
        contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), colors='w',linestyles='dotted')
        title(r'$\Delta_{\theta\theta}$')
        xlim(0., 20.)
        ylim(-10., 10.)
        xlabel(r'$\varpi$')
        ylabel(r'$z$')
        bhole(rhor)
        subplot(336)
        contourf(r2*sin(h2), r2*cos(h2), old_div(dxy[1,2],vtrace),levels=vlevs,norm=norm)
        contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), colors='w',linestyles='dotted')
        title(r'$\Delta_{\theta\varphi}$')
        xlim(0., 20.)
        ylim(-10., 10.)
        xlabel(r'$\varpi$')
        ylabel(r'$z$')
        bhole(rhor)
        subplot(337)
        contourf(r2*sin(h2), r2*cos(h2), old_div(dxy[2,0],vtrace),levels=vlevs,norm=norm)
        contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), colors='w',linestyles='dotted')
        title(r'$\Delta_{\varphi r}$')
        xlim(0., 20.)
        ylim(-10., 10.)
        xlabel(r'$\varpi$')
        ylabel(r'$z$')
        bhole(rhor)
        subplot(338)
        contourf(r2*sin(h2), r2*cos(h2), old_div(dxy[2,1],vtrace),levels=vlevs,norm=norm)
        contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), colors='w',linestyles='dotted')
        title(r'$\Delta_{\varphi\theta}$')
        xlim(0., 20.)
        ylim(-10., 10.)
        xlabel(r'$\varpi$')
        ylabel(r'$z$')
        bhole(rhor)
        subplot(339)
        contourf(r2*sin(h2), r2*cos(h2), old_div(dxy[2,2],vtrace),levels=vlevs,norm=norm)
        contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3), colors='w',linestyles='dotted')
        title(r'$\Delta_{\varphi\varphi}$')
        xlim(0., 20.)
        ylim(-10., 10.)
        xlabel(r'$\varpi$')
        ylabel(r'$z$')
        bhole(rhor)
        fig.set_size_inches(12, 12)
        fig.tight_layout(pad=1.0,h_pad=0.5, w_pad=0.5)
        savefig(dire+'/dmatrix.eps')
        close()

    if(dmatrix&False):
        # the tetrad has reasonable physical sense only if u^h << u^r,phi
        drrdisk=old_div((dxy[0,0]*wdisk*rho).mean(axis=1),(wdisk*rho).mean(axis=1))
        dhhdisk=old_div((dxy[1,1]*wdisk*rho).mean(axis=1),(wdisk*rho).mean(axis=1))
        dppdisk=old_div((dxy[2,2]*wdisk*rho).mean(axis=1),(wdisk*rho).mean(axis=1))
        drpdisk=old_div((dxy[0,2]*wdisk*rho).mean(axis=1),(wdisk*rho).mean(axis=1))
        drhdisk=old_div((dxy[0,1]*wdisk*rho).mean(axis=1),(wdisk*rho).mean(axis=1))
        dhpdisk=old_div((dxy[1,2]*wdisk*rho).mean(axis=1),(wdisk*rho).mean(axis=1))
        dhpdiskplus=old_div((dxy[1,2]*wdisk*cos(h2)*rho).mean(axis=1),(wdisk*rho).mean(axis=1))
        #    dhpdiskplus=(dxy[1,2]*wdisk*cos(h2)*rho).mean(axis=1)/(wdisk*rho).mean(axis=1)
        drhdiskplus=old_div((dxy[1,0]*wdisk*cos(h2)*rho).mean(axis=1),(wdisk*rho).mean(axis=1))
        #    drhdiskplus=(dxy[1,0]*wdisk*cos(h2)*rho).mean(axis=1)/(wdisk*rho).mean(axis=1)

        dtot=drrdisk+dhhdisk+dppdisk

        clf()
        plot(r, old_div(drrdisk,dtot), color='k')
        plot(r, old_div(dhhdisk,dtot), color='g')
        plot(r, old_div(dppdisk,dtot), color='r')
        plot(r, old_div(drpdisk,dtot), color='r', linestyle='dotted')
        plot(r, old_div(drhdisk,dtot), color='g', linestyle='dotted')
        plot(r, old_div(dhpdisk,dtot), color='orange', linestyle='dotted')
        plot(r, old_div(dhpdiskplus,dtot), color='orange', linestyle='dashed')
        plot(r, old_div(drhdiskplus,dtot), color='g', linestyle='dashed')
        plot(r*0.+h.Risco(a), arange(nr)/double(nr-1)*2.-1., color='k', linestyle='dotted')
        xlabel(r'$r$')
        ylabel(r'$\Delta_{ik} / \Delta_{\rm tot}$')
        xscale('log')
        xlim(1,20)
        #    ylim(-1e-2,1e-2)
        savefig(dire+'/dmatrix_rslice.eps')

        drrvert=old_div((dxy[0,0]*rrange*rho).mean(axis=0),(rrange*rho).mean(axis=0))
        dhhvert=old_div((dxy[1,1]*rrange*rho).mean(axis=0),(rrange*rho).mean(axis=0))
        dppvert=old_div((dxy[2,2]*rrange*rho).mean(axis=0),(rrange*rho).mean(axis=0))
        drpvert=old_div((dxy[0,2]*rrange*rho).mean(axis=0),(rrange*rho).mean(axis=0))
        drhvert=old_div((dxy[0,1]*rrange*rho).mean(axis=0),(rrange*rho).mean(axis=0))
        dhpvert=old_div((dxy[1,2]*rrange*rho).mean(axis=0),(rrange*rho).mean(axis=0))
        dvertot=drrvert+dhhvert+dppvert
        clf()
        plot(cos(th), old_div(drrvert,dvertot), color='k')
        plot(cos(th), old_div(dhhvert,dvertot), color='g')
        plot(cos(th), old_div(dppvert,dvertot), color='r')
        plot(cos(th), old_div(drpvert,dvertot), color='r', linestyle='dotted')
        plot(cos(th), old_div(drhvert,dvertot), color='k', linestyle='dotted')
        plot(cos(th), old_div(dhpvert,dvertot), color='orange', linestyle='dotted')
        xlabel(r'$\cos \theta$')
        ylabel(r'$\Delta_{ik} / \Delta_{\rm tot}$')
        savefig(dire+'/dmatrix_thslice.eps')

# reading R\Theta file in ascii ; plotting the results of framerip
def ascframe(prefix='dumps/dump000', xmax=20.):

    rfile=prefix+'_r.dat'
    hfile=prefix+'_h.dat'
    rhofile=prefix+'_rho.dat'
    pfile=prefix+'_p.dat'
    pmfile=prefix+'_pm.dat'
    uufile=prefix+'_uu.dat'
    udfile=prefix+'_ud.dat'
    bfile=prefix+'_b.dat'
    orifile=prefix+'_ori.dat'

    nr, nh, nphi, a, t=dinforead(prefix)

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

    # density:
    frho=open(rhofile, 'r')
    s=str.split(str.strip(frho.readline()))
    rho=[]
    while(s):
        rho.append(s[0])
        s=str.split(str.strip(frho.readline()))
    frho.close()
    rho=asarray(rho, dtype=double)
    rho=reshape(rho, [nr, nh])

    # pressure:
    fp=open(pfile, 'r')  ;  fpm=open(pmfile, 'r')
    s=str.split(str.strip(fp.readline()))
    sm=str.split(str.strip(fpm.readline()))
    p=[] ; pm=[]
    while(s):
        p.append(s[0]) ;    pm.append(sm[0])
        s=str.split(str.strip(fp.readline()))
        sm=str.split(str.strip(fpm.readline()))
    fp.close() ;   fpm.close()
    p=asarray(p, dtype=double)  ;  pm=asarray(pm, dtype=double)
    p=reshape(p, [nr, nh])   ;  pm=reshape(pm, [nr, nh]) 

    # velocity field:
    fuu=open(uufile, 'r')
    s=str.split(str.strip(fuu.readline()))
    ur=[] ; uh=[] ; omega=[]
    while(s):
        ur.append(s[1])
        uh.append(s[2])
        omega.append(old_div(double(s[3]),double(s[0])))
        s=str.split(str.strip(fuu.readline()))
    fuu.close()
    ur=asarray(ur, dtype=double) ;   uh=asarray(uh, dtype=double)
    ur=reshape(ur, [nr, nh]) ; uh=reshape(uh, [nr, nh])
    omega=asarray(omega, dtype=double)
    omega=reshape(omega, [nr, nh])

    # origin variables:
    fori=open(orifile, 'r')
    s=str.split(str.strip(fori.readline()))
    orr=[] ; orth=[] ; orphi=[]
    while(s):
        orr.append(s[0]) ;  orth.append(s[1]) ;  orphi.append(s[2])
        s=str.split(str.strip(fori.readline()))
    fori.close()
    orr=reshape(asarray(orr, dtype=double), [nr, nh]) 
    orth=reshape(asarray(orth, dtype=double), [nr, nh]) 
    orphi=reshape(asarray(orphi, dtype=double), [nr, nh]) 
    # magnetic field (the last component is A_\phi)
    fb=open(bfile, 'r')
    s=str.split(str.strip(fb.readline()))
    aphi=[] 
    while(s):
        aphi.append(s[3])
        s=str.split(str.strip(fb.readline()))
    fb.close()
    aphi=reshape(asarray(aphi, dtype=double), [nr, nh]) 
    print("size(aphi) = "+str(shape(aphi)))
    
    rhor=1.+sqrt(1.-a**2)

    cmap = plt.get_cmap('jet')
    ono=30
    lmin=-5.
    lmax=1.
    lrholevs=(lmax-lmin)*arange(ono)/double(ono)+lmin
    norm = BoundaryNorm(lrholevs, ncolors=cmap.N, clip=True)
    x=r2*sin(h2) ; y=r2*cos(h2)
    clf()
    fig=figure()
    contourf(x, y, log10(rho+1e-3),levels=lrholevs,norm=norm,cmap=cmap)
    contour(x, y, aphi, colors='k')
    xlim(0.,xmax)
    ylim(-xmax/4.,xmax/2.)
    bhole(rhor)
    # need to put time in dinfo!
    title('t='+str(t)+' ('+prefix+')')
    savefig(prefix+'_rho.eps')
    savefig(prefix+'_rho.png')
    close()
    nxx=10
    rlevs=xmax*np.arange(nxx)/np.double(nxx) ; thlevs=np.pi*np.arange(nxx)/np.double(nxx)
    clf()
    contourf(x, y, rho,cmap=cmap)
    contour(x, y, r2, colors='w', levels=rlevs)
    contour(x, y, h2, colors='w', levels=thlevs)
    contour(x, y, orr, colors='k', levels=rlevs)
    contour(x, y, orth, colors='k', levels=thlevs)
    plt.xlim(0., xmax) ; plt.ylim(-xmax/4., xmax/2.)
    plt.savefig(prefix+"_ori.png")
    close()
    lmin=0.
    lmax=3.
    lbetalevs=(lmax-lmin)*arange(ono)/double(ono)+lmin
    norm = BoundaryNorm(lbetalevs, ncolors=cmap.N, clip=True)
    clf()
    fig=figure()
    contourf(r2*sin(h2), r2*cos(h2), log10(old_div(p,pm)),levels=lbetalevs,norm=norm,cmap=cmap)
    contour(r2*sin(h2), r2*cos(h2), aphi, colors='k')
    xlim(0.,xmax)
    ylim(-xmax/2.,xmax/2.)
    bhole(rhor)
    # need to put time in dinfo!
    title('t='+str(t)+' ('+prefix+')')
    savefig(prefix+'_beta.eps')
    savefig(prefix+'_beta.png')
    close()
    vlevs=(arange(ono)/double(ono)*2.-1.)*0.1
    wv=where((r2<10.)&(r2>5.)&(fabs(cos(h2))<0.25))
#    vlevs[0]=ur[wv].min()*1.5
#    vlevs[ono-1]=ur[wv].max()*1.5
    norm = BoundaryNorm(vlevs, ncolors=cmap.N, clip=True)
    clf()
    fig=figure()
    contourf(r2*sin(h2), r2*cos(h2), ur,levels=vlevs,norm=norm)
    colorbar()
    contour(r2*sin(h2), r2*cos(h2), ur,levels=[0.], color='w', linestyles='dotted')
    contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3),levels=lrholevs, colors='w')
    xlim(0.,xmax)
    ylim(-xmax/2.,xmax/2.)
    bhole(rhor)
    title('t='+str(t))
    savefig(prefix+'_ur.eps')
    savefig(prefix+'_ur.png')
    close()
    clf()
    fig=figure()
    contourf(r2*sin(h2), r2*cos(h2), omega,levels=vlevs,norm=norm)
    colorbar()
    contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3),levels=lrholevs, colors='w')
    xlim(0.,xmax)
    ylim(-xmax/2.,xmax/2.)
    bhole(rhor)
    title('t='+str(t))
    savefig(prefix+'_o.eps')
    savefig(prefix+'_o.png')
    close()

    nx=20 ; ny=20
    xmin=0. ; ymin=-xmax/4. ; ymax=xmax/2.
    xflow=(xmax-xmin)*(arange(nx)+0.5)/double(nx)+xmin ;   yflow=(ymax-ymin)*(arange(ny)+0.5)/double(ny)+ymin
    x2,y2=meshgrid(xflow, yflow)
    xgrid=(r2*sin(h2)).flatten()  ;  ygrid=(r2*cos(h2)).flatten()
    vxflow=(ur*sin(h2)+uh*cos(h2)).flatten()
    vyflow=(-uh*sin(h2)+ur*cos(h2)).flatten()
    vx=griddata(list(zip(xgrid, ygrid)),vxflow, (x2,y2),method='nearest')
    vy=griddata(list(zip(xgrid, ygrid)),vyflow, (x2,y2),method='nearest')

    vmin=1e-8 # sqrt((vx**2+vy**2)).min()*9.
    vmax=1.0 # sqrt((vx**2+vy**2)).max()*1.1
    vlevs=log10((old_div(vmax,vmin))**(old_div(arange(20),double(19)))*vmin)
    vlevs[0]=-30.
    norm = BoundaryNorm(vlevs, ncolors=cmap.N, clip=True)
    clf()
    fig=figure()
    contourf(r2*sin(h2), r2*cos(h2), log10(sqrt(ur**2+uh**2)),levels=vlevs,norm=norm)
    colorbar()
    title('t='+str(t))
    contour(r2*sin(h2), r2*cos(h2), log10(rho+1e-3),levels=lrholevs, colors='w')
    streamplot(xflow, yflow, vx, vy,color='k')
    xlim(xmin, xmax)
    ylim(ymin, ymax)
    xlabel(r'$\varpi$')
    ylabel(r'$z$')
    bhole(rhor)
    plot([arms(a), arms(a)], [-1., 1.], color='w') #, linestyle='dotted')
    fig.set_size_inches(10, 6)
    savefig(prefix+'_stream.eps')
    savefig(prefix+'_stream.png')
    close()

# reading R\Phi file in ascii 
def eqframe(prefix):

    rfile=prefix+'_eq_r.dat'
    phifile=prefix+'_eq_phi.dat'
    rhofile=prefix+'_eq_rho.dat'
    pfile=prefix+'_eq_p.dat'
    pmfile=prefix+'_eq_pm.dat'
    uufile=prefix+'_eq_uu.dat'
    udfile=prefix+'_eq_ud.dat'

    nr, nh, nphi, a, t=dinforead(prefix)

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
    # azimuthal angle mesh:
    fh=open(phifile, 'r')
    s=str.split(str.strip(fh.readline()))
    phi=[]
    while(s):
        phi.append(s[0])
        s=str.split(str.strip(fh.readline()))
    fh.close()
    phi=asarray(phi, dtype=double)
#    nh=size(phi)

    # 2d-grid (order??)
    h2,r2=meshgrid(phi,r)
    print(shape(r2))
    print(nr, nphi)

    # density:
    frho=open(rhofile, 'r')
    s=str.split(str.strip(frho.readline()))
    rho=[]
    while(s):
        rho.append(s[0])
        s=str.split(str.strip(frho.readline()))
    frho.close()
    rho=asarray(rho, dtype=double)
#    print shape(rho)
    rho=reshape(rho, [nr, nphi])
    # pressure(s):
    fp=open(pfile, 'r') ;   fpm=open(pmfile, 'r')
    s=str.split(str.strip(fp.readline()))  ;  sm=str.split(str.strip(fpm.readline()))
    p=[] ; pm=[]
    while(s):
        p.append(s[0])     ;   pm.append(sm[0])
        s=str.split(str.strip(fp.readline())) ;  sm=str.split(str.strip(fpm.readline()))
    fp.close() ; fpm.close()
    p=asarray(p, dtype=double) ; pm=asarray(pm, dtype=double)
    #    print shape(rho)
    p=reshape(p, [nr, nphi])  ;  pm=reshape(pm, [nr, nphi])

    # velocity field:
    fuu=open(uufile, 'r')
    s=str.split(str.strip(fuu.readline()))
    ur=[] ; omega=[]
    while(s):
        ur.append(s[1])
        omega.append(old_div(double(s[3]),double(s[0])))
        s=str.split(str.strip(fuu.readline()))
    fuu.close()
    ur=asarray(ur, dtype=double)
    ur=reshape(ur, [nr, nphi])
    omega=asarray(omega, dtype=double)
    omega=reshape(omega, [nr, nphi])

    rhor=1.+sqrt(1.-a**2)

    # density variations:
    rhomean=rho.mean(axis=1)
    drho=zeros([nr, nphi], dtype=double)
    for k in arange(nr):
        drho[k,:]=old_div(rho[k,:],rhomean[k])-1.
    drholevs=levels=(old_div(arange(40),double(20.))-0.5)*5.
    clf()
    fig=figure()
    contourf(r2*sin(phi), r2*cos(phi), drho,levels=drholevs)
    xlim(0.,20.)
    ylim(-10.,10.)
    bhole(rhor)
    title('deviations from mean density profile, t='+str(t))
    #    axis('equal')
    savefig(prefix+'_eq_rho.eps')
    savefig(prefix+'_eq_rho.png')
    close()
    clf()
    contourf(r2*sin(phi), r2*cos(phi), ur*rho,nlevels=30)
    colorbar()
    xlim(0.,20.)
    ylim(-10.,10.)
    bhole(rhor)
    title(r'$\rho u^r$, t='+str(t))
    #    axis('equal')
    savefig(prefix+'_eq_ur.eps')
    savefig(prefix+'_eq_ur.png')
    close()
    # beta (magnetization) plot
    beta1=0.1 ; beta2=1000. ; nbeta=20
    betalevs=log10((old_div(beta2,beta1))**(old_div(arange(nbeta),double(nbeta-1)))*beta1)
    clf()
    contourf(r2*sin(phi), r2*cos(phi), log10(old_div(p,pm)),levels=betalevs)
    colorbar()
    contour(r2*sin(phi), r2*cos(phi), drho, colors='k',levels=drholevs, lineswidth=1)
    xlim(0.,15.)
    ylim(-5.,10.)
    bhole(rhor)
    title(r'$\lg\beta$, t='+str(t))
    #    axis('equal')
    savefig(prefix+'_eq_beta.eps')
    savefig(prefix+'_eq_beta.png')
    

def dumpmovie():
    dire='dumps/'
    n1=0
    n2=340
    for k in n1+arange(n2-n1+1):
        prefix=dire+rk.dumpname(k)
        print(prefix)
#        eqframe(prefix)
        ascframe(prefix)
        
    # then: ffmpeg -framerate 15 -pattern_type glob -i 'titan2/dump???_rho.png' -b 4096k titanic2.mp4
    # ffmpeg -framerate 15 -pattern_type glob -i 'titan2/dump???_rho.png' -b 4096k titanic2_rho_1.mp4
    # ffmpeg -framerate 15 -pattern_type glob -i 'titan2/dump????_rho.png' -b 4096k titanic2_rho_2.mp4
    # ffmpeg -framerate 15 -pattern_type glob -i 'titan2/dump???_eq_beta.png' -b 4096k titanic2_beta1.mp4
    # ffmpeg -framerate 15 -pattern_type glob -i 'titan2/dump????_eq_beta.png' -b 4096k titanic2_beta2.mp4

# making a plot of a PDS produced by powerstack
def readsp():

    qq = np.loadtxt('pdsout.dat')
    kr=qq[:,0] ; kh=qq[:,1] ; kphi=qq[:,2] ;  pds3d=qq[:,3] ; dpds3d=qq[:,4] 
    krun=unique(kr) ; khun=unique(kh) ; kpun=unique(kphi)
    nr=np.size(krun) ; nh=np.size(khun) ; nphi=np.size(kpun)
    pds3d=np.reshape(pds3d, [nr,nh,nphi])
    dpds3d=np.reshape(dpds3d, [nr,nh,nphi])
    kr=np.reshape(kr, [nr,nh,nphi])
    kh=np.reshape(kh, [nr,nh,nphi])
    kphi=np.reshape(kphi, [nr,nh,nphi])
    krun=kr[:,0,0] ; khun=kh[0,:,0] ; kpun=kphi[0,0,:]
    print(shape(pds3d))
    # rk.uread('pdsout.dat')
    # plotting
    pds2d=pds3d.mean(axis=-1) # average over phi
    dpds2d=dpds3d.mean(axis=-1) # average over phi
    print(shape(pds2d))
    pdsr=pds2d.mean(axis=1) ; dpdsr=dpds2d.mean(axis=1)
    print(shape(pdsr))
    pdsh=pds2d.mean(axis=0) ; dpdsh=dpds2d.mean(axis=0)
    print(shape(pdsh))
    if(nphi>1):
        pdsphi=pds3d.mean(axis=0).mean(axis=0) ; dpdsphi=dpds3d.mean(axis=0).mean(axis=0)
        
    plt.clf()
    plt.plot(krun, pdsr, '.k')
    plt.errorbar(krun, pdsr, yerr=dpdsr, color='k', fmt='.')
    plt.plot(khun, pdsh, '.g')
    plt.errorbar(khun, pdsh, yerr=dpdsh, color='g', fmt='.')
    if(nphi>1):
        plt.plot(kpun, pdsphi, '.r')
        plt.errorbar(kpun, pdsphi, yerr=dpdsphi, color='r', fmt='.')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('pdss.png')
    
    
    
    
