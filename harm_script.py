import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl

from matplotlib import rc
from matplotlib.patches import Ellipse
from scipy.interpolate import interp1d
from matplotlib.gridspec import GridSpec
from matplotlib import cm,ticker
from numpy import sin, cos, tan, pi
#matplotlib.use('Agg') #so that it does ok with graphics in batch mode

#choose Computer Modern Roman fonts by default
mpl.rcParams['font.serif'] = 'cmr10'
mpl.rcParams['font.sans-serif'] = 'cmr10'

#font = { 'size'   : 20}
#rc('font', **font)
rc('xtick', labelsize=20) 
rc('ytick', labelsize=20) 
#rc('xlabel', **font) 
#rc('ylabel', **font) 

#Uncomment the following if you want to use LaTeX in figures 
rc('font',**{'family':'serif','serif':['Times']})
rc('mathtext',fontset='cm')
rc('mathtext',rm='stix')
rc('text', usetex=True)
# #add amsmath to the preamble
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amssymb,amsmath}"] 

# legend = {'fontsize': 20}
# rc('legend',**legend)
axes = {'labelsize': 20}
rc('axes', **axes)
rc('mathtext',fontset='cm')
#use this, but at the expense of slowdown of rendering
#rc('text', usetex=True)
# #add amsmath to the preamble
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amssymb,amsmath}"] 
import pdb

import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from numpy import ma
import matplotlib.colors as colors
#use_math_text = True

from reader import *

def mathify_axes_ticks(ax,fontsize=20,xticks=None,yticks=None):
    if xticks is None:
        xticks = ax.get_xticks()
    if yticks is None:
        yticks = ax.get_yticks()
    if ax.get_xscale() != 'log': ax.set_xticklabels([(r'$%g$' % lab) for lab in xticks])
    if ax.get_yscale() != 'log': ax.set_yticklabels([(r'$%g$' % lab) for lab in yticks])
    if fontsize is not None:
        if ax.get_xscale() != 'log':
            for label in ax.get_xticklabels():
                label.set_fontsize(fontsize)
        if ax.get_yscale() != 'log':
            for label in ax.get_yticklabels():
                label.set_fontsize(fontsize)


def convert_to_single_file(startn=0,endn=-1,ln=10,whichi=0,whichn=1,**kwargs):
    which = kwargs.pop("which","convert_file")
    rg("gdump")
    flist1 = np.sort(glob.glob( os.path.join("dumps/", "dump[0-9][0-9][0-9]_0000") ) )
    flist2 = np.sort(glob.glob( os.path.join("dumps/", "dump[0-9][0-9][0-9][0-9]_0000") ) )
    flist1.sort()
    flist2.sort()
    flist = np.concatenate((flist1,flist2))
    firsttime = 1
    for fldname in flist:
        #find the index of the file
        fldindex = np.int(fldname.split("_")[0].split("p")[-1])
        if fldindex < startn:
            continue
        if endn>=0 and fldindex >= endn:
            break
        if fldindex % whichn != whichi:
            #do every whichn'th snapshot starting with whichi'th snapshot
            continue
        #print( "Reading " + fldname + " ..." )
        fname = "dump%03d" % fldindex
        if os.path.isfile( fname ):
            print("File %s exists, skipping..." % fname)
            continue
        if not os.path.isfile( fname ):
            rd(fname)
        

def ellk(a,r):
    ekval = ek(a,r)
    lkval = lk(a,r)
    return(lkval/ekval)

def ek(a,r):
    #-u_t, I suspect
    ek = (r**2-2*r+a*r**0.5)/(r*(r**2-3*r+2*a*r**0.5)**0.5)
    return(ek)

def lk(a,r):
    udphi = r**0.5*(r**2-2*a*r**0.5+a**2)/(r*(r**2-3*r+2*a*r**0.5)**0.5)
    return( udphi )

def Risco(ain):
    eps = np.finfo(np.float64).eps
    a = np.minimum(ain,1.)
    Z1 = 1 + (1. - a**2)**(1./3.) * ((1. + a)**(1./3.) + (1. - a)**(1./3.))
    Z2 = (3*a**2 + Z1**2)**(1./2.)
    risco = 3 + Z2 - np.sign(a)* ( (3 - Z1)*(3 + Z1 + 2*Z2) )**(1./2.)
    return(risco)

def Ebind(r,a):
    #1+u_t, I suspect
    Eb = 1 - (r**2-2*r+a*r**0.5)/(r*(r**2-3*r+2*a*r**0.5)**0.5)
    return( Eb )

def etaNT(a):
    return( Ebindisco(a) )

def Ebindisco(a):
    eps = np.finfo(np.float64).eps
    a0 = 0.99999 #1.-1e8*eps
    if a > a0: 
        a = a0
        Eb = Ebind( Risco(a), a )
        return((a-a0)/(1.-a0)*(1.-3.**(-0.5)) + (1.-a)/(1.-a0)*Eb)
    Eb = Ebind( Risco(a), a)
    #Eb = (1.-3.**(-0.5))*a**2
    return( Eb )

def mkmov_simple(starti=0,endi=400):
    for i in xrange(starti,endi+1):
        rd("dump%03d" % i);
        aphi=psicalc()
        if i == starti: amax = aphi.max()
        cs, cb = plco(np.log10(rho),levels=np.linspace(-8,0,100),isfilled=1,k=0,xy=1,xmax=10,ymax=5,dobh=1,cb=1,extend="both",pretty=1)
        ax = plt.gca()
        ax.set_xlabel(r"$R\ [r_g]$",fontsize=20,labelpad=-5)
        ax.set_ylabel(r"$z\ [r_g]$",fontsize=20,labelpad=-5)
        cb.ax.set_xlabel(r"$\log\rho$",fontsize=20,ha="left")
        plc(aphi,levels=np.linspace(-amax,amax,10)[1:-1],colors="white",linewidths=2,xy=-1)
        print i;
        plt.title("t=%.4g"%np.round(t)); 
        plt.draw();
        plt.savefig("frame%03d.png"%i)

def convert_wrapper(**kwargs):
    if len(sys.argv[2:])==2 and sys.argv[2].isdigit() and sys.argv[3].isdigit():
        whichi = int(sys.argv[2])
        whichn = int(sys.argv[3])
    else:
        print( "Usage: %s %s <whichi> <whichn>" % (sys.argv[0], sys.argv[1]) )
        return
    convert_to_single_file(whichi = whichi, whichn = whichn, **kwargs)

def mkmov_wrapper(**kwargs):
    if len(sys.argv[2:])==2 and sys.argv[2].isdigit() and sys.argv[3].isdigit():
        whichi = int(sys.argv[2])
        whichn = int(sys.argv[3])
    else:
        print( "Usage: %s %s <whichi> <whichn>" % (sys.argv[0], sys.argv[1]) )
        return
    mkmov(whichi = whichi, whichn = whichn, **kwargs)

def mkmov(startn=0,endn=-1,ln=10,whichi=0,whichn=1,**kwargs):
    which = kwargs.pop("which","mkfrm8panel")
    dosavefig = kwargs.pop("dosavefig",1)
    print("Doing %s movie" % which)
    rg("gdump")
    #compute the total magnetic flux at t = 0
    rd("dump000")
    aphi=psicalc()
    aphimax = aphi.max()
    #construct file list
    flist1 = np.sort(glob.glob( os.path.join("dumps/", "dump[0-9][0-9][0-9]") ) )
    flist2 = np.sort(glob.glob( os.path.join("dumps/", "dump[0-9][0-9][0-9][0-9]") ) )
    flist1.sort()
    flist2.sort()
    flist = np.concatenate((flist1,flist2))
    if len(flist) == 0:
        flist1 = np.sort(glob.glob( os.path.join("dumps/", "dump[0-9][0-9][0-9]_0000") ) )
        flist2 = np.sort(glob.glob( os.path.join("dumps/", "dump[0-9][0-9][0-9][0-9]_0000") ) )
        flist1.sort()
        flist2.sort()
        flist = np.concatenate((flist1,flist2))
    firsttime = 1
    dpi = 135
    for fldname in flist:
        #find the index of the file
        fldindex = np.int(fldname.split("_")[0].split("p")[-1])
        if fldindex < startn:
            continue
        if endn>=0 and fldindex >= endn:
            break
        if fldindex % whichn != whichi:
            #do every whichn'th snapshot starting with whichi'th snapshot
            continue
        if dosavefig:
            fname = "%s%04d.png" % (which,fldindex)
            if os.path.isfile( fname ):
                print("File %s exists, skipping..." % fname)
                continue
        #print( "Reading " + fldname + " ..." )
        rd("dump%03d" % fldindex);
        if which == "mkfrmsimple":
            if firsttime:
                firsttime = 0
                fig = plt.figure(figsize=(12,8))
                plt.clf()
            mkfrmsimple(fig=fig,aphimax = aphimax)
        else:
            print("Unknown movie type: %s" % which)
            return
        print fldindex
        plt.draw()
        if dosavefig:
            plt.savefig(fname,dpi = dpi)

            
#############
def mkfrmsimple(fig=None,aphimax=None,lnx=100,lny=100,vmin=-10,vmax=1,fntsize=20,asp=1.):
    if fig is None: fig = plt.gcf();
    aphi = psicalc() #vpot[3].mean(-1)
    if aphimax is None: aphimax = aphi.max()
    #ax.set_aspect(asp)
    res,cb=plco(lrho,xy=1,xmax=lnx,ymax=lny,symmx=1,
                isfilled=1,cb=1,pretty=1,
                levels=np.linspace(vmin,vmax,100),
                extend="both",cbxla=r"$\ \ \ \ \ \ \ \ \log_{10}\rho$")
    plt.xlabel(r"$x\ [r_g]$",fontsize=fntsize)
    plt.ylabel(r"$z\ [r_g]$",fontsize=fntsize,labelpad=-30)
    ax = plt.gca()
    #cmap = cm.jet
    #label = r"$\log_{10}\rho$"
    #cx1,cb1 = mkvertcolorbar(ax,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="right",
    #            label=label,ticks=tcks,fntsize=fntsize,cmap=cmap,extend="both")
    plc(aphi/aphimax,symmx=1,xy=-1,levels=np.linspace(0.,1.,20)[1:],colors="black",linewidths=1.)
    plt.title(r"$t=%g$" % int(t+0.5), fontsize=fntsize)
    plt.xlim(-lnx,lnx)
    plt.ylim(-lny,lny)
    mathify_axes_ticks(ax)
    
def mkvertcolorbar(ax,fig,vmin=0,vmax=1,label=None,ylabel=None,ticks=None,fntsize=20,cmap=mpl.cm.jet,gap=0.03,width=0.02,extend="neither",loc="right"):
    box = ax.get_position()
    #pdb.set_trace()
    # cpos = [box.x0,box.y0+box.height+0.05,box.width,0.03]
    locs = loc.split()
    loc0 = locs[0]
    if len(locs)>1:
        loc1 = locs[1]
    else:
        loc1 = None
    if loc0 == "left":
        cpos = box.x0-gap-width,box.y0,width,box.height
    elif loc0 == "right":
        cpos = box.x0+box.width+gap,box.y0,width,box.height
    elif loc0 == "top":
        if loc1 == "right":
            cpos = box.x0+box.width*0.55,box.y0+box.height+gap,box.width*0.45,width
        elif loc1 == "left":
            cpos = box.x0+box.width*0.0,box.y0+box.height+gap,box.width*0.45,width
        else:
            cpos = box.x0,box.y0+box.height+gap,box.width,width
    ax1 = fig.add_axes(cpos)
    #cmap = mpl.cm.jet
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    if loc0 == "left" or loc0 == "right":
        ori = "vertical"
    else:
        ori = "horizontal"
    if ticks is not None:
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                        norm=norm,
                                        orientation=ori,
                                        ticks=ticks,
                                        extend=extend)
    else:
        cb1 = mpl.colorbar.ColorbarBase(ax1,  cmap=cmap,
                                        norm=norm,
                                        orientation=ori,
                                        extend=extend)
    if loc0 == "top" or loc0 == "bottom":
        cb1.ax.xaxis.set_ticks_position(loc0)
        mathify_axes_ticks(cb1.ax,fontsize=fntsize,xticks=ticks)
    elif loc0 == "left" or loc0 == "right":
        cb1.ax.yaxis.set_ticks_position(loc0)
        mathify_axes_ticks(cb1.ax,fontsize=fntsize,yticks=ticks)
    if label is not None:
        ax1.set_xlabel(label,fontsize=fntsize)
    if ylabel is not None:
        ax1.set_ylabel(ylabel,fontsize=fntsize)
    for label in ax1.get_xticklabels() + ax1.get_yticklabels():
        label.set_fontsize(fntsize)
    return ax1,cb1

def Qmri(dir=2):
    """
    APPROXIMATELY Computes number of theta cells resolving one MRI wavelength
    """
    global bu,rho,uu,_dx2,_dx3
    #cvel()
    #corrected this expression to include both 2pi and dxdxp[3][3]
    #also corrected defition of va^2 to contain bsq+gam*ug term
    #need to figure out how to properly measure this in fluid frame
    vaudir = np.abs(bu[dir])/np.sqrt(rho+bsq+gam*ug)
    omega = dxdxp[3][3]*uu[3]/uu[0]+1e-15
    lambdamriudir = 2*np.pi * vaudir / omega
    if dir == 2:
        res=lambdamriudir/_dx2
    elif dir == 3:
        res=lambdamriudir/_dx3
    return(res)

def goodlabs(fntsize=20):
    ax = plt.gca()
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(fntsize)



def iofr(rval):
    rval = np.array(rval)
    if np.max(rval) < r[0,0,0]:
        return 0
    res = interp1d(r[:,0,0], ti[:,0,0], kind='linear', bounds_error = False, fill_value = 0)(rval)
    if len(res.shape)>0 and len(res)>0:
        res[rval<r[0,0,0]]*=0
        res[rval>r[nx-1,0,0]]=res[rval>r[nx-1,0,0]]*0+nx-1
    else:
        res = np.float64(res)
    return(np.floor(res+0.5).astype(int))


def plco(myvar,**kwargs):
    global r,h,ph
    plt.clf()
    return plc(myvar,**kwargs)

def plc(myvar,**kwargs): #plc
    global r,h,ph
    #xcoord = kwargs.pop('x1', None)
    #ycoord = kwargs.pop('x2', None)
    if(np.min(myvar)==np.max(myvar)):
        print("The quantity you are trying to plot is a constant = %g." % np.min(myvar))
        return
    cb = kwargs.pop('cb', False)
    nc = kwargs.pop('nc', 15)
    k = kwargs.pop('k',0)
    mirrorx = kwargs.pop('mirrorx',0)
    mirrory = kwargs.pop('mirrory',0)
    symmx = kwargs.pop('symmx',0)
    #cmap = kwargs.pop('cmap',cm.jet)
    isfilled = kwargs.pop('isfilled',False)
    xy = kwargs.pop('xy',0)
    xcoord = kwargs.pop("xcoord",None)
    ycoord = kwargs.pop("ycoord",None)
    lin = kwargs.pop('lin',0)
    xmax = kwargs.pop('xmax',10)
    ymax = kwargs.pop('ymax',5)
    cbxlabel = kwargs.pop('cbxla',None)
    cbylabel = kwargs.pop('cbyla',None)
    fntsize = kwargs.pop("fntsize",20)
    cbgoodticks = kwargs.pop("cbgoodticks",1)
    xlabel = kwargs.pop("xla",None)
    ylabel = kwargs.pop("yla",None)
    dobh = kwargs.pop("dobh",1)
    pretty = kwargs.pop("pretty",0)
    ax = kwargs.pop("ax",None)
    cbticks = kwargs.pop("cbticks",None)
    domathify = kwargs.pop("domathify",0)
    if np.abs(xy)==1:
        if xcoord is None: xcoord = r * np.sin(h)
        if ycoord is None: ycoord = r * np.cos(h)
        if mirrory: ycoord *= -1
        if mirrorx: xcoord *= -1
    if xcoord is not None and ycoord is not None:
        xcoord = xcoord[:,:,None] if xcoord.ndim == 2 else xcoord[:,:,k:k+1]
        ycoord = ycoord[:,:,None] if ycoord.ndim == 2 else ycoord[:,:,k:k+1]
    if np.abs(xy)==1 and symmx:
        if myvar.ndim == 2:
            myvar = myvar[:,:,None] if myvar.ndim == 2 else myvar[:,:,k:k+1]
            myvar=np.concatenate((myvar[:,::-1],myvar),axis=1)
            xcoord=np.concatenate((-xcoord[:,::-1],xcoord),axis=1)
            ycoord=np.concatenate((ycoord[:,::-1],ycoord),axis=1)
        else:
            if myvar.shape[-1] > 1: 
                symmk = (k+nz/2)%nz 
            else: 
                symmk = k
            myvar=np.concatenate((myvar[:,ny-1:ny,k:k+1],myvar[:,::-1,symmk:symmk+1],myvar[:,:,k:k+1]),axis=1)
            xcoord=np.concatenate((xcoord[:,ny-1:ny,k:k+1],-xcoord[:,::-1],xcoord),axis=1)
            ycoord=np.concatenate((ycoord[:,ny-1:ny,k:k+1],ycoord[:,::-1],ycoord),axis=1)
    elif np.abs(xy) == 2 and symmx:
        #if fracphi == 0.5 done in a robust way
        if get_fracphi() < 0.75:
            r1 = np.concatenate((r,r,r[...,0:1]),axis=2)
            ph1 = np.concatenate((ph,ph+np.pi,ph[...,0:1]+2*np.pi),axis=2)
            myvar = np.concatenate((myvar,myvar,myvar[...,0:1]),axis=2)
        else:
            r1 = np.concatenate((r,r[...,0:1]),axis=2)
            ph1 = np.concatenate((ph,ph[...,0:1]+2*np.pi),axis=2)
            myvar = np.concatenate((myvar,myvar[...,0:1]),axis=2)
        xcoord=(r1*cos(ph1))[:,ny/2,:,None]
        ycoord=(r1*sin(ph1))[:,ny/2,:,None]
        myvar = myvar[:,ny/2,:,None]
    else:
        myvar = myvar[:,:,None] if myvar.ndim == 2 else myvar[:,:,k:k+1]
    if lin:
        xcoord = r
        ycoord = h
    if ax is None:
        ax = plt.gca()
    if  xcoord is None or ycoord is None:
        if isfilled:
            res = ax.contourf(myvar[:,:,0].transpose(),nc,**kwargs)
        else:
            res = ax.contour(myvar[:,:,0].transpose(),nc,**kwargs)
    else:
        if isfilled:
            res = ax.contourf(xcoord[:,:,0],ycoord[:,:,0],myvar[:,:,0],nc,**kwargs)
        else:
            res = ax.contour(xcoord[:,:,0],ycoord[:,:,0],myvar[:,:,0],nc,**kwargs)
    if xy>0 and not symmx:
        ax.set_xlim(0,xmax)
        ax.set_ylim(-ymax,ymax)
    if xy> 0 and symmx:
        ax.set_xlim(-xmax,xmax)
        ax.set_ylim(-ymax,ymax)
    if xlabel is not None:
        ax.set_xlabel(xlabel,fontsize=fntsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel,fontsize=fntsize)
    if pretty:
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(fntsize)
            if domathify: mathify_axes_ticks(ax,fontsize=fntsize)
    if cb: #use color bar
        cb = plt.colorbar(res,ax=ax)
        if pretty and cbgoodticks and cbticks is None:
            vmin = cb.vmin
            vmax = cb.vmax
            #this returns incorrect ticks! so ignore it
            #ticks = cb.ax.get_yticks()
            #nticks = len(ticks)
            #if not too many ticks, then pretty them up
            rvmin = np.round(vmin)
            rvmax = np.round(vmax)
            if rvmin == vmin and rvmax == vmax and vmax-vmin <= 10:
                ticks = np.arange(rvmin,rvmax+1)
                cb.set_ticks(ticks)
                mathify_axes_ticks(cb.ax,fontsize=fntsize,yticks=ticks)
            elif rvmin == vmin and rvmax == vmax and vmax-vmin <= 20:
                ticks = np.arange(rvmin,rvmax+1)[::2]
                cb.set_ticks(ticks)
                mathify_axes_ticks(cb.ax,fontsize=fntsize,yticks=ticks)
        if cbticks is not None:
            cb.set_ticks(cbticks)
            mathify_axes_ticks(cb.ax,fontsize=fntsize,yticks=cbticks)
        if cbxlabel is not None:
            cb.ax.set_xlabel(cbxlabel,fontsize=fntsize)
        if cbxlabel is not None:
            cb.ax.set_xlabel(cbxlabel,fontsize=fntsize)
        if cbylabel is not None:
            cb.ax.set_ylabel(cbylabel,fontsize=fntsize)
        if pretty:
            for label in cb.ax.get_yticklabels():
                label.set_fontsize(fntsize)
    if xy and dobh and "rhor" in globals(): 
        el = Ellipse((0,0), 2*rhor, 2*rhor, facecolor='k', alpha=1)
        art=ax.add_artist(el)
        art.set_zorder(20)
    if cb:
        return res, cb
    else:
        return res

def faraday():
    global omegaf1, omegaf2
    if 'omegaf1' in globals():
        del omegaf1
    if 'omemaf2' in globals():
        del omegaf2
    omegaf1=fFdd(0,1)/fFdd(1,3)
    omegaf2=fFdd(0,2)/fFdd(2,3)

def Tcalcud():
    global Tud, TudEM, TudMA
    global mu, sigma
    global enth
    global unb, isunbound
    pg = (gam-1)*ug
    w=rho+ug+pg
    eta=w+bsq
    if 'Tud' in globals():
        del Tud
    if 'TudMA' in globals():
        del TudMA
    if 'TudEM' in globals():
        del TudEM
    if 'mu' in globals():
        del mu
    if 'sigma' in globals():
        del sigma
    if 'unb' in globals():
        del unb
    if 'isunbound' in globals():
        del isunbound
    Tud = np.zeros((4,4,nx,ny,nz),dtype=np.float32,order='F')
    TudMA = np.zeros((4,4,nx,ny,nz),dtype=np.float32,order='F')
    TudEM = np.zeros((4,4,nx,ny,nz),dtype=np.float32,order='F')
    for kapa in np.arange(4):
        for nu in np.arange(4):
            if(kapa==nu): delta = 1
            else: delta = 0
            TudEM[kapa,nu] = bsq*uu[kapa]*ud[nu] + 0.5*bsq*delta - bu[kapa]*bd[nu]
            TudMA[kapa,nu] = w*uu[kapa]*ud[nu]+pg*delta
            #Tud[kapa,nu] = eta*uu[kapa]*ud[nu]+(pg+0.5*bsq)*delta-bu[kapa]*bd[nu]
            Tud[kapa,nu] = TudEM[kapa,nu] + TudMA[kapa,nu]
    mu = -Tud[1,0]/(rho*uu[1])
    sigma = TudEM[1,0]/TudMA[1,0]
    enth=1+ug*gam/rho
    unb=enth*ud[0]
    isunbound=(-unb>1.0)


def aux():
    faraday()
    Tcalcud()

if __name__ == "__main__":
    if False:
        #1D plot example
        plt.clf()
        rg("gdump")
        rd("dump000")
        plt.plot(r[:,ny/2,0],rho[:,ny/2,0])
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("r")
        plt.ylabel("rho")
    if False:
        #2D plot example
        plt.clf()
        rg("gdump")
        rd("dump000")
        #R-z plot of the logarithm of density distribution
        plc(r,np.log10(rho),cb=True,xy=1,xmax=100,ymax=50)


def bhole():

    ax = plt.gca()
    el = Ellipse((0,0), 2*rhor, 2*rhor, facecolor='k', alpha=1)
    art=ax.add_artist(el)
    art.set_zorder(20)
    plt.draw()



def testfail(fldname = "dump000"):
    try: 
        rd(fldname)
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)


def get_sorted_file_list(prefix="dump"):
    flist0 = np.sort(glob.glob( os.path.join("dumps/", "%s[0-9][0-9][0-9]"%prefix) ) )
    flist1 = np.sort(glob.glob( os.path.join("dumps/", "%s[0-9][0-9][0-9][0-9]"%prefix) ) )
    flist2 = np.sort(glob.glob( os.path.join("dumps/", "%s[0-9][0-9][0-9][0-9][0-9]"%prefix) ) )
    flist0.sort()
    flist1.sort()
    flist2.sort()
    flist = np.concatenate((flist0,flist1,flist2))
    return flist
    


def fFdd(i,j):
    if i==0 and j==1:
        fdd =  gdet*(uu[2]*bu[3]-uu[3]*bu[2]) # f_tr
    elif i==1 and j==0:
        fdd = -gdet*(uu[2]*bu[3]-uu[3]*bu[2]) # -f_tr
    elif i==0 and j==2:
        fdd =  gdet*(uu[3]*bu[1]-uu[1]*bu[3]) # f_th
    elif i==2 and j==0:
        fdd = -gdet*(uu[3]*bu[1]-uu[1]*bu[3]) # -f_th
    elif i==0 and j==3:
        fdd =  gdet*(uu[1]*bu[2]-uu[2]*bu[1]) # f_tp
    elif i==3 and j==0:
        fdd = -gdet*(uu[1]*bu[2]-uu[2]*bu[1]) # -f_tp
    elif i==1 and j==3:
        fdd =  gdet*(uu[2]*bu[0]-uu[0]*bu[2]) # f_rp = gdet*B2
    elif i==3 and j==1:
        fdd = -gdet*(uu[2]*bu[0]-uu[0]*bu[2]) # -f_rp = gdet*B2
    elif i==2 and j==3:
        fdd =  gdet*(uu[0]*bu[1]-uu[1]*bu[0]) # f_hp = gdet*B1
    elif i==3 and j==2:
        fdd = -gdet*(uu[0]*bu[1]-uu[1]*bu[0]) # -f_hp = gdet*B1
    elif i==1 and j==2:
        fdd =  gdet*(uu[0]*bu[3]-uu[3]*bu[0]) # f_rh = gdet*B3
    elif i==2 and j==1:
        fdd = -gdet*(uu[0]*bu[3]-uu[3]*bu[0]) # -f_rh = gdet*B3
    else:
        fdd = np.zeros_like(uu[0])
    return fdd

delta = lambda kapa,nu: (kapa==nu)
fTudEM = lambda kapa,nu: bsq*uu[kapa]*ud[nu] + 0.5*bsq*delta(kapa,nu) - bu[kapa]*bd[nu]
fTudMA = lambda kapa,nu: (rho+gam*ug)*uu[kapa]*ud[nu]+(gam-1)*ug*delta(kapa,nu)
fTud = lambda kapa,nu: fTudEM(kapa,nu) + fTudMA(kapa,nu)
fRud = lambda kapa,nu: 4./3.*Erf*uradu[kapa]*uradd[nu]+1./3.*Erf*delta(kapa,nu)

def odot(a,b):
    """ Outer product of two vectors a^mu b_nu"""
    #the shape of the product is (4,4,nx,ny,max(a.nz,b.nz))
    outer_product = np.zeros(np.concatenate((np.array((4,4)),amax(a[0].shape,b[0].shape))),dtype=np.float32,order='F')
    for mu in np.arange(4):
        for nu in np.arange(4):
            outer_product[mu,nu] = a[mu]*b[nu]
    return(outer_product)


def amax(arg1,arg2):
    return(np.maximum(arg1,arg2))

def amin(arg1,arg2):
    return(np.minimum(arg1,arg2))


#############################
#
# End of movie making
#
#############################
        
    
if __name__ == "__main__":
    if len(sys.argv)>1:
        if sys.argv[1].startswith("mkfrm"):
            mkmov_wrapper(which=sys.argv[1])
        elif sys.argv[1].startswith("convertdump"):
            convert_wrapper(which=sys.argv[1])
        else:
            print( "Unknown command %s" % sys.argv[1] )
    else:
        if False:
            #1D plot example
            plt.clf()
            rg("gdump")
            rd("dump000")
            plt.loglog(r[:,ny/2,0],rho[:,ny/2,0])
            plt.xlabel("r")
            plt.ylabel("rho")
            plt.savefig("dumptest.eps")
        if True:
            #2D plot example
            plt.clf()
            rg("gdump")
            rd("dump000")
            #R-z plot of the logarithm of density distribution
            plc(np.log10(rho),cb=True,xy=1,xmax=50,ymax=30, isfilled=True)
            aphi=psicalc()
            plc(aphi, xy=1,xmax=50,ymax=30, colors='k')
            plt.savefig("inittest.eps")

