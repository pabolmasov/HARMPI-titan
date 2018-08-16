from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import sys
import glob
import os

from numpy import sin, cos, tan, pi
import numpy as np

# makes a dump-file name from its number
# probably, could be done in one line...
def dumpname(n):
    s=str(n)
    if(n>0):
        ls=int(np.floor(np.log10(np.double(n))))
    else:
        ls=0
    if(ls<2):
        s='0'*(2-ls)+s
    s='dump'+s
    return s

#read in a dump file
def rd(dump):
    read_file(dump,type="dump")


#read in a grid file
def rg(dump):
    read_file(dump,type="gdump")

#read in a grid file
def rg2(dump):
    read_file(dump,type="gdump2",noround=True)

#high-level function that reads either MPI or serial gdump's
def read_file(dump,type=None,savedump=True,saverdump=False,noround=False):
    if type is None:
        if dump.startswith("dump"):
            type = "dump"
            print("Reading a dump file %s ..." % dump)
        elif dump.startswith("gdump2"):
            type = "gdump2"
            print("Reading a gdump2 file %s ..." % dump)
        elif dump.startswith("gdump"):
            type = "gdump"
            print("Reading a gdump file %s ..." % dump)
        elif dump.startswith("rdump"):
            type = "rdump"
            print("Reading a rdump file %s ..." % dump)
        elif dump.startswith("fdump"):
            type = "fdump"
            print("Reading a fdump file %s ..." % dump)
        else:
            print("Couldn't guess dump type; assuming it is a data dump")
            type = "dump"
    #normal dump
    if os.path.isfile( "dumps/" + dump ):
        headerline = read_header("dumps/" + dump, returnheaderline = True, issilent=False)
        gd = read_body("dumps/" + dump,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G,noround=1)
        if noround:
            res = data_assign(         gd,type=type,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G)
        else:
            res = data_assign(myfloat(gd),type=type,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G)
        return res
    #MPI-type dump that is spread over many files
    else:
        flist = np.sort(glob.glob( "dumps/" + dump + "_[0-9][0-9][0-9][0-9]" ))
        if len(flist) == 0:
            print( "Could not find %s or its MPI counterpart" % dump )
            return
        sys.stdout.write( "Reading %s (%d files)" % (dump, len(flist)) )
        sys.stdout.flush()
        ndots = 10
        dndot = old_div(len(flist),ndots)
        if dndot == 0: dndot = 1
        for i,fname in enumerate(flist):
            #print( "Reading file %d out of %d..." % (i,len(flist)) )
            #header for each file might be different, so read each
            header = read_header(fname,issilent=1)
            if header is None:
                print( "Error reading header of %s, aborting..." % fname )
                return
            lgd = read_body(fname,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G)
            #this gives an array of dimensions (-1,N1,N2,N3)+potentially ghost cells
            if 0 == i:
                #create full array: of dimensions, (-1,nx,ny,nz)
                fgd = np.zeros( (lgd.shape[0], nx+2*N1G, ny+2*N2G, nz+2*N3G), dtype=np.float32)
            if not type == "rdump":
                #construct full indices: ti, tj, tk
                #fti,ftj,ftk = mgrid[0:nx,0:ny,0:nz]
                lti,ltj,ltk = lgd[0:3,:,:].view();
                lti = np.int64(lti)
                ltj = np.int64(ltj)
                ltk = np.int64(ltk)
                fgd[:,lti+N1G,ltj+N2G,ltk+N3G] = lgd[:,:,:,:]
            else:
                print(starti,startj,startk)
                fgd[:,starti:starti+N1+2*N1G,startj:startj+N2+2*N2G,startk:startk+N3+2*N3G] = lgd[:,:,:,:]
            del lgd
            if i%dndot == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
        res = data_assign(fgd,type=type,nx=nx+2*N1G,ny=ny+2*N2G,nz=nz+2*N3G)
        if savedump:
            #if the full dump file does not exist, create it
            dumpfullname = "dumps/" + dump
            if (type == "dump") and not os.path.isfile(dumpfullname):
                sys.stdout.write("Saving full dump to %s..." % dumpfullname)
                sys.stdout.flush()
                header[1] = header[4] #N1 = nx
                header[2] = header[5] #N2 = ny
                header[3] = header[6] #N3 = nz
                fout = open( dumpfullname, "wb" )
                #join header items with " " (space) as a glue
                #see http://stackoverflow.com/questions/12377473/python-write-versus-writelines-and-concatenated-strings
                #write it out with a new line char at the end
                fout.write(b' '.join(header) + b'\n')
                fout.flush()
                os.fsync(fout.fileno())
                #reshape the dump content
                gd1 = fgd.transpose(1,2,3,0)
                gd1.tofile(fout)
                fout.close()
                print( " done!" )
                if res is not None:
                    return res
        return res

#read in a header
def read_header(dump,issilent=True,returnheaderline=False):
    global t,nx,ny,nz,N1,N2,N3,N1G,N2G,N3G,starti,startj,startk,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,ti,tj,tk,x1,x2,x3,r,h,ph,gcov,gcon,gdet,drdx,gn3,gv3,guu,gdd,dxdxp, games, startx1, startx2, startx3, tf, NPR, DOKTOT, BL
    global fractheta
    global fracphi
    global rbr
    global npow2
    global cpow2
    #read image
    fin = open( dump, "rb" )
    headerline = fin.readline()
    header = headerline.split()
    nheadertot = len(header)
    fin.close()
    if not dump.startswith("dumps/rdump"):
        if not issilent: print( "dump header: len(header) = %d" % len(header) )
        nheader = 45
        n = 0
        t = myfloat(np.float64(header[n])); n+=1
        #per tile resolution
        N1 = int(header[n]); n+=1
        N2 = int(header[n]); n+=1
        N3 = int(header[n]); n+=1
        #total resolution
        nx = int(header[n]); n+=1
        ny = int(header[n]); n+=1
        nz = int(header[n]); n+=1
        #numbers of ghost cells
        N1G = int(header[n]); n+=1
        N2G = int(header[n]); n+=1
        N3G = int(header[n]); n+=1
        startx1 = myfloat(float(header[n])); n+=1
        startx2 = myfloat(float(header[n])); n+=1
        startx3 = myfloat(float(header[n])); n+=1
        _dx1=myfloat(float(header[n])); n+=1
        _dx2=myfloat(float(header[n])); n+=1
        _dx3=myfloat(float(header[n])); n+=1
        tf=myfloat(float(header[n])); n+=1
        nstep=myfloat(float(header[n])); n+=1
        a=myfloat(float(header[n])); n+=1
        gam=myfloat(float(header[n])); n+=1
        cour=myfloat(float(header[n])); n+=1
        DTd=myfloat(float(header[n])); n+=1
        DTl=myfloat(float(header[n])); n+=1
        DTi=myfloat(float(header[n])); n+=1
        DTr=myfloat(float(header[n])); n+=1
        DTr01=myfloat(float(header[n])); n+=1
        dump_cnt=myfloat(float(header[n])); n+=1
        image_cnt=myfloat(float(header[n])); n+=1
        rdump_cnt=myfloat(float(header[n])); n+=1
        rdump01_cnt=myfloat(float(header[n])); n+=1
        dt=myfloat(float(header[n])); n+=1
        lim=myfloat(float(header[n])); n+=1
        failed=myfloat(float(header[n])); n+=1
        Rin=myfloat(float(header[n])); n+=1
        Rout=myfloat(float(header[n])); n+=1
        hslope=myfloat(float(header[n])); n+=1
        R0=myfloat(float(header[n])); n+=1
        NPR=int(header[n]); n+=1
        DOKTOT=int(header[n]); n+=1
        fractheta = myfloat(header[n]); n+=1
        fracphi   = myfloat(header[n]); n+=1
        rbr       = myfloat(header[n]); n+=1
        npow2     = myfloat(header[n]); n+=1
        cpow2     = myfloat(header[n]); n+=1
        BL = myfloat(header[n]); n+=1
    else:
        print("rdump header")
        nheader = 46
        n = 0
        #per tile resolution
        N1 = int(header[n]); n+=1
        N2 = int(header[n]); n+=1
        N3 = int(header[n]); n+=1
        #total resolution
        nx = int(header[n]); n+=1
        ny = int(header[n]); n+=1
        nz = int(header[n]); n+=1
        #numbers of ghost cells
        N1G = int(header[n]); n+=1
        N2G = int(header[n]); n+=1
        N3G = int(header[n]); n+=1
        #starting indices
        starti = int(header[n]); n+=1
        startj = int(header[n]); n+=1
        startk = int(header[n]); n+=1
        t = myfloat(header[n]); n+=1
        tf = myfloat(header[n]); n+=1
        nstep = int(header[n]); n+=1
        a = myfloat(header[n]); n+=1
        gam = myfloat(header[n]); n+=1
        game = myfloat(header[n]); n+=1
        game4 = myfloat(header[n]); n+=1
        game5 = myfloat(header[n]); n+=1
        cour = myfloat(header[n]); n+=1
        DTd = myfloat(header[n]); n+=1
        DTl = myfloat(header[n]); n+=1
        DTi = myfloat(header[n]); n+=1
        DTr = myfloat(header[n]); n+=1
        DTr01 = myfloat(header[n]); n+=1
        dump_cnt = myfloat(header[n]); n+=1
        image_cnt = myfloat(header[n]); n+=1
        rdump_cnt = myfloat(header[n]); n+=1
        rdump01_cnt=myfloat(float(header[n])); n+=1
        dt = myfloat(header[n]); n+=1
        lim = myfloat(header[n]); n+=1
        failed = myfloat(header[n]); n+=1
        Rin = myfloat(header[n]); n+=1
        Rout = myfloat(header[n]); n+=1
        hslope = myfloat(header[n]); n+=1
        R0 = myfloat(header[n]); n+=1
        fractheta = myfloat(header[n]); n+=1
        fracphi = myfloat(header[n]); n+=1
        rbr = myfloat(header[n]); n+=1
        npow2 = myfloat(header[n]); n+=1
        cpow2 = myfloat(header[n]); n+=1
        tdump = myfloat(header[n]); n+=1
        trdump = myfloat(header[n]); n+=1
        timage = myfloat(header[n]); n+=1
        tlog  = myfloat(header[n]); n+=1
    if n != nheader or n != nheadertot:
        print("Wrong number of elements in header: nread = %d, nexpected = %d, nototal = %d: incorrect format?"
              % (n, nheader, nheadertot) )
        return headerline
    if returnheaderline:
        return headerline
    else:
        return header
            
def read_body(dump,nx=None,ny=None,nz=None,noround=False):
        fin = open( dump, "rb" )
        header = fin.readline()
        if dump.startswith("dumps/rdump"):
            dtype = np.float64
            body = np.fromfile(fin,dtype=dtype,count=-1)
            gd = body.view().reshape((nx,ny,nz,-1), order='C')
            if noround:
                gd=gd.transpose(3,0,1,2)
            else:
                gd=myfloat(gd.transpose(3,0,1,2))
        elif dump.startswith("dumps/gdump2"):
            dtype = np.float64
            body = np.fromfile(fin,dtype=dtype,count=-1)
            gd = body.view().reshape((nx,ny,nz,-1), order='C')
            if noround:
                gd=gd.transpose(3,0,1,2)
            else:
                gd=myfloat(gd.transpose(3,0,1,2))
        elif dump.startswith("dumps/fdump"):
            dtype = np.int64
            body = np.fromfile(fin,dtype=dtype,count=-1)
            gd = body.view().reshape((-1,nz,ny,nx), order='F')
            gd=myfloat(gd.transpose(0,3,2,1))
        else:
            dtype = np.float32
            body = np.fromfile(fin,dtype=dtype,count=-1)
            gd = body.view().reshape((-1,nz,ny,nx), order='F')
            gd=myfloat(gd.transpose(0,3,2,1))
        return gd

def data_assign(gd,type=None,**kwargs):
    if type is None:
        print("Please specify data type")
        return
    if type == "gdump":
        gdump_assign(gd,**kwargs)
        return None
    elif type == "gdump2":
        gdump2_assign(gd,**kwargs)
        return None
    elif type == "dump":
        dump_assign(gd,**kwargs)
        return None
    elif type == "rdump":
        gd = rdump_assign(gd,**kwargs)
        return gd
    elif type == "fdump":
        gd = fdump_assign(gd,**kwargs)
        return gd
    else:
        print("Unknown data type: %s" % type)
        return gd
    
def gdump_assign(gd,**kwargs):
    global t,nx,ny,nz,N1,N2,N3,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,ti,tj,tk,x1,x2,x3,r,h,ph,gcov,gcon,gdet,drdx,gn3,gv3,guu,gdd,dxdxp, games
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    ti,tj,tk,x1,x2,x3,r,h,ph = gd[0:9,:,:].view();  n = 9
    gv3 = gd[n:n+16].view().reshape((4,4,nx,ny,nz),order='F').transpose(1,0,2,3,4); n+=16
    gn3 = gd[n:n+16].view().reshape((4,4,nx,ny,nz),order='F').transpose(1,0,2,3,4); n+=16
    gcov = gv3
    gcon = gn3
    guu = gn3
    gdd = gv3
    gdet = gd[n]; n+=1
    drdx = gd[n:n+16].view().reshape((4,4,nx,ny,nz),order='F').transpose(1,0,2,3,4); n+=16
    dxdxp = drdx
    if n != gd.shape[0]:
        print("rd: WARNING: nread = %d < ntot = %d: incorrect format?" % (n, gd.shape[0]) )
        return 1
    return 0

def gdump2_assign(gd,**kwargs):
    global t,nx,ny,nz,N1,N2,N3,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,ti,tj,tk,x1,x2,x3,gdet,games,rf1,hf1,phf1,rf2,hf2,phf2,rf3,hf3,phf3,rcorn,hcord,phcorn,re1,he1,phe1,re2,he2,phe2,re3,he3,phe3
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    ti,tj,tk,x1,x2,x3 = gd[0:6,:,:].view();  n = 6
    rf1,hf1,phf1,rf2,hf2,phf2,rf3,hf3,phf3 = gd[0:9,:,:].view();  n += 9
    rcorn,hcord,phcorn,rcent,hcent,phcen = gd[0:6,:,:].view();  n += 6
    re1,he1,phe1,re2,he2,phe2,re3,he3,phe3 = gd[0:9,:,:].view();  n += 9
    gdet = gd[n]; n+=1
    if n != gd.shape[0]:
        print("rd: WARNING: nread = %d < ntot = %d: incorrect format?" % (n, gd.shape[0]) )
        return 1
    return 0

#read in a dump file
def dump_assign(gd,**kwargs):
    global t,nx,ny,nz,_dx1,_dx2,_dx3,gam,hslope,a,R0,Rin,Rout,ti,tj,tk,x1,x2,x3,r,h,ph,rho,ug,vu,B,pg,cs2,Sden,U,gdetB,divb,uu,ud,bu,bd,v1m,v1p,v2m,v2p,gdet,bsq,gdet,alpha,rhor, ktot, pg, origin_r, origin_th, origin_phi
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    ti,tj,tk,x1,x2,x3,r,h,ph,rho,ug = gd[0:11,:,:].view(); n = 11
    pg = (gam-1)*ug
    lrho=np.log10(rho)
    vu=np.zeros_like(gd[0:4])
    B=np.zeros_like(gd[0:4])
    vu[1:4] = gd[n:n+3]; n+=3
    B[1:4] = gd[n:n+3]; n+=3
    origin_r = gd[n]/rho ; n+=1
    origin_th = gd[n]/rho ; n+=1
    origin_phi = gd[n]/rho ; n+=1
    #if total entropy equation is evolved (on by default)
    if DOKTOT == 1:
      ktot = gd[n]; n+=1
    divb = gd[n]; n+=1
    uu = gd[n:n+4]; n+=4
    ud = gd[n:n+4]; n+=4
    bu = gd[n:n+4]; n+=4
    bd = gd[n:n+4]; n+=4
    bsq = mdot(bu,bd)
    v1m,v1p,v2m,v2p,v3m,v3p=gd[n:n+6]; n+=6
    gdet=gd[n]; n+=1
    rhor = 1+(1-a**2)**0.5
    if "guu" in globals():
        #lapse
        alpha = (-guu[0,0])**(-0.5)
    if n != gd.shape[0]:
        print("rd: WARNING: nread = %d < ntot = %d: incorrect format?" % (n, gd.shape[0]) )
        return 1
    return 0

def rdump_assign(gd,**kwargs):
    global t,nx,ny,nz,_dx1,_dx2,_dx3,gam,hslope,a,R0,Rin,Rout,ti,tj,tk,x1,x2,x3,r,h,ph,rho,ug,vu,B,pg,cs2,Sden,U,gdetB,divb,uu,ud,bu,bd,v1m,v1p,v2m,v2p,gdet,bsq,gdet,alpha,rhor, ktot, Ttot, game, qisosq, pflag, qisodotb, kel, uelvar, Tel4, Tel5,Teldis, Tels, kel4, kel5,ugel,ugeldis, ugcon, sel, ugscon, ugel4, ugel5,stot, uelvar, Telvar, Tsel, sel, ugels, games, phi, keldis, phihat,csphib,lrho
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    n = 0
    rho = gd[n]; n+=1
    ug = gd[n]; n+=1
    vu=np.zeros_like(gd[0:4])
    B=np.zeros_like(gd[0:4])
    vu[1:4] = gd[n:n+3]; n+=3
    B[1:4] = gd[n:n+3]; n+=3
    # if n != gd.shape[0]:
    #     print("rd: WARNING: nread = %d < ntot = %d: incorrect format?" % (n, gd.shape[0]) )
    #     return 1
    return gd

def fdump_assign(gd,**kwargs):
    global t,nx,ny,nz,_dx1,_dx2,_dx3,gam,hslope,a,R0,Rin,Rout,ti,tj,tk,x1,x2,x3,r,h,ph,rho,ug,vu,B,pg,cs2,Sden,U,gdetB,divb,uu,ud,bu,bd,v1m,v1p,v2m,v2p,gdet,bsq,gdet,alpha,rhor, ktot, Ttot, game, qisosq, pflag, qisodotb, kel, uelvar, Tel4, Tel5,Teldis, Tels, kel4, kel5,ugel,ugeldis, ugcon, sel, ugscon, ugel4, ugel5,stot, uelvar, Telvar, Tsel, sel, ugels, games, phi, keldis, phihat,csphib,lrho,fail
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    fail = gd
    return gd



def mdot(a,b):
    """
    Computes a contraction of two tensors/vectors.  Assumes
    the following structure: tensor[m,n,i,j,k] OR vector[m,i,j,k], 
    where i,j,k are spatial indices and m,n are variable indices. 
    """
    if (a.ndim == 3 and b.ndim == 3) or (a.ndim == 4 and b.ndim == 4):
          c = (a*b).sum(0)
    elif a.ndim == 5 and b.ndim == 4:
          c = np.empty(np.maximum(a[:,0,:,:,:].shape,b.shape),dtype=b.dtype)
          for i in range(a.shape[0]):
                c[i,:,:,:] = (a[i,:,:,:,:]*b).sum(0)
    elif a.ndim == 4 and b.ndim == 5:
          c = np.empty(np.maximum(b[0,:,:,:,:].shape,a.shape),dtype=a.dtype)
          for i in range(b.shape[1]):
                c[i,:,:,:] = (a*b[:,i,:,:,:]).sum(0)
    elif a.ndim == 5 and b.ndim == 5:
          c = np.empty((a.shape[0],b.shape[1],a.shape[2],a.shape[3],max(a.shape[4],b.shape[4])),dtype=a.dtype)
          for i in range(c.shape[0]):
                for j in range(c.shape[1]):
                      c[i,j,:,:,:] = (a[i,:,:,:,:]*b[:,j,:,:,:]).sum(0)
    elif a.ndim == 5 and b.ndim == 6:
          c = np.empty((a.shape[0],b.shape[1],b.shape[2],max(a.shape[2],b.shape[3]),max(a.shape[3],b.shape[4]),max(a.shape[4],b.shape[5])),dtype=a.dtype)
          for mu in range(c.shape[0]):
              for k in range(c.shape[1]):
                  for l in range(c.shape[2]):
                      c[mu,k,l,:,:,:] = (a[mu,:,:,:,:]*b[:,k,l,:,:,:]).sum(0)
    else:
           raise Exception('mdot', 'wrong dimensions')
    return c


def psicalc(B1=None):
    """
    Computes the field vector potential
    """
    global B
    if B1 is None: B1 = B[1]
    daphi = -(gdet*B1).mean(-1)*_dx2
    aphi=daphi[:,::-1].cumsum(axis=1)[:,::-1]
    aphi-=0.5*daphi #correction for half-cell shift between face and center in theta
    return(aphi)


def myfloat(f,acc=1):
    """ acc=1 means np.float32, acc=2 means np.float64 """
    if acc==1:
        return( np.float32(f) )
    else:
        return( np.float64(f) )

def get_fracphi():
    fracphi = dxdxp[3,3,0,0,0]*_dx3*nz/(2*np.pi)
    return( fracphi )

##########################
# reading a four-column text file:
def uread(infile, dims):
    fin=open(infile, 'r')
    s=str.split(str.strip(fin.readline()))
    u0=int(s[0]) ; u1=int(s[1]) ; u2=int(s[2]) ; u3=int(s[3]) 
    while(s):
        u0=int(s[0]) ; u1=int(s[1]) ; u2=int(s[2]) ; u3=int(s[3]) 
        s=str.split(str.strip(fin.readline()))
    fin.close()
    u0=reshape(asarray(u0, dtype=double), dims)
    u1=reshape(asarray(u1, dtype=double), dims)
    u2=reshape(asarray(u2, dtype=double), dims)
    u3=reshape(asarray(u3, dtype=double), dims)
    return u0, u1, u2, u3
