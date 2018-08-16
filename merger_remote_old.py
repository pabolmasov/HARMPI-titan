from numpy import *
import numpy as np
import time
import os
from os.path import *
from os import fsync
# import harm_script as h
import glob
import sys

# working directory
dire=''

# the challenge here is to get without plotting and ipython libraries

# harm_script routines:
#read in a dump file
def rd(dump):
    read_file(dump,type="dump")

#read in a grid file
def rg(dump):
    read_file(dump,type="gdump")

#read in a grid file
def rg2(dump):
    read_file(dump,type="gdump2",noround=True)

def psicalc(B1=None):
    """
    Computes the field vector potential
    """
    global B
    if B1 is None: B1 = B[1]
    daphi = -(gdet*B1).mean(-1)*_dx2
#    print shape(daphi)
    aphi=daphi[:,::-1]
    aphi=aphi.cumsum(axis=1)[:,::-1]
    aphi-=0.5*daphi #correction for half-cell shift between face and center in theta
    return(aphi)


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
        headerline = read_header("dumps/" + dump, returnheaderline = True)
        gd = read_body("dumps/" + dump,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G,noround=1)
        if noround:
            res = data_assign(gd,type=type,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G)
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
        dndot = len(flist)/ndots
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
                print starti,startj,startk
                fgd[:,starti:starti+N1+2*N1G,startj:startj+N2+2*N2G,startk:startk+N3+2*N3G] = lgd[:,:,:,:]
            del lgd
            if i%dndot == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
        res = data_assign(fgd,type=type,nx=nx+2*N1G,ny=ny+2*N2G,nz=nz+2*N3G)
        if savedump:
            #if the full dump file does not exist, create it
            dumpfullname = "dumps/" + dump
            if (type == "dump" or type == "gdump" or type == "rdump") and not os.path.isfile(dumpfullname):
                sys.stdout.write("Saving full dump to %s..." % dumpfullname)
                sys.stdout.flush()
                header[1] = header[4] #N1 = nx
                header[2] = header[5] #N2 = ny
                header[3] = header[6] #N3 = nz
                fout = open( dumpfullname, 'wb' )
                #join header items with " " (space) as a glue
                #see http://stackoverflow.com/questions/12377473/python-write-versus-writelines-and-concatenated-strings
                #write it out with a new line char at the end
                fout.write(" ".join(header) + "\n")
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
    print(dump)
    fin = open(dump, 'rb')   # here be integer
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
        fin = open( dump, 'rb' )
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
    global t,nx,ny,nz,_dx1,_dx2,_dx3,gam,hslope,a,R0,Rin,Rout,ti,tj,tk,x1,x2,x3,r,h,ph,rho,ug,vu,B,pg,cs2,Sden,U,gdetB,divb,uu,ud,bu,bd,v1m,v1p,v2m,v2p,gdet,bsq,gdet,alpha,rhor, ktot, Ttot, game, qisosq, pflag, qisodotb, kel, uelvar, Tel4, Tel5,Teldis, Tels, kel4, kel5,ugel,ugeldis, ugcon, sel, ugscon, ugel4, ugel5,stot, uelvar, Telvar, Tsel, sel, ugels, games, phi, keldis, phihat,csphib,lrho, Tnuc, Tnuc_cgs, etae, flr
    global kel4a, kel4b, kel4c, kel4d, kel4e, ugel4a, ugel4b, ugel4c, ugel4d, ugel4e
    global Tel4a, Tel4b, Tel4c, Tel4d, Tel4e 
    global vpot
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    ti,tj,tk,x1,x2,x3,r,h,ph,rho,ug = gd[0:11,:,:].view(); n = 11
    lrho=np.log10(rho)
    vu=np.zeros_like(gd[0:4])
    B=np.zeros_like(gd[0:4])
    vu[1:4] = gd[n:n+3]; n+=3
    B[1:4] = gd[n:n+3]; n+=3
    #if electrons are evolved
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

def myfloat(f,acc=1):
    """ acc=1 means np.float32, acc=2 means np.float64 """
    if acc==1:
        return( np.float32(f) )
    else:
        return( np.float64(f) )

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
#############################################################
# tetrad components... and all the grid parameters should be set as well
def tetrad_t(uumean, udmean):
    return -udmean[0]/drdx[0,0], -udmean[1]/drdx[1,1], 0.*udmean[0]/drdx[2,2], -udmean[3]/drdx[3,3]

def tetrad_r(uumean, udmean):
    gg=guu[0,3]**2-guu[0,0]*guu[3,3] #1.+uumean[3]*udmean[3]
    hh=guu[3,3]+uumean[3]*uumean[3]
    s=sqrt(fabs(gg*hh*guu[1,1]))
    return guu[3,3]*uumean[1]/s/drdx[0,0], gg*udmean[0]/s/drdx[1,1], 0.*uumean[0]/drdx[2,2], -guu[0,3]*uumean[1]/s/drdx[3,3]

def tetrad_h(uumean, udmean):
    return 0.*uumean[0]/drdx[0,0], 0.*uumean[0]/drdx[1,1], sqrt(gdd[2,2])/drdx[2,2], 0.*uumean[0]/drdx[3,3]

def tetrad_p(uumean, udmean):
    hh=sqrt(guu[3,3]+uumean[3]*uumean[3])
    return uumean[3]*udmean[0]/hh/drdx[0,0], uumean[3]*udmean[1]/hh/drdx[1,1], 0.*uumean[0]/drdx[2,2], (1.+uumean[3]*udmean[3])/hh/drdx[3,3]

# outputting all the basic information about the simulation setup to a file
def dumpinfo(prefix):
    
    rg("gdump")
    rd(prefix)
    fout=open(prefix+"_dinfo.dat", "w")

    fout.write(str(nx)+" "+str(ny)+" "+str(nz)+"\n")
    print "nr x ntheta x nphi = "+str(nx)+" "+str(ny)+" "+str(nz)
    
    fout.write(str(a)+"\n")
    print "Kerr parameter a = "+str(a)

    fout.write(str(t)+"\n")
    print "time "+str(t)
    fout.close()

# writes out evolution of some global parameters
def glevol(nmax, rref):
#    global rho, uu
   # rg("gdump")

    # rref is reference radius
    run=unique(r)
    nr=run[where(run<rref)].argmax()

#    maccre=zeros(nmax, dtype=double)

    fmdot=open("merge_mevol"+str(rref)+".dat", "w")

    for k in arange(nmax):
        fname=dumpname(k)
	print "reading "+str(fname)
        rd(fname)
	print "rho is "+str(shape(rho))
        #        Tcalcud()
        # accretion rate at rref
	rhom=rho.mean(axis=2)
	uum=uu.mean(axis=3)
	maccre=-trapz((rhom*uum[1]*(uum[1]<0.)*sqrt(gdet/gcov[1,1]))[nr,:], x=h[nr,:,0])*_dx2*_dx3
        mwind=-trapz((rhom*uum[1]*(uum[1]>0.)*sqrt(gdet/gcov[1,1]))[nr,:], x=h[nr,:,0])*_dx2*_dx3
        laccre=-trapz((rho*fabs(ud[3]/ud[0])*uu[1]*(uu[1]<0.)*sqrt(gdet/gcov[1,1]))[nr,:], x=h[nr,:,0])*_dx2*_dx3
	lwind=-trapz((rho*fabs(ud[3]/ud[0])*uu[1]*(uu[1]>0.)*sqrt(gdet/gcov[1,1]))[nr,:], x=h[nr,:,0])*_dx2*_dx3
        # maccre=-((rho*uu[1]*sqrt(gdet/gcov[1,1]))[nr,:,:]).sum()*_dx2*_dx3
        fmdot.write(str(t)+" "+str(maccre)+" "+str(mwind)+" "+str(laccre/maccre)+" "+str(lwind/mwind)+"\n")
    
    fmdot.close()

# makes a dump-file name from its number
# probably, could be done in one line...
def dumpname(n):
    s=str(n)
    if(n>0):
        ls=int(floor(log10(double(n))))
    else:
        ls=0
    if(ls<2):
        s='0'*(2-ls)+s
    s='dump'+s
    return s

def mint(rref):
    print "mint:"
    print "shape(rho) = "+str(shape(rho))
    print "shape(uu1) = "+str(shape(uu[1]))
    print "shape gdet = "+str(shape(gdet))
    hun=unique(h)
    run=unique(r)
    nr=run[where(run<rref)].argmax()
    indd=squeeze((rho*uu[1]*(uu[1]<0.)*sqrt(gdet/gcov[1,1]))[nr,:,:])
    maccre=-trapz((indd).mean(axis=-1), x=hun)*_dx2*_dx3
    indd=squeeze((rho*uu[1]*(uu[1]>0.)*sqrt(gdet/gcov[1,1]))[nr,:,:])
    mwind=trapz((indd).mean(axis=-1), x=hun)*_dx2*_dx3
    indd=squeeze((rho*uu[1]*(uu[1]<0.)*fabs(ud[3]/drdx[3,3])*sqrt(gdet/gcov[1,1]))[nr,:,:])
    laccre=-trapz((indd).mean(axis=-1), x=hun)*_dx2*_dx3
    indd=squeeze((rho*uu[1]*(uu[1]>0.)*fabs(ud[3]/drdx[3,3])*sqrt(gdet/gcov[1,1]))[nr,:,:])
    lwind=-trapz((indd).mean(axis=-1), x=hun)*_dx2*_dx3

    return maccre, mwind, laccre/maccre, lwind/mwind

def readndump(n1, n2, rref=5.0):

#    run=unique(r)
#    nr=run[where(run<rref)].argmax()

    rg("gdump")

    if (n2<n1):
        print("readndump: invalid file number range")
        exit()

    fmdot=open("merge_mevol"+str(rref)+".dat", "w")

    nframes=n2-n1+1
    n=n1+arange(nframes)

    for k in n:
        fname=dumpname(k)
        rd(fname)
        #        print(shape(rho))
        Tcalcud()
#	rhoaver=rho
	p=(gam-1.)*ug
	magp=bsq/2.
#        if(pflag):
#            uaver=(gam-1.)*ug
#	    if(withmag):
#	    	uaver+=bsq/2.
#        else:
#            uaver=rho
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

            tudem=TudEM ; tudma=TudMA
            pmean=(gam-1.)*ug
	   # unorm=uaver
	    magp_mean=magp
	    aphi=psicalc()
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

            pmean+=(gam-1.)*ug
	    magp_mean+=magp
	    aphi+=psicalc()
#	    unorm+=uaver
	    tudem+=TudEM ; tudma+=TudMA
	# mass accretion rate estimate
#        uum=uu.mean(axis=3)

#        rhom=rho.mean(axis=2)
#        uum=uu.mean(axis=3)
#        maccre=-trapz(((rho*uu[1]*(uu[1]<0.)*sqrt(gdet/gcov[1,1]))[nr,:]).mean(axis=1), x=h[nr,:,0])*_dx2*_dx3
#        mwind=trapz(((rho*uu[1]*(uu[1]>0.)*sqrt(gdet/gcov[1,1]))[nr,:]).mean(axis=1), x=h[nr,:,0])*_dx2*_dx3
#        laccre=-trapz(((rho*fabs(ud[3])*uu[1]*(uu[1]<0.)*sqrt(gdet/gcov[1,1]))[nr,:]).mean(axis=1), x=h[nr,:,0])*_dx2*_dx3
#        lwind=-trapz(((rho*fabs(ud[3])*uu[1]*(uu[1]>0.)*sqrt(gdet/gcov[1,1]))[nr,:]).mean(axis=1), x=h[nr,:,0])*_dx2*_dx3
        # maccre=-((rho*uu[1]*sqrt(gdet/gcov[1,1]))[nr,:,:]).sum()*_dx2*_dx3
	maccre, mwind, laccre, lwind = mint(rref)
        fmdot.write(str(t)+" "+str(maccre)+" "+str(mwind)+" "+str(laccre/maccre)+" "+str(lwind/mwind)+"\n")
    fmdot.close()
    # velocity normalization:
    uu0*=drdx[0,0]/rhomean ; uur*=drdx[1,1]/rhomean ; uuh*=drdx[2,2]/rhomean ; uup*=drdx[3,3]/rhomean
    ud0/=drdx[0,0]*rhomean ; udr/=drdx[1,1]*rhomean ; udh/=drdx[2,2]*rhomean ; udp/=drdx[3,3]*rhomean
    puu0*=drdx[0,0]/pmean ; puur*=drdx[1,1]/pmean ; puuh*=drdx[2,2]/pmean ; puup*=drdx[3,3]/pmean
    pud0/=drdx[0,0]*pmean ; pudr/=drdx[1,1]*pmean ; pudh/=drdx[2,2]*pmean ; pudp/=drdx[3,3]*pmean
    mpuu0*=drdx[0,0]/magp_mean ; mpuur*=drdx[1,1]/magp_mean ; mpuuh*=drdx[2,2]/magp_mean ; mpuup*=drdx[3,3]/magp_mean
    mpud0/=drdx[0,0]*magp_mean ; mpudr/=drdx[1,1]*magp_mean ; mpudh/=drdx[2,2]*magp_mean ; mpudp/=drdx[3,3]*magp_mean

#    uu0/=unorm ; uur/=unorm ; uuh/=unorm ; uup/=unorm
#    ud0/=unorm ; udr/=unorm ; udh/=unorm ; udp/=unorm
    # averaging the density:
    rhomean=rhomean/double(nframes)
    rhodisp=rhosqmean/double(nframes)-rhomean**2
    pmean=pmean/double(nframes)
    magp_mean=magp_mean/double(nframes)
    tudem/=double(nframes) ; tudma/=double(nframes)
    aphi/=double(nframes)
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
	fout.write(str(ph[0,0,kz])+'\n')
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
    
    rg("gdump")

    # velocities:
    uufile='merge_uu.dat'
    udfile='merge_ud.dat'
    fuu=open(uufile, 'r')
    fud=open(udfile, 'r')
#    s=str.split(str.strip(fuu.readline()))

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
        fname=dumpname(k)
        rd(fname)
        if(k==n1):
            rhomean=rho
            # velocity components:
            uu0=uu[0]*drdx[0,0]-uumean[0] ; ud0=ud[0]/drdx[0,0]-udmean[0]
            uur=uu[1]*drdx[1,1]-uumean[1] ; udr=ud[1]/drdx[1,1]-udmean[1]
            uuh=uu[2]*drdx[2,2]-uumean[2] ; udh=ud[2]/drdx[2,2]-udmean[2]
            uup=uu[3]*drdx[3,3]-uumean[3] ; udp=ud[3]/drdx[3,3]-udmean[3]
            tuur=(uu0*tr[0]+uur*tr[1]+uuh*tr[2]+uup*tr[3])
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
    rd(fname)

    run=unique(r)
    hun=unique(h)
    phun=unique(ph)

    indisk=double(abs(cos(h))<htor)
    innorm=indisk.mean(axis=1)

    ugmean=(ug*indisk).mean(axis=1)/innorm
    uu[0]*=drdx[0,0]  ; uu[1]*=drdx[1,1]  ; uu[2]*=drdx[2,2]  ; uu[3]*=drdx[3,3]  # to physical
    ud[0]/=drdx[0,0]  ; ud[1]/=drdx[1,1]  ; ud[2]/=drdx[2,2]  ; ud[3]/=drdx[3,3]  # to physical

    uumean0=(uu[0]*indisk*rho).mean(axis=1)/(rho*indisk).mean(axis=1)
    uumeanr=(uu[1]*indisk*rho).mean(axis=1)/(rho*indisk).mean(axis=1)
    uumeanh=(uu[2]*indisk*rho).mean(axis=1)/(rho*indisk).mean(axis=1)
    uumeanp=(uu[3]*indisk*rho).mean(axis=1)/(rho*indisk).mean(axis=1)
    uumean0_p=(uu[0]*indisk*ug).mean(axis=1)/(ug*indisk).mean(axis=1)
    uumeanr_p=(uu[1]*indisk*ug).mean(axis=1)/(ug*indisk).mean(axis=1)
    uumeanh_p=(uu[2]*indisk*ug).mean(axis=1)/(ug*indisk).mean(axis=1)
    uumeanp_p=(uu[3]*indisk*ug).mean(axis=1)/(ug*indisk).mean(axis=1)
    udmean0=(ud[0]*indisk*rho).mean(axis=1)/(rho*indisk).mean(axis=1)
    udmeanr=(ud[1]*indisk*rho).mean(axis=1)/(rho*indisk).mean(axis=1)
    udmeanh=(ud[2]*indisk*rho).mean(axis=1)/(rho*indisk).mean(axis=1)
    udmeanp=(ud[3]*indisk*rho).mean(axis=1)/(rho*indisk).mean(axis=1)
    udmean0_p=(ud[0]*indisk*ug).mean(axis=1)/(ug*indisk).mean(axis=1)
    udmeanr_p=(ud[1]*indisk*ug).mean(axis=1)/(ug*indisk).mean(axis=1)
    udmeanh_p=(ud[2]*indisk*ug).mean(axis=1)/(ug*indisk).mean(axis=1)
    udmeanp_p=(ud[3]*indisk*ug).mean(axis=1)/(ug*indisk).mean(axis=1)

    rhomean=(rho*indisk).mean(axis=1)/innorm
    pmag=(bsq*indisk).mean(axis=1)/innorm/2.

    # r mesh:
    fout=open(fname+'_eq_r.dat', 'w')
    for kx in arange(nx):
        if(kx%alifactor==0):
            fout.write(str(run[kx])+'\n')
    fout.close()
    fout=open(fname+'_eq_phi.dat', 'w')
    # phi mesh
    for kz in arange(nz):
	#        if(kz%alifactor==0):
        fout.write(str(phun[kz])+'\n')
    fout.close()

    foutrho=open(fname+'_eq_rho.dat', 'w') # density
    foutp=open(fname+'_eq_p.dat', 'w') # gas pressure
    foutpm=open(fname+'_eq_pm.dat', 'w') # magnetic pressure
    foutu=open(fname+'_eq_uu.dat', 'w') # contravariant velocities
    foutd=open(fname+'_eq_ud.dat', 'w') # covariant velocities
    foutup=open(fname+'_eq_puu.dat', 'w') # pressure-averaged contravariant velocities
    foutdp=open(fname+'_eq_pud.dat', 'w') # pressure-averaged covariant velocities

    # rho on the mesh
    for kx in arange(nx):
        if(kx%alifactor==0):
            for kz in arange(nz):
                foutrho.write(str(rhomean[kx,kz])+'\n')
                foutp.write(str(ugmean[kx,kz]*(gam-1.))+'\n')
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
    os.system('tar -cf '+fname+'_eq.tar '+fname+'_dinfo.dat '+fname+'_eq_*.dat')

# reading one frame and making it lighter and 2D and saving as an ascii
def framerip(fname, alifactor=1):
    global rho, ug, uu, ud, gam, B, aphi

    # alifactor = alias factor, every alifactorth point in r and th will be outputted, phi averaged over
#    rg("gdump")

    rd(fname)
	#    print(rho)
    run=unique(r)
    hun=unique(h)

    ug=ug.mean(axis=2)
#    uu[0]*=drdx[0,0]  ; uu[1]*=drdx[1,1]  ; uu[2]*=drdx[2,2]  ; uu[3]*=drdx[3,3]  # to physical
    uu0=(uu[0]*rho*drdx[0,0]).mean(axis=2)/rho.mean(axis=2)
    uur=(uu[1]*rho*drdx[1,1]).mean(axis=2)/rho.mean(axis=2)
    uuh=(uu[2]*rho*drdx[2,2]).mean(axis=2)/rho.mean(axis=2)
    uup=(uu[3]*rho*drdx[3,3]).mean(axis=2)/rho.mean(axis=2)

#   uu[0]*=drdx[0,0]  ; uu[1]*=drdx[1,1]  ; uu[2]*=drdx[2,2]  ; uu[3]*=drdx[3,3]  # to physical
#    ud[0]/=drdx[0,0]  ; ud[1]/=drdx[1,1]  ; ud[2]/=drdx[2,2]  ; ud[3]/=drdx[3,3]  # to physical

#    ud=(ud*rho).mean(axis=3)/rho.mean(axis=2)
#    ud[0]/=drdx[0,0]  ; ud[1]/=drdx[1,1]  ; ud[2]/=drdx[2,2]  ; ud[3]/=drdx[3,3]  # to physical
    ud0=(ud[0]*rho/drdx[0,0]).mean(axis=2)/rho.mean(axis=2)
    udr=(ud[1]*rho/drdx[1,1]).mean(axis=2)/rho.mean(axis=2)
    udh=(ud[2]*rho/drdx[2,2]).mean(axis=2)/rho.mean(axis=2)
    udp=(ud[3]*rho/drdx[3,3]).mean(axis=2)/rho.mean(axis=2)

    rhom=rho.mean(axis=2)
    pmag=(bsq).mean(axis=2)/2.
    aphi=psicalc()
    B[1]=B[1]*drdx[1,1] ; B[2]=B[2]*drdx[2,2] ; B[3]=B[3]*drdx[3,3]
    B=B.mean(axis=3)
#    aphi=aphi.mean(axis=2)

    # r mesh:
    fout=open(fname+'_r.dat', 'w')
    for kx in arange(nx):
        if(kx%alifactor==0):
            fout.write(str(run[kx])+'\n')
    fout.close()
    fout=open(fname+'_h.dat', 'w')
    # theta mesh
    for ky in arange(ny):
        if(ky%alifactor==0):
            fout.write(str(hun[ky])+'\n')
    fout.close()
    
    foutrho=open(fname+'_rho.dat', 'w')
    foutp=open(fname+'_p.dat', 'w')
    foutpm=open(fname+'_pm.dat', 'w')
    foutu=open(fname+'_uu.dat', 'w')
    foutd=open(fname+'_ud.dat', 'w')
    foutb=open(fname+'_b.dat', 'w')
    # rho on the mesh
    for kx in arange(nx):
        if(kx%alifactor==0):
            for ky in arange(ny):
                if(ky%alifactor==0):
                    foutrho.write(str(rhom[kx,ky])+'\n')
                    foutp.write(str(ug[kx,ky]*(gam-1.))+'\n')
                    foutpm.write(str(pmag[kx,ky]/2.)+'\n')
                    foutu.write(str(uu0[kx,ky])+' '+str(uur[kx,ky])+' '+str(uuh[kx,ky])+' '+str(uup[kx,ky])+'\n')
                    foutd.write(str(ud0[kx,ky])+' '+str(udr[kx,ky])+' '+str(udh[kx,ky])+' '+str(udp[kx,ky])+'\n')
		    foutb.write(str(B[1, kx,ky])+' '+str(B[2, kx,ky])+' '+str(B[3, kx,ky])+' '+str(aphi[kx,ky])+'\n')
    foutb.close()
    foutrho.close()
    foutp.close()
    foutpm.close()
    foutu.close()
    foutd.close()
    os.system('tar -cf '+fname+'.tar '+fname+'_*.dat')

def defaultrun():

    dire=''
    rref=10.
 #   os.chdir(dire)
#    dumpinfo()
    rg("gdump")
    run=unique(r)
    nr=run[where(run<rref)].argmax()

    fout=open("dumps_mevol.dat", "w")

    nmin=0
    nmax=70
#    glevol(nmax,5.)
#    glevol(nmax,10.)
    for k in arange(nmax-nmin+1)+nmin:
	dumpinfo(dumpname(k))
	framerip(dumpname(k), alifactor=3)
	maccre, mwind, laccre, lwind = mint(rref)
	fromabove(dumpname(k), alifactor=3)
	fout.write(str(t)+" "+str(maccre)+" "+str(mwind)+" "+str(laccre)+" "+str(lwind)+"\n")

    fout.close()

def corveerun(details=False):
    rg("gdump")

    if(details):
# transition regime
    	readndump(800,1000, pflag=True)
    	os.system('cp merge_r.dat mergeA_r.dat')
    	os.system('cp merge_h.dat mergeA_h.dat')
    	os.system('cp merge_uu.dat mergeA_puu.dat')
    	os.system('cp merge_ud.dat mergeA_pud.dat')
    	os.system('cp merge_rho.dat mergeA_rho.dat')
    	os.system('cp merge_p.dat mergeA_p.dat')
    	readndump(800,1000)
    	os.system('cp merge_uu.dat mergeA_uu.dat')
    	os.system('cp merge_ud.dat mergeA_ud.dat')
# maximum
    	readndump(1000,1175, pflag=True)
    	os.system('cp merge_r.dat mergeB_r.dat')
    	os.system('cp merge_h.dat mergeB_h.dat')
    	os.system('cp merge_rho.dat mergeB_rho.dat')
    	os.system('cp merge_p.dat mergeB_p.dat')
    	os.system('cp merge_uu.dat mergeB_puu.dat')
    	os.system('cp merge_ud.dat mergeB_pud.dat')
    	readndump(1000,1175)
    	os.system('cp merge_uu.dat mergeB_uu.dat')
   	os.system('cp merge_ud.dat mergeB_ud.dat')
# low mdot
    	readndump(1180,1280, pflag=True)
    	os.system('cp merge_r.dat mergeC_r.dat')
    	os.system('cp merge_h.dat mergeC_h.dat')
    	os.system('cp merge_rho.dat mergeC_rho.dat')
    	os.system('cp merge_p.dat mergeC_p.dat')
    	os.system('cp merge_uu.dat mergeC_puu.dat')
    	os.system('cp merge_ud.dat mergeC_pud.dat')
    	readndump(1180,1280)
    	os.system('cp merge_uu.dat mergeC_uu.dat')
   	os.system('cp merge_ud.dat mergeC_ud.dat')
# maximum
    	readndump(1300,1600, pflag=True)
    	os.system('cp merge_r.dat mergeD_r.dat')
    	os.system('cp merge_h.dat mergeD_h.dat')
    	os.system('cp merge_rho.dat mergeD_rho.dat')
    	os.system('cp merge_p.dat mergeD_p.dat')
    	os.system('cp merge_uu.dat mergeD_puu.dat')
    	os.system('cp merge_ud.dat mergeD_pud.dat')
    	readndump(1300,1600)
    	os.system('cp merge_uu.dat mergeD_uu.dat')
    	os.system('cp merge_ud.dat mergeD_ud.dat')
# total quasi-steady accretion phase
    readndump(1000,1905)
#    corvee(1000, 1905)
#    os.system('cp merge_uu.dat merge_puu.dat')
#    os.system('cp merge_ud.dat merge_pud.dat')
#    os.system('cp merge_u3d.dat merge_pu3d.dat')
#    readndump(1000,1999, pflag=True,withmag=True)
#    os.system('cp merge_uu.dat merge_mpuu.dat')
#    os.system('cp merge_ud.dat merge_mpud.dat')
#    os.system('cp merge_u3d.dat merge_mpu3d.dat')
    os.system('tar -cf mergereads.tar merge_*.dat')

#    readndump(1000,1905)

#    corvee(200,775) 
#    os.chdir('..')

# rg('gdump')
# print ph

defaultrun()
# corveerun()

