# genray_plots.py'
# Plots genray.nc (output data file produced by GENRAY)

# Yuri Petrov, Bob Harvey   CompX   2013

# Needed in working directory: genray.nc and eqdsk (if available)
# Rename your particular equilibrium-B file to eqdsk.

from numpy import *
from mpl_toolkits.mplot3d import Axes3D

from pylab import *
from matplotlib import rc 
from matplotlib.pyplot import cm,figure,axes,plot,xlabel,ylabel,title,savefig,show

import os
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import time
import pylab as pylab
import scipy.io.netcdf as nc

#-----------------------------------------------
# NetCDF issues: machine-dependent
# Try netcdf=4 to envoke netCDF4,
# Or try netcdf=2 to work with older netCDF.
netcdf=4
#-----------------------------------------------
if netcdf==4: from netCDF4 import Dataset # YuP
#-----------------------------------------------

#------ For Linux users: un-comment  matplotlib.rc('text', usetex = True)
#matplotlib.interactive(True) # no plots on screen
matplotlib.interactive(False) # with plots on screen
#Render with externally installed LateX (BH200816):
matplotlib.rc('text', usetex = True) # Does not work on PC version


#-------- Specify for plots:
fnt  = 11 #13 #9   # Font size for axis numbers (see 'param=' below) 
linw = 1.  # LineWidth for contour plots
Ncont= 50  # Number of contour levels for PSI (pol.flux)
# For plots of omega/omega_c in (R,Z) plane ("genray_wwc_inRZ.png"),
# specify:
level_mn= 1 #23 # Lowest level of omega/omega_c to be plotted
level_mx= 9 #45 # Highest level of omega/omega_c to be plotted
arr_len=5. # Specify arrow length for plots of (Nr,Nz) at start point (disabled)
nskip=1  # How many end-points to skip, for plotting (to avoid end-point jumps)
rays_style=0 # 0 Rays in (R,Z) are plotted with constant line thickness.
             # 1 Rays in (R,Z) are plotted with variable width proportional to 
             #    remaining power in ray channel (may become too thin/invisible)
             # 2 Rays are plotted with markers which size corresponds to
             #    absorbed power.
isp_wc=1  #  Species number for which you want to plot omega/omega_c res.layers
          #  1- electrons, 2 or higher - ions
             
# Constants
pi=3.14159265358979
clight= 2.99792458e10   # speed of light [cm/s]
e     = 4.8032e-10      # e-charge [cgs]
e_mass= 9.10938291e-28  # e-mass   [gram]
p_mass= 1.67262158e-24  # proton mass    [gram]
ergtkev=1.6022e-09
mp=p_mass/e_mass

# Make a guess about Min/Max of plasma; 
# Set to 0 for automatic setting
# (In this case Rmin,Rmax,Zmin,Zmax will be found from eqdsk, if present
#  and will be used for plots)
#Rmin_plot=5 #120 # 400.
#Rmax_plot=90 #240 # 850. # just any guess [cm] 
#Zmin_plot=-70 #-100 #-400.
#Zmax_plot=+70 #+100 #500. # just any guess [cm] 
Rmin_plot=0.
Rmax_plot=0. 
Zmin_plot=0.
Zmax_plot=0.

# Additional structures (ports, coils) to be plotted in (R,Z) section
istruct=0 # To skip, set istruct=0
Rstruct1= np.array([156.,156.,175.,175.,156.])  # PF coil 1
Zstruct1= np.array([ 26., 45., 45., 26., 26.])  # PF coil 1
Rstruct2= np.array([156.,156.,175.,175.,156.])  # PF coil 2
Zstruct2= np.array([-27.,-46.,-46.,-27.,-27.])  # PF coil 2
# Window ports:
Rstruct3= np.array([184.,184.,195.,195.,184.])
Zstruct3= np.array([44.,114.,114.,44.,44.])
Rstruct4= np.array([184.,184.,195.,195.,184.])
Zstruct4= np.array([-25.,25.,25.,-25.,-25.])
#should be 4 structures only; some will be plotted with mirror image,too

e0 = time.time()  # elapsed time since the epoch
c0 = time.clock() # total cpu time spent in the script so far

#------------------------------------------
#eqdsk_name='eqdsk'
#eqdsk_name='g1100901026.01060'  # LHCD_Sugiyama  test
#eqdsk_name='eqdsk_NSTX'
#eqdsk_name='eqdsk_ITER_sym'
#eqdsk_name='eqdsk_H24'

#eqdsk_name='140358R01_350ms'  # for PPPL/IPS/Poli test
#eqdsk_name='12028Z33_ps.geqdsk' # PPPL Jin Chen

#eqdsk_name='g010004.01020_21t' # see /GENRAY_wk/test_ysbae_problem/

#eqdsk_name='G_STEADY_EC1F_ITER_LR_01000.TXT'
#eqdsk_name='G_HYBRID_LH03_ITER_LR_00400.TXT'

#eqdsk_name='g7777.002001'  # KAERI/S.Ho Kim  LH case 2017 07/19
#eqdsk_name='g032974.02300'  # ST QUEST 07/20 2017 f=28e9

#eqdsk_name='Scen4_bn2.57_129x129' # /test1/
#eqdsk_name='equilib.dat_HPRT_test_case_061219' # /test2/
#eqdsk_name='g113544.00325_mod' # /test3/
#eqdsk_name='g113544.00325_mod' # /test4/  (OXB, i_ox=2)
#eqdsk_name='g521022.01000'     # /test5/  (canonical 2004 ITER 1ray)
#eqdsk_name='g106270.02500'      # /test6/ multi-rays EC
#eqdsk_name='g1060728011.01100'  # /test7/ (LHCD)
#eqdsk_name='equilib.dat_EC_ITER_Central_CD' #[rename to equilib.dat] /test8/
#eqdsk_name='g1060728011.01100'  # /test9/
#eqdsk_name='scenario_c1_gc_irod_110_rout_146_cm' # test 10.3

#eqdsk_name='eqdsk_pegasus'
#eqdsk_name='equilib_diiid.dat' # for Lohr/fast-run-test
#eqdsk_name='g130608.00352.LRDFIT04.mds.corrected.qscale_1.00000' #test_NSTX_HHFW
#eqdsk_name='eqdsk_240rays'

#eqdsk_name='eqdsk_EAST'
#eqdsk_name='Equilibriumfile_3e19_B_1.5T_Te_1.7keV.dat' #see /issue_ADITYA/
#eqdsk_name='Equilibriumfile_7e19_B_1.5T_Te_3.4keV.dat'

#----eqdsk_name='Pegasus_eqdsk2' #'Pegasus_eqdsk'

#----- UKAEA MAST runs --------------------------------
#eqdsk_name='C1a.geqdsk'   # 2020-07-22   case 110_110 (unscaled B field)
# Rescaled Btor and Bpol (see /EQDSK/) 
#eqdsk_name='C1a.geqdsk_80_110' # Both are rescaled by 80/110.
#eqdsk_name='C1a.geqdsk_120_110' # Both are rescaled by 120/110.
#  2020-09-04  sensitivity scans 
#  Use files like scenario_c1_gc_irod_110_rout_146_cm (no extension)
#  which are adjusted from scenario_c1_gc_irod_110_rout_146_cm.eqdsk
#eqdsk_name='scenario_c1_gc_irod_110_rout_146_cm' 
#eqdsk_name='scenario_c1_gc_irod_110_rout_143_cm' 
#eqdsk_name='scenario_c1_gc_irod_110_rout_140_cm' 
#eqdsk_name='scenario_c1_gc_irod_110_rout_137_cm' 
#eqdsk_name='scenario_c1_gc_irod_110_rout_134_cm' 
#eqdsk_name='scenario_c1_gc_irod_110_rout_132_cm' 
#eqdsk_name='scenario_c1_gc_irod_110_rout_130_cm' 
#eqdsk_name='scenario_c1_gc_irod_110_rout_128_cm' 
#eqdsk_name='scenario_c1_gc_irod_110_rout_126_cm'

eqdsk_name='g_000000.00000_anv1'  # issue_OXB_high_freq


#------------------------------------------
# Open the genray netcdf file; specify name:
filenm='genray.nc' # This is a common name of the file generated by GENRAY 
#filenm='genray_iter5_060828.1.nc' # LHCD in ITER (also used for alpha heating)
#filenm='iter5_060828.nc'  # LH

#filenm='iter_ec1f_helicon.nc'           # 500MHz z=Zax
####filenm='iter_ec1f_helicon500shift1m.nc' # 500MHz shifted 1m down
#filenm='iter_ec1f_helicon800.nc'        # 800MHz z=Zax
#filenm='iter_ec1f_helicon800shift1m.nc' # 800MHz shifted 1m down

#filenm='iter_lh60s_helicon500.nc'        # 500MHz no shift
#filenm='iter_lh60s_helicon500shift1m.nc' #500MHz with antenna vert. shift 1m down
#filenm='iter_lh60s_helicon800.nc'        # 800MHz no shift
#filenm='iter_lh60s_helicon800shift1m.nc' #800MHz with antenna vert. shift 1m down

#filenm='genray_kstar_lhsw.nc'  # KAERI/S.Ho Kim  LH case 2017 07/19  

#filenm='Pegasus_rayech.nc'
#filenm='Pegasus_genray.nc'




#==================================================================
def read_vector(flnm,ndim,nums_per_line):
#     global nlines  #Can be used to return number of lines.
    """
    Reads delimited items from an open file, flnm.
    Returns them in a vector of input length ndim.
    nums_per_line is (assumed constant) number of
    items per line (except last line may be shorter).
    BH2009.
    """
    a=np.ones(ndim)
    nlines=int(ndim/nums_per_line)
    if nlines*nums_per_line != ndim: nlines=nlines+1
    #print nlines
    for i in range(nlines):
        ibegin=i*nums_per_line
        iend=min(ibegin+nums_per_line,ndim)
        #print ibegin,iend
        a[ibegin:iend]=np.array(flnm.readline().split(),float)
    return a
#==================================================================
# From https://stackoverflow.com/questions/19390895/
#              matplotlib-plot-with-variable-line-width
def plot_widths(xs, ys, widths, ax=None, color='b', xlim=None, ylim=None,
                **kwargs):
    if not (len(xs) == len(ys) == len(widths)):
        raise ValueError('xs, ys, and widths must have identical lengths')
    #fig = None
    #if ax is None:
    #    fig, ax = plt.subplots(1)

    segmentx, segmenty = [xs[0]], [ys[0]]
    current_width = widths[0]
    for ii, (x, y, width) in enumerate(zip(xs, ys, widths)):
        segmentx.append(x)
        segmenty.append(y)
        if (width != current_width) or (ii == (len(xs) - 1)):
            plot(segmentx, segmenty, linewidth=current_width, color=color,
                    **kwargs)
            segmentx, segmenty = [x], [y]
            current_width = width
    #if xlim is None:
    #    xlim = [min(xs), max(xs)]
    #if ylim is None:
    #    ylim = [min(ys), max(ys)]
    #ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    #return ax if fig is None else fig
        
        
#******************************************************************
#******************************************************************
i_eqdsk=0  # it will be set to 1 if eqdsk is found
nlimiter=0 # it will be reset to other value, if data on limiter is avail. 
# Min/Max of plasma:
Rmin= 0.0 #initialize [cm]
Rmax= 0.0 #initialize [cm]
Zmin= 0.0 #initialize [cm] 
Zmax= 0.0 #initialize [cm]

try:# Try reading eqdsk file, hoping it exists in the working directory.
    # Copy equilibrium B data file into working directory under name 'eqdsk'
    eqdsk=open(eqdsk_name,'r')
    print('reading eqdsk file:', eqdsk_name)
    lines=eqdsk.readlines()
    #fortran format 110  format(6a8,4i4), picking off 2nd and 3rd slots of 4i4.
    nr=lines[0][52:56]; nz=lines[0][56:60]
    nr=int(nr); nz=int(nz)
    try:
        nv=int(lines[0][60:64])
    except:   #line[60:64] not characters representing integers
        nv=nr
    #A difficulty making reading of eqdsk files difficult is that
    #when a number is procedded by a "-" sign, then there is no
    #blank space separating adjacent numbers.
    #For remaining lines after 1st down to possible occurance of
    #an "&", fortran format is format(5e16.9) [Except possible 
    #sequence of integers extending less that 1st 16 columns.
    #Therefore, Add blanks in columns 17,33,49,65, and use
    #above read_vector() function.
    print('nr=',nr, '  nz=',nz)
    #Determine line number of "&", if exists:
    #(lines[iamp] will be line before "&".)
    iamp=0
    for line in lines[1:]:
        if line.strip().find('&') == -1:
            iamp=iamp+1
        else:
            break

    #Introduce spaces in lines[1:iamp], skipping 1st and & lines.
    #This should be ok for standard eqdsks through reading xlimiter,ylimiter.
    #BH noticed a nonstandard format (i5,e16.9) further on.
    cols=list(range(16,65,16)); cols.reverse()  #cols=[64, 48, 32, 16
    b=' '
    for i in range(1,iamp+1):
        for j in range(4):
            lines[i]=lines[i][0:cols[j]]+b+lines[i][cols[j]:]

    # Open tmp file, put adjusted lines in it, 
    # then read using previously constructed read_vector.
    tmp=open('tmpfile','w')
    tmp.writelines(lines)
    tmp.close()
    tmp=open('tmpfile','r')
    tmp.readline()  #space down one line
    #---
    rbox,zbox,radmaj,rboxdst,ymideqd=read_vector(tmp,5,5)
    print('rbox=',rbox)
    raxis,zaxis,psimag,psilim,btor=read_vector(tmp,5,5)
    print('raxis=',raxis)
    toteqd,psimx1,psimx2,xax1,xax2=read_vector(tmp,5,5)
    zax1,zax2,psisep,xsep,zsep=read_vector(tmp,5,5)
    fpsiar=read_vector(tmp,nv,5)
    prar=read_vector(tmp,nv,5)
    ffpar=read_vector(tmp,nv,5)
    ppar=read_vector(tmp,nv,5)
    epsi=read_vector(tmp,nr*nz,5)
    epsi.resize((nz,nr))
    qar=read_vector(tmp,nv,5)
    ncontr, nlimiter=read_vector(tmp,2,2)
    print('ncontr, nlimiter=', ncontr, nlimiter)
    ncontr=int(ncontr); nlimiter=int(nlimiter)
    #---
    # raxis,zaxis      are the major and vertical height of magnetic axis.
    # psimag,psilim    are the poloidal flux function values at the
    #                  magnetic axis and the last closed flux surface
    #                   (touching the limiter or the separatrix).
    # btor             is the vacuum toroidal magnetic field at radmaj.
    # toteqd           is the toroidal current.
    #
    # feqd(nnr),pres(nnr), usually, although dimension is
    #                      specified by nnv if it is present.
    # fpsiar = r * B_phi  at nnr equispaced points in psi
    #          from psimag to psilim.
    # prar   = are the pressure values at the same points.
    # ffpar  = fpsiar_prime, [i.e., d f/d (psi)] at the same points.
    # ppeqd  =  p_prime [i.e., d prar/d (psi)] at the same points.
    # epsi(nnr,nnz) are the psi values on the nnr * nnzprint 'nv=',nv
    #               equispaced grid.
    # qar(nnr)  gives safety factor q on the equispaced psi grid.
    #
    # The following quantities are given in some eqdsks, but are not
    # presently used in cql3d:
    # nlimit,nves  are the numbers of point at limiters and
    #              vacuum vessel wall.
    # rlimit(nlimit),zlimit(nlimit):
    # rlimit,zlimit      is the r,z location of the limiters wall.
    #       rves(nves),zves(nves):
    # rves,zves          is the r,z location of the limiters
    #                    vacuum vessel wall.
    print('========================================')
    print('nr=',nr)
    print('nz=',nz)
    print('nv=',nv)
    print('psimag=', psimag)
    print('psilim=', psilim)
    print('zaxis [m]=', zaxis)
    print('raxis [m]=', raxis)
    print('radmaj [m]=', radmaj)
    print('btor [Tesla]=',btor)
    print('fpsiar shape =', fpsiar.shape)
    print('epsi shape =', epsi.shape)
    print('nlimiter=',nlimiter)
    print('ncontr=',ncontr)
    print('read_eqdsk done')
    print('========================================')
    # convert to cgs:
    rbox=rbox*1.e+2 
    zbox=zbox*1.e+2
    rboxdst=rboxdst*1.e+2
    radmaj=radmaj*1.e+2
    toteqd=toteqd*3.e+9       
    R_axis  = raxis*1.e+2 # cgs
    Z_axis  = zaxis*1.e+2 # cgs
    fpsiar=fpsiar*1.e+6   # r * B_phi  [cgs now]
    epsi=  epsi*1.e+8     # cgs now
    psilim= psilim*1.e8  
    psimag= psimag*1.e8  
    btor=btor*1e4 # vacuum toroidal magnetic field at radmaj [Gauss]
    # Min/Max of plasma:
    #Rmin= rboxdst
    #Rmax= Rmin + rbox
    #Zmin=(zaxis -zbox*.5) # in [cm]
    #Zmax= Zmin + zbox     # in [cm]
    # Form R,Z grids
    dzz=zbox/(nz-1) # cgs
    drr=rbox/(nr-1) # cgs
    er=np.zeros(nr)
    ez=np.zeros(nz)
    er[0]=rboxdst    # cgs
    ez[0]= zaxis -zbox*.5
    for nn in range(1,nr,1):  # nn goes from 1 to nr-1
        er[nn]=er[nn-1]+drr   # R-grid [cm] 
    for nn in range(1,nz,1):  # nn goes from 1 to nz-1
        ez[nn]=ez[nn-1]+dzz   # Z-grid [cm]
    R,Z = np.meshgrid(er,ez) # 2D grids [cm]    
    #........................................................
    # Form the equally spaced psi array/grid.
    # psimag < psilim;  epsi has min. at m.axis
    delpsi= (psilim-psimag)/(nv-1)  
    psiar=np.zeros(nv)
    for ix in range(0,nv,1):  # ix goes from 0 to nv-1
        psiar[ix]= psimag + ix*delpsi # [psimag; psilim]
    #.........................................................
    # Grad(psi) components:
    grpsi_r, grpsi_z = np.gradient(epsi,dzz,drr)
    grpsi= sqrt(grpsi_r*grpsi_r + grpsi_z*grpsi_z)
    Bpol= zeros((nz,nr))
    Btor= zeros((nz,nr))
    Bpmn=1.e10 # will be found below
    Bpmx=0.    # will be found below
    # Define |Bpol|  and  Btor  [Gauss]
    for ir in range(0,nr,1):  # ir goes from 0 to nr-1
        for iz in range(0,nz,1):  # iz goes from 0 to nz-1
            Bpol[iz,ir]= grpsi[iz,ir]/er[ir]
            if psilim>psimag:
                if epsi[iz,ir]<psilim:
                    iv= int((epsi[iz,ir]-psimag)/delpsi)
                    #if iz==64: print epsi[iz,ir], iv, psiar[iv],psiar[iv+1]
                    Btor[iz,ir]=fpsiar[iv]/er[ir]
                    Bpmn=min(Bpol[iz,ir],Bpmn)
                    Bpmx=max(Bpol[iz,ir],Bpmx)
                else:
                    Btor[iz,ir]=btor*radmaj/er[ir]
            else: # psilim<psimag
                if epsi[iz,ir]>psilim:
                    iv= int((epsi[iz,ir]-psimag)/delpsi)
                    #if iz==64: print epsi[iz,ir], iv, psiar[iv],psiar[iv+1]
                    Btor[iz,ir]=fpsiar[iv]/er[ir]
                    Bpmn=min(Bpol[iz,ir],Bpmn)
                    Bpmx=max(Bpol[iz,ir],Bpmx)
                else:
                    Btor[iz,ir]=btor*radmaj/er[ir]            
    B=sqrt(Btor*Btor + Bpol*Bpol) # Total B as a function of (ir,iz)
    #print"shape(B)=" , shape(B)
    if ncontr>4:
        rzcontr=read_vector(tmp,2*ncontr,5)
        rzcontr.resize((ncontr,2))
        rcontr=rzcontr[:,0]
        zcontr=rzcontr[:,1]
        print 'min/max of rcontr=',np.min(rcontr),np.max(rcontr)
        print 'min/max of zcontr=',np.min(zcontr),np.max(zcontr)        
    if nlimiter>4:
        try:
            rzlimiter=read_vector(tmp,2*nlimiter,5)
            #print('rzlimiter=',rzlimiter)
            rzlimiter.resize((nlimiter,2))
            rlimiter=rzlimiter[:,0]
            zlimiter=rzlimiter[:,1]
            print 'min/max of rlimiter=',np.min(rlimiter),np.max(rlimiter)
            print 'min/max of zlimiter=',np.min(zlimiter),np.max(zlimiter)        
        except:  
            print('possibly corrupted data on limiter')
            nlimiter=0
        
    eqdsk.close()  # close 'eqdsk' file
    tmp.close()
    i_eqdsk=1
    
except IOError:
    print("i_eqdsk=",i_eqdsk)
    
if i_eqdsk==0: # eqdsk file is not readable
    nlimiter=0
    ncontr=0
    print("===========================================================")
    print("eqdsk file is not in the working directory or not readable.")
    print("===========================================================")

#******************************************************************
#******************************************************************







#Input netcdf file into a structure:
if netcdf==2: 
    dat=nc.netcdf_file(filenm,'r')
 
#------------YuP:
if netcdf==4: 
    dat= Dataset(filenm, 'r', format='NETCDF4')

print(dat.file_format) # print which format was used for the genray.nc file



print('The genray file, ',filenm,', contains:')
print('========================================')
print("The global attributes: ",dat.dimensions.keys())    
print("File contains variables: ",dat.variables.keys())
print('========================================')

# The name of eqdsk file used for equilibrium B-field:
eqdskin=dat.variables['eqdskin']
print("eqdskin=",eqdskin[:])

n_eqdsk=0 # initialize
if i_eqdsk==0: # no data was read from eqdsk 
    # some eqdsk data:
    # Note: this data did not exist in early genray.nc:
    try:
        try:
            eqdsk_r=dat.variables['eqdsk_r']
        except:
            print('No data on eqdsk in genray.nc')
            n_eqdsk=0
        else:
            n_eqdsk=1
            eqdsk_r=dat.variables['eqdsk_r']
            print('eqdsk_r is: ', eqdsk_r.long_name, eqdsk_r.shape)
            eqdsk_z=dat.variables['eqdsk_z']
            print('eqdsk_z is: ', eqdsk_z.long_name, eqdsk_z.shape)
            eqdsk_psi=dat.variables['eqdsk_psi']
            print('eqdsk_psi is: ', eqdsk_psi.long_name, eqdsk_psi.shape)
    finally:
        print('----------------------------------------')


n_wall=0 # initialize

# Note: r_wall array may not exist in genray.nc:
try:
    try:
        r_wall=dat.variables['r_wall']
    except:
        print('No data on r_wall, z_wall in genray.nc')
        n_wall=0
    else:
        r_wall=dat.variables['r_wall']
        z_wall=dat.variables['z_wall']
        n_wall=r_wall.shape[0]
        print('r_wall is: ', r_wall.long_name, r_wall.shape)
        print('z_wall is: ', z_wall.long_name, z_wall.shape)
        print('n_wall=',n_wall)
finally:
    print('----------------------------------------')
    
    
# Note: powden_cl array did not exist in early genray.nc:
try:
    try:
        powden_cl=dat.variables['powden_cl']
    except:
        print('No data on powden_cl in genray.nc')
        n_powden_cl=0
    else:
        n_powden_cl=1
        powden_cl=dat.variables['powden_cl']
        print('powden_cl is: ',powden_cl.units,powden_cl.shape)
finally:
    print('----------------------------------------')



# O-X modes Transmission data: Only present for i_ox=2 runs:
i_ox=0 # initialize
try:
    try:
        i_ox_conversion=dat.variables['i_ox_conversion']
    except:
        print('No data on O-X conversion in genray.nc')
        i_ox=0
    else:
        i_ox=2
        i_ox_conversion=dat.variables['i_ox_conversion'] #
        print 'i_ox_conversion is:',i_ox_conversion.long_name
        transm_ox=dat.variables['transm_ox'] #
        print 'transm_ox is: ', transm_ox.long_name, transm_ox.shape
        cnpar_ox=dat.variables['cnpar_ox'] #
        print 'cnpar_ox is: ', cnpar_ox.long_name, cnpar_ox.shape
        cn_b_gradpsi=dat.variables['cn_b_gradpsi'] #
        print 'cn_b_gradpsi is: ', cn_b_gradpsi.long_name, cn_b_gradpsi.shape
        cn_par_optimal=dat.variables['cn_par_optimal'] #
        print 'cn_par_optimal is: ', cn_par_optimal.long_name
        print cn_par_optimal[:]
finally:
    print '----------------------------------------'


# Define min/max for major radius range:
R_LCFS_min=0.0 #initialize
R_LCFS_max=0.0 #initialize
Z_LCFS_min=0.0 #initialize
Z_LCFS_max=0.0 #initialize
Rlimiter_min=0.0 #initialize
Rlimiter_max=0.0 #initialize
Zlimiter_min=0.0 #initialize
Zlimiter_max=0.0 #initialize
R_wall_min=0.0
R_wall_max=0.0
Z_wall_min=0.0
Z_wall_max=0.0
if Rmax_plot==0:
    if n_wall >4:  # data on wall is present
        R_wall_max=np.amax(r_wall)*100
        R_wall_min=np.amin(r_wall)*100
        Z_wall_max=np.amax(z_wall)*100
        Z_wall_min=np.amin(z_wall)*100
        print('Data on wall: R_wall_min,R_wall_max=',R_wall_min,R_wall_max)
        Rmin=min(Rmin,R_wall_min)
        Rmax=max(Rmax,R_wall_max)
        Zmin=min(Zmin,Z_wall_min)
        Zmax=max(Zmax,Z_wall_max)
    if ncontr >4:  # data on LCFS is present
        R_LCFS_max=np.amax(rcontr)*100
        R_LCFS_min=np.amin(rcontr)*100
        Z_LCFS_max=np.amax(zcontr)*100
        Z_LCFS_min=np.amin(zcontr)*100
        print('Data on LCFS (ncontr): R_LCFS_min,R_LCFS_max=',R_LCFS_min,R_LCFS_max)
        Rmin=min(Rmin,R_LCFS_min)
        Rmax=max(Rmax,R_LCFS_max)
        Zmin=min(Zmin,Z_LCFS_min)
        Zmax=max(Zmax,Z_LCFS_max)
    if nlimiter >4: # data on limiter is present
        Rlimiter_max=np.amax(rlimiter)*100
        Rlimiter_min=np.amin(rlimiter)*100
        Zlimiter_max=np.amax(zlimiter)*100
        Zlimiter_min=np.amin(zlimiter)*100
        print('Data on limiter: Rlimiter_min,Rlimiter_max=',Rlimiter_min,Rlimiter_max)
        Rmin=min(Rmin,Rlimiter_min)
        Rmax=max(Rmax,Rlimiter_max)
        Zmin=min(Zmin,Zlimiter_min)
        Zmax=max(Zmax,Zlimiter_max)
    if istruct==1:
        Rstruct1_max=np.amax(Rstruct1) # already in cm
        Rstruct2_max=np.amax(Rstruct2)
        Rstruct3_max=np.amax(Rstruct3)
        Rstruct4_max=np.amax(Rstruct4)
        Rstruct_max= np.max([Rstruct1_max,Rstruct2_max,Rstruct3_max,Rstruct4_max])
        #Rmin=min(Rmin,Rstruct_min)
        Rmax=max(Rmax,Rstruct_max)
        #Zmin=min(Zmin,Zstruct_min)
        #Zmax=max(Zmax,Zstruct_max)
    if n_wall+ncontr+nlimiter<5:  #nothing is defined; get limits from eqdsk:
        if n_eqdsk==1:
            Rmax=np.amax(eqdsk_r)*100
            Rmin=np.amin(eqdsk_r)*100 # meters -> cm
            Zmax=np.amax(eqdsk_z)*100
            Zmin=np.amin(eqdsk_z)*100 # meters -> cm
            print 'getting limits from eqdsk_r; Rmin,Rmax[cm]=',Rmin,Rmax
        else: 
            Rmin= 50.
            Rmax= 300. # just any guess [cm] 
            Zmin=-200.
            Zmax= 200. # just any guess [cm] 

    
print('----------------------------------------')
if Rmin_plot==0:
    xmin=Rmin  # Limits for plots
else:
    xmin=Rmin_plot
if Rmax_plot==0:            
    xmax=Rmax  # 
else:
    xmax=Rmax_plot
if Zmin_plot==0:
    zmin=Zmin  # Limits for plots
else:
    zmin=Zmin_plot
if Zmax_plot==0:            
    zmax=Zmax  # 
else:
    zmax=Zmax_plot

Rmin_plot=xmin
Rmax_plot=xmax
Zmin_plot=zmin
Zmax_plot=zmax
print('Rmin_plot,Rmax_plot=', Rmin_plot,Rmax_plot)
print('Zmin_plot,Zmax_plot=', Zmin_plot,Zmax_plot)
print('----------------------------------------')

freqcy=dat.variables['freqcy']
#BH  print 'freqcy =',freqcy.long_name, freqcy[:], freqcy.units
print('freqcy =',freqcy.long_name, freqcy.getValue(), freqcy.units)
#BH f=freqcy[0]
f=freqcy.getValue()

mass=dat.variables['dmas']
print('mass=', mass.long_name, mass[:], mass.units)

charge=dat.variables['charge']
print('charge=', charge.long_name, charge[:], charge.units)

rho_bin_center=dat.variables['rho_bin_center']
print('rho_bin_center is: ', rho_bin_center.long_name, rho_bin_center.shape)
rho_bin_center= np.asarray(rho_bin_center)
print('rho_bin_center=',rho_bin_center)

rho_bin=dat.variables['rho_bin']
print('rho_bin is: ', rho_bin.long_name, rho_bin.shape)

densprof=dat.variables['densprof']
print('densprof is: ', densprof.long_name, densprof.shape)

temprof=dat.variables['temprof']
print('temprof is: ', temprof.long_name, temprof.shape)

Nsp=temprof[:,0].size  # number of species
print('Number of species: Nsp=',Nsp)
# Identify species:
isp_name=[]
for isp in range(0,Nsp,1):
    msme= np.asscalar(mass[isp]/mass[0])      # ms/me ratio
    Zs= np.asscalar(charge[isp]/charge[0])    # Zs/e
    #isp_name[i]='0' # initialize
    #print 'isp=',isp,' mass=',mass[isp],' Z=',charge[isp],' isp_name=',isp_name
    if (np.absolute(Zs-1.)<0.1) & (np.absolute(msme-1.)<0.1):
        isp_name.append('e')
    if (np.absolute(Zs-1.)<0.1) & (np.absolute(msme-mp)<10):  # mp is 1836 
        isp_name.append('H')
    if (np.absolute(Zs-1.)<0.1) & (np.absolute(msme-2*mp)<20):
        isp_name.append('D')
    if (np.absolute(Zs-1.)<0.1) & (np.absolute(msme-3*mp)<20):
        isp_name.append('T')
    if (np.absolute(Zs-2.)<0.1) & (np.absolute(msme-3*mp)<100):
        isp_name.append('He3')
    if (np.absolute(Zs-2.)<0.1) & (np.absolute(msme-4*mp)<100):
        isp_name.append('He4')
    if (np.absolute(Zs-3.)<0.1) :
        isp_name.append('Li+3')
    if (np.absolute(Zs-4.)<0.1) & (np.absolute(msme-12*mp)<100):
        isp_name.append('C+4')
    if (np.absolute(Zs-5.)<0.1) & (np.absolute(msme-12*mp)<100):
        isp_name.append('C+5')
    if (np.absolute(Zs-6.)<0.1) & (np.absolute(msme-12*mp)<100):
        isp_name.append('C+6')
    # OTHER NAMES CAN BE ADDED HERE.
    print 'isp=',isp,' m=',mass[isp],' Z=',charge[isp],' isp_name=',isp_name[isp]
    #txt1=r"$m/m_e =$"+r"$%1.0f$" %(mass[i])
    #txt2=r"$q/e =$"+r"$%1.0f$" %(charge[i])
print isp_name
    


istart=dat.variables['istart']
istart=np.asscalar(istart.getValue())
print 'istart=',istart

iray_status_nc= dat.variables['iray_status_nc']
print('iray_status_nc.shape', iray_status_nc.shape)
print(iray_status_nc.long_name)
print('iray_status_nc=', iray_status_nc[:])
print(iray_status_nc.long_name1)
print(iray_status_nc.long_name2)
print(iray_status_nc.long_name3)
print(iray_status_nc.long_name4)
print(iray_status_nc.long_name5)
print(iray_status_nc.long_name6)
print(iray_status_nc.long_name7)
print(iray_status_nc.long_name8)
print(iray_status_nc.long_name9)
print(iray_status_nc.long_name10)
print(iray_status_nc.long_name11)
print(iray_status_nc.long_name12)
print(iray_status_nc.long_name13)
print(iray_status_nc.long_name14)
print(iray_status_nc.long_name15)
print(iray_status_nc.long_name16)
print(iray_status_nc.long_name17)
print(iray_status_nc.long_name18)
print(iray_status_nc.long_name19)
#print(iray_status_nc.long_name20)

# Averaged current densities
s_cur_den_parallel=dat.variables['s_cur_den_parallel']
print('s_cur_den_parallel is: ', s_cur_den_parallel.long_name, s_cur_den_parallel.shape)
s_cur_den_parallel=np.asarray(s_cur_den_parallel)
print('s_cur_den_parallel=',s_cur_den_parallel)

s_cur_den_onetwo=dat.variables['s_cur_den_onetwo']
print('s_cur_den_onetwo is: ', s_cur_den_onetwo.long_name, s_cur_den_onetwo.shape)
s_cur_den_toroidal=dat.variables['s_cur_den_toroidal']
print('s_cur_den_toroidal is: ', s_cur_den_toroidal.long_name, s_cur_den_toroidal.shape)
s_cur_den_poloidal=dat.variables['s_cur_den_poloidal']
print('s_cur_den_poloidal is: ', s_cur_den_poloidal.long_name, s_cur_den_poloidal.shape)

# Power profiles [erg/sec/cm^3]
powden=dat.variables['powden']
print('powden is: ',  powden.units, powden.shape)
powden_e=dat.variables['powden_e']
print('powden_e is: ', powden_e.units,powden_e.shape)
powden_i=dat.variables['powden_i']
print('powden_i is: ', powden_i.units, powden_i.shape)

# For iabsorp.eq.3 only:
#powden_s=dat.variables['powden_s']
#print 'powden_s is: ', powden_s.long_name, powden_s.shape

# Total power [erg/sec]
power_total=dat.variables['power_total']
print('power_total is: ', power_total.long_name, power_total.shape)
power_inj_total=dat.variables['power_inj_total']
print('power_inj_total is: ', power_inj_total.long_name, power_inj_total.shape)
powtot_e=dat.variables['powtot_e']
print('powtot_e is: ', powtot_e.long_name, powtot_e.shape)
powtot_i=dat.variables['powtot_i']
print('powtot_i is: ', powtot_i.long_name, powtot_i.shape)
powtot_cl=dat.variables['powtot_cl']
print('powtot_cl is: ', powtot_cl.long_name, powtot_cl.shape)
# For iabsorp.eq.3 only:
#powtot_s=dat.variables['powtot_s']
#print 'powtot_s is: ', powtot_s.long_name, powtot_s.shape

ws=dat.variables['ws'] #   pol.distance along ray [cm]
ws_min=np.min(ws)
ws_max=np.max(ws)
ws_step= np.ceil((ws_max*1.1-ws_min)/4)  # step for xticks() labeling ;
                     #    ()/4 means 4 labels only
print('ws is: ', ws.long_name, ws.shape)
print('min/max of ws [cm]', ws_min,ws_max)

delpwr=dat.variables['delpwr']
print('delpwr is: ', delpwr.long_name, delpwr.shape)

spsi=dat.variables['spsi'] # normalized small radius=rho given by indexrho
print('spsi is: ', spsi.long_name, spsi.shape)
print('min/max of spsi', np.min(spsi),np.max(spsi))


# RAY trajectories
wr=dat.variables['wr']
print('wr is: ', wr.long_name, wr.shape)
wz=dat.variables['wz']
print('wz is: ', wz.long_name, wz.shape)
wphi=dat.variables['wphi']
print('wphi is: ', wphi.long_name, wphi.shape)
# Number of elements for each ray:
nrayelt=dat.variables['nrayelt']
print('nrayelt is: ', nrayelt.long_name, nrayelt.shape)
# Refractive indices along rays:
wnper=dat.variables['wnper']
print('wnper is: ', wnper.long_name, wnper.shape)
wnpar=dat.variables['wnpar']
print('wnpar is: ', wnpar.long_name, wnpar.shape)
wn_r=dat.variables['wn_r']
print('wn_r is: ', wn_r.long_name, wn_r.shape)
wn_z=dat.variables['wn_z']
print('wn_z is: ', wn_z.long_name, wn_z.shape)
wn_phi=dat.variables['wn_phi']
print('wn_phi is: ', wn_phi.long_name, wn_phi.shape)
# E-wave-field along rays
cwexde=dat.variables['cwexde']
print('cwexde is: ', cwexde.long_name, cwexde.shape)
cweyde=dat.variables['cweyde']
print('cweyde is: ', cweyde.long_name, cweyde.shape)
cwezde=dat.variables['cwezde']
print('cwezde is: ', cwezde.long_name, cwezde.shape)
# fluxn = Power Flux along rays:
fluxn=dat.variables['fluxn']
print('fluxn is: ', fluxn.long_name, fluxn.shape, fluxn.units)
# Total magnetic field along rays:
sbtot=dat.variables['sbtot']
print('sbtot is: ',sbtot.long_name, sbtot.shape, sbtot.units)
# density along rays [1/cm^3]
sene=dat.variables['sene']
print('sene is: ', sene.long_name, sene.shape, sene.units)
# Ki along rays
salphal=dat.variables['salphal']
print('salphal is: ', salphal.long_name, salphal.shape)

# Vgroup normalized to c:
vgr_z=dat.variables['vgr_z']
print('vgr_z is: ', vgr_z.long_name, vgr_z.shape)
vgr_r=dat.variables['vgr_r']
print('vgr_r is: ', vgr_r.long_name, vgr_r.shape)
vgr_phi=dat.variables['vgr_phi']
print('vgr_phi is: ', vgr_phi.long_name, vgr_phi.shape)

sb_z=dat.variables['sb_z']
sb_phi=dat.variables['sb_phi']
print 'sb_z[0,0]=  ', sb_z[0,0], ' Gauss'
print 'sb_phi[0,0]=', sb_phi[0,0], ' Gauss'

z_starting=dat.variables['z_starting'] #starting Z[m], each ray
print('z_starting is: ', z_starting.long_name, z_starting.shape)
z_starting=np.asarray(z_starting)
print 'z_starting[m]=',z_starting

r_starting=dat.variables['r_starting'] #starting R[m], each ray
print('r_starting is: ', r_starting.long_name, r_starting.shape)
r_starting=np.asarray(r_starting)
print 'r_starting[m]=',r_starting

phi_starting=dat.variables['phi_starting'] #starting phi[degree], each ray
print('phi_starting is: ', phi_starting.long_name, phi_starting.shape)
phi_starting=np.asarray(phi_starting)
print 'phi_starting[degree]=',phi_starting

if istart==1:
    alphast=dat.variables['alphast'] #aiming toroidal angle[degree] measured from 
    # R-vector through source;  each ray
    print('alphast is: ', alphast.long_name, alphast.shape)
    alphast=np.asarray(alphast)
    print 'alphast[degree]=',alphast

    betast=dat.variables['betast'] #aiming poloidal angle[degree] measured from 
    # z=constant plane, pos above plane, neg below;  each ray
    print('betast is: ', betast.long_name, betast.shape)
    betast=np.asarray(betast)
    print 'betast[degree]=',betast

# Number of rays
Nrays=wz[:,0].size  
print('Number of rays: Nrays=',Nrays)
    
print('----------------------------------------')
#BH print power_inj_total.long_name, power_inj_total[0], power_inj_total.units
#BH print power_total.long_name, power_total[0], power_total.units
#BH print powtot_e.long_name,  powtot_e[0],  powtot_e.units
#BH print powtot_i.long_name,  powtot_i[0],  powtot_i.units
#BH print powtot_cl.long_name, powtot_cl[0], powtot_cl.units
print(power_inj_total.long_name, power_inj_total.getValue(), power_inj_total.units)
print(power_total.long_name, power_total.getValue(), power_total.units)
print(powtot_e.long_name,  powtot_e.getValue(),  powtot_e.units)
print(powtot_i.long_name,  powtot_i.getValue(),  powtot_i.units)
print(powtot_cl.long_name, powtot_cl.getValue(), powtot_cl.units)
print('----------------------------------------')

# Convert erg/sec to kW
power_total=     np.asscalar(power_total[:])/1.e10    
power_inj_total= np.asscalar(power_inj_total[:])/1.e10
powtot_e=        np.asscalar(powtot_e[:])/1.e10
powtot_i=        np.asscalar(powtot_i[:])/1.e10
powtot_cl=       np.asscalar(powtot_cl[:])/1.e10


#===================  PLOTS =============================================
#set fonts and line thicknesses
# If this part with Params causes problems, comment it out.
# (May depend on Python version)
mpl.rcParams['axes.linewidth']=linw
mpl.rcParams['lines.linewidth']=linw
mpl.rcParams['axes.labelsize']=fnt+4
mpl.rcParams['font.size']=fnt+4
mpl.rcParams['legend.fontsize']=fnt
mpl.rcParams['xtick.labelsize']=fnt
mpl.rcParams['ytick.labelsize']=fnt
#mpl.rcParams['xtick.linewidth']=linw
#mpl.rcParams['ytick.linewidth']=linw
mpl.rcParams['font.weight']='regular'
#mpl.rcParams['format']='%1.1e'
mpl.rcParams['font.size']=fnt+2  # set font size for text in mesh-plots



#----------------------------------------------------------------T(rho), n(rho)
fig0=plt.figure()
# Plot with different thickness of lines, depending on species number
linw_mx= 4*linw # largest thickness
linw_mn= 0.5    # smallest
if Nsp>1:
    dlinw= (linw_mx-linw_mn)/(Nsp-1)
else:
    dlinw=0
plt.subplot(221)
# plt.hold(True)
#plt.title('$ T_e, T_i $')
#plt.xlabel(r'$\rho$')
plt.title('$ T $   $ (keV) $',fontsize=fnt+4,y=1.03)
plt.xlabel(r'$\rho$',fontsize=fnt+4)
plt.grid(True)
for i in range(0,Nsp,1):
    if remainder(i,4)==0: col='k'
    if remainder(i,4)==1: col='b'
    if remainder(i,4)==2: col='r'
    if remainder(i,4)==3: col='g'    
    plt.plot(rho_bin[:],temprof[i,:],color=col,linewidth=linw_mx-i*dlinw)

plt.subplot(223) #-------------------------
dx=0.7/Nsp
for i in range(0,Nsp,1):
    if remainder(i,4)==0: col='k'
    if remainder(i,4)==1: col='b'
    if remainder(i,4)==2: col='r'
    if remainder(i,4)==3: col='g'  
    txt1=r"$m/m_e =$"+r"$%1.0f$" %(mass[i])
    txt2=r"$q/e =$"+r"$%1.0f$" %(charge[i])
    txt3=r"$T_0[keV] =$"+r"$%1.3f$" %(temprof[i,0])
    txt4=r"$n_0[cm^{-3}] =$"+r"$%1.5e$" %(densprof[i,0])
    plt.plot(0.01, 0.7-dx*i, 'o', color=col)
    plt.text(0.04, 0.7-dx*i, txt1 , va='center',fontsize=fnt+4)
    plt.text(0.61, 0.7-dx*i, txt2 , va='center',fontsize=fnt+4) 
    plt.text(0.93, 0.7-dx*i, txt3 , va='center',fontsize=fnt+4) 
    plt.text(1.55, 0.7-dx*i, txt4 , va='center',fontsize=fnt+4) 
    plt.axis([0., 1., 0., 1.])
    plt.axis('off')

plt.subplot(222)
# plt.hold(True)
#plt.title('$ n_e, n_i $')
plt.xlabel(r'$\rho$',fontsize=fnt+4)
plt.title('$ n $   $ (cm^{-3}) $',fontsize=fnt+4,y=1.03)
plt.grid(True)
for i in range(0,Nsp,1):
    if remainder(i,4)==0: col='k'
    if remainder(i,4)==1: col='b'
    if remainder(i,4)==2: col='r'
    if remainder(i,4)==3: col='g'    
    plt.plot(rho_bin[:],densprof[i,:],color=col,linewidth=linw_mx-i*dlinw)   

savefig('genray_profiles_T-n.png',format='png') # try pdf,eps,ps
show() 
#print 'rho_bin=',rho_bin
#print 'temprof[0,:]=',temprof[0,:]
#stop



#--------------------------------------------------------------------------
fig4=plt.figure()  # RAYS in top-view (R-phi or X-Y)
ax = plt.subplot(1, 1, 1)
ax.set_aspect(1.0)
# plt.hold(True)
locs, labels = plt.xticks() 
plt.setp(labels, rotation=90)
plt.title('$Rays$  $X(t),Y(t)$',y=1.03)
plt.xlabel('$X$  $(cm)$')
plt.ylabel('$Y$  $(cm)$')
plt.grid(True)

# Define boundary for top-view (R-phi) plots:
phi=arange(301.)*2*pi/300.  # tor.angle in rad.
bndry_in_X= Rmin_plot*cos(phi)
bndry_in_Y= Rmin_plot*sin(phi)
bndry_out_X= Rmax_plot*cos(phi)
bndry_out_Y= Rmax_plot*sin(phi)
print 'Rmax_plot=',Rmax_plot

if n_wall>4:  # plot walls, if any:
    plt.plot(np.max(r_wall)*100*cos(phi),np.max(r_wall)*100*sin(phi),'k',linewidth=linw*2)
    plt.plot(np.min(r_wall)*100*cos(phi),np.min(r_wall)*100*sin(phi),'k',linewidth=linw*2)
if ncontr>4:  # plot LCFS, if any
    plt.plot(np.max(rcontr)*100*cos(phi),np.max(rcontr)*100*sin(phi),'g',linewidth=linw*2)
    plt.plot(np.min(rcontr)*100*cos(phi),np.min(rcontr)*100*sin(phi),'g',linewidth=linw*2)
if nlimiter>4: # plot limiter surface, if any
    plt.plot(np.max(rlimiter)*100*cos(phi),np.max(rlimiter)*100*sin(phi),'b',linewidth=linw*2)
    plt.plot(np.min(rlimiter)*100*cos(phi),np.min(rlimiter)*100*sin(phi),'b',linewidth=linw*2)
plt.plot(R_axis*cos(phi),R_axis*sin(phi),'k--',linewidth=linw*0.5)

for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    print(i,Nm)
    if Nm>0:
        X=multiply(cos(wphi[i,0:Nm]),wr[i,0:Nm])
        Y=multiply(sin(wphi[i,0:Nm]),wr[i,0:Nm])
        plt.plot(X,Y,color=col,linewidth=linw)  
        # plot arrow for refractive vector (Nr,Nz) at starting point: 
        #Nr0=wn_r[i,0]
        #Ny0=wn_y[i,0]
        #N0=sqrt(Nr0*Nr0+Ny0*Ny0)
        #plt.arrow(wr[i,0],wy[i,0],(Nr0/N0)*arr_len,(Ny0/N0)*arr_len,'->') 
        plt.plot(X[0],Y[0],'k.') # small circle at launching point
plt.plot(bndry_in_X,  bndry_in_Y,  'b',linewidth=linw*2)
plt.plot(bndry_out_X, bndry_out_Y, 'b',linewidth=linw*2)

if (sb_phi[0,0]<0):   
    # Bphi<0 plot marker pointing DOWN
    #plt.plot(R_LCFS_mn,0,'r^',markersize=10)    
    text(R_LCFS_min,0,r'$\downarrow$ $B_{\phi}$')
if (sb_phi[0,0]>0):   
    # Bphi>0 plot marker pointing UP
    #plt.plot(R_LCFS_mn,0,'rv',markersize=10)    
    text(R_LCFS_min,0,r'$\uparrow$ $B_{\phi}$')
  
#plt.axis('off')
plt.savefig('genray_rays_in_R-phi.png') 
show()


#--------------------------------------------------------------------------
fig3=plt.figure()  # w/wc resonances and Pol.Flux in cross-section view R-Z
ax = plt.subplot(1, 1, 1)
ax.set_aspect(1.0)
ax.axis([xmin,xmax,zmin,zmax])
# plt.hold(True)
plt.xlabel('$R$  $(cm)$')
locs, labels = plt.xticks() 
plt.setp(labels, rotation=90)
plt.ylabel('$Z$  $(cm)$')
plt.grid(True)

if n_wall>4:  # plot walls, if any:
    plt.plot(multiply(r_wall,100),multiply(z_wall,100),'k',linewidth=linw*2)
if ncontr>4:  # plot LCFS, if any
    plt.plot(multiply(rcontr,100),multiply(zcontr,100),'g',linewidth=linw*2)
if nlimiter>4: # plot limiter surface, if any
    plt.plot(multiply(rlimiter,100),multiply(zlimiter,100),'b',linewidth=linw*2)
    
if n_eqdsk==1: # plot contour lines of poloidal_flux/2pi
    R,Z = np.meshgrid(eqdsk_r, eqdsk_z)
    PSI = eqdsk_psi
    CS=plt.contour( R*100,Z*100,PSI,levels=Ncont,linewidths=1,cmap=plt.cm.get_cmap("hot"))
    #CB=plt.colorbar(orientation='vertical', shrink=0.5, format='%.2e')
if i_eqdsk==1: # Plot resonance layers for electrons or ions
    if n_eqdsk==0: # There was no data on PSI. Use data from eqdsk
        Req,Zeq = np.meshgrid(er,ez) # 2D grids [cm]
        PSI = epsi
        #CS=plt.contour( Req,Zeq,PSI,Ncont,linewidths=1,cmap=plt.cm.hot)
        #CB=plt.colorbar(orientation='vertical', shrink=0.5, format='%.2e')
        # Levels of PSI (to mark psimag and psilim)
        levels2=np.arange(psimag,psilim,(psilim-psimag)/30)
        CS=plt.contour(Req,Zeq,PSI,levels=levels2,linewidth=2,cmap=plt.cm.get_cmap("jet"))
    # B(r,z) was defined from eqdsk data
    #print shape(B), shape(R), n_eqdsk
    WWce0= f/abs(28e5*btor)  # omega/omega_ce at mag.axis
    print('omega/omega_ce at mag.axis :' , WWce0)
    print('Nsp=',Nsp)
    if (Nsp==1 or WWce0>1.0):  # Plot res.layers for electrons
        isp=0  # electrons
        level0= ceil(WWce0)
    if Nsp>1:
        #isp=1  # Plot res.layers for ions#1 if present
        isp=isp_wc-1 # ion (isp_wc is the genray index; isp is the python index)
        qm=  charge[isp]/mass[isp]  # q/m     
        WWci0= f/abs(28e5*qm*btor)  # omega/omega_ci
        print('isp=',isp)
        print('omega/omega_ci at mag.axis :' , WWci0)
        level0= ceil(WWci0)
    qm=  charge[isp]/mass[isp]  # q/m     
    WWc= f/abs(28e5*qm*B)  # omega/omega_c
    #level_mn= 1 #max(1,level0-50)
    #level_mx= 40 #level_mn+100
    levels=np.arange(level_mn,level_mx,1)
    CS=plt.contour(R,Z,WWc,levels=levels,linewidths=3,cmap=plt.cm.get_cmap("jet"))
    plt.clabel(CS,inline=1, fmt='%1.0f',fontsize=fnt+1, colors='k',rightside_up=False)
    #CB=plt.colorbar(orientation='vertical', shrink=0.4,format='%.2f') #colorbar
    plot(R_axis,Z_axis,'k+') # plot '+' at magnetic axis
    txt1= r"  $m/m_e=$"+r"$%1.0f$" %(mass[isp])
    txt2= r"  $q/e=$"+r"$%1.0f$" %(charge[isp])
    plt.title('$\omega / \omega_{c}$'+'  $for$'+txt1+txt2, y=1.03)

plt.savefig('genray_wwc_inRZ.png') 
show()



#--------------------------------------------------------------------------
fig0=plt.figure() # RAYS in cross-sectional view R-Z
ax = plt.subplot(111)
ax.set_aspect(1.0)
# plt.hold(True)    
ax.axis([xmin,xmax,zmin,zmax])
plt.title('$Rays$  $R(t),Z(t)$',y=1.03)
plt.xlabel('$R$  $(cm)$',y=1.02)
locs, labels = plt.xticks() 
plt.setp(labels, rotation=90)
plt.ylabel('$Z$  $(cm)$')
plt.grid(True)

if i_eqdsk==1: # Plot resonance layers for electrons or ions
    CS=plt.contour(R,Z,WWc,levels=levels,linewidths=2,cmap=plt.cm.get_cmap("jet"))
    plt.clabel(CS,inline=1, fmt='%1.0f',fontsize=fnt+1, colors='k',rightside_up=False)

if n_wall>4:  # plot walls, if any:
    plt.plot(multiply(r_wall,100),multiply(z_wall,100),'k',linewidth=linw*2)
if ncontr>4:  # plot LCFS, if any
    plt.plot(multiply(rcontr,100),multiply(zcontr,100),'g',linewidth=linw*2)
if nlimiter>4: # plot limiter surface, if any
    plt.plot(multiply(rlimiter,100),multiply(zlimiter,100),'b',linewidth=linw*2)    
    
if n_wall+ncontr+nlimiter <5: # no data on walls or LCFS; Try plotting PSI
    if n_eqdsk==1:
        R,Z = np.meshgrid(eqdsk_r, eqdsk_z)
        PSI = np.array(eqdsk_psi)
        CS=plt.contour( R*100,Z*100,PSI,Ncont,linewidths=1,cmap=plt.cm.get_cmap("hot"))
    if n_eqdsk==0: # There was no data on PSI. Use data from eqdsk
        Req,Zeq = np.meshgrid(er,ez) # 2D grids [cm]
        PSI = epsi
        #CS=plt.contour( Req,Zeq,PSI,Ncont,linewidths=1,cmap=plt.cm.hot)
        #CB=plt.colorbar(orientation='vertical', shrink=0.5, format='%.2e')
        # Levels of PSI (to mark psimag and psilim)
        levels2=np.arange(psimag,psilim*1.2,(psilim*1.2-psimag)/30)
        CS=plt.contour(Req,Zeq,PSI,levels=levels2,linewidth=2,cmap=plt.cm.get_cmap("jet"))

for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    #print i,Nm
    if Nm>0:
        if rays_style==0:
            #--->USUAL rays (fast plotting):
            plt.plot(wr[i,0:Nm],wz[i,0:Nm],color=col,linewidth=linw)
        if rays_style==1:
            #--->Rays with width corresponding to remaining power:
            width=5.0*(delpwr[i,0:Nm]/delpwr[i,0])  #**0.5
            plot_widths(wr[i,0:Nm],wz[i,0:Nm],width,color=col)
        if rays_style==2:
            #--->Rays with width corr. to absorbed power:
            width=(delpwr[i,0:Nm]/delpwr[i,0])
            width_diff=abs(np.diff(width)) #diff(width) is neg when power is reduced
            dwidth=100.0*width_diff/np.max(width_diff) 
            plt.scatter(wr[i,0:Nm],wz[i,0:Nm],\
              marker='o',edgecolor="red",facecolor="none",s=dwidth) # as markers
        # plot arrow for refractive vector (Nr,Nz) at starting point: 
        #Nr0=wn_r[i,0]
        #Nz0=wn_z[i,0]
        #N0=sqrt(Nr0*Nr0+Nz0*Nz0)
        #plt.arrow(wr[i,0],wz[i,0],(Nr0/N0)*arr_len,(Nz0/N0)*arr_len,'->') 
        # plot small circle at the launching point:
    if Nm>0:
        # plot small circle at the launching point:
        plt.plot(wr[i,0],wz[i,0],'k.')
        #print i+1, ' R0=',wr[i,0],  '   Z0=',wz[i,0]
if R_axis>0: plot(R_axis,Z_axis,'k+') # plot '+' at magnetic axis

R_LCFS_mn= np.min(rcontr)*100 # cm
if (sb_z[0,0]>0)&(wr[0,0]>R_axis):   
    # Bz>0 at R>R_axis which means at R<R_axis plot marker pointing DOWN
    #plt.plot(R_LCFS_mn,Z_axis,'rv',markersize=10)    
    text(R_LCFS_mn,Z_axis,r'$\downarrow$ $B_p$')
if (sb_z[0,0]<0)&(wr[0,0]>R_axis):   
    # Bz<0 at R>R_axis which means at R<R_axis plot marker pointing UP
    #plt.plot(R_LCFS_mn,Z_axis,'r^',markersize=10)    
    text(R_LCFS_mn,Z_axis,r'$\uparrow$ $B_p$')

# Additional structures (ports, coils) 
if istruct==1: 
    plot(Rstruct1, Zstruct1,'k')
    plot(Rstruct2, Zstruct2,'k') 
    plot(Rstruct3, Zstruct3,'k')
    plot(Rstruct3,-Zstruct3,'k') # mirror image
    plot(Rstruct4, Zstruct4,'k')
plt.savefig('genray_rays_inRZ.png') 
show()
#stop

#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
if i_ox==2:  # only for O-X transmission case
    fig8=plt.figure()
    plt.subplot(321) #-------------------------
    plt.hold(True)
    plt.grid(True)
    #plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
    plt.ylabel('$1-converted;$ $0-not$')
    plt.title('$Indicator$ $of$ $OX$ $conversion$')
    plt.ylim((0.,1.1))
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        plt.plot(wnpar[i,0],i_ox_conversion[i],'o',color=col)
        
    plt.subplot(322) #-------------------------
    plt.hold(True)
    plt.grid(True)
    #plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
    #text(wnpar[0,0],1.02,'   $pwr(Xmode)/pwr(Omode)$')
    plt.title('$Transm.coef.=P(X)/P(O)$')
    plt.ylim((0.,1.1))
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        plt.plot(wnpar[i,0],transm_ox[i],'o',color=col)

    plt.subplot(323) #-------------------------
    plt.hold(True)
    plt.grid(True)
    #plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
    plt.title('$N_{\perp,tang}$ $at$ $OX$ $tran.$',y=0.8)
    y_mn= min(0.,np.amin(cn_b_gradpsi[:]))
    y_mx= max(0.,np.amax(cn_b_gradpsi[:]))
    dy= 0.35*(y_mx-y_mn)
    plt.ylim((y_mn-dy,y_mx+dy))
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        plt.plot(wnpar[i,0],cn_b_gradpsi[i],'o',color=col)

    plt.subplot(324) #-------------------------
    plt.hold(True)
    plt.grid(True)
    #plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
    plt.title('$N_{||}$ $and$ $N_{||opt}$ $at$ $OX$ $tran.$',y=0.8)
    y_mn= min(0.,np.amin(cnpar_ox[:]))
    y_mx= max(0.,np.amax(cnpar_ox[:]))
    dy= 0.35*(y_mx-y_mn)
    plt.ylim((y_mn-dy,y_mx+dy))
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        plt.plot(wnpar[i,0],cnpar_ox[i],'o',color=col)
        if cnpar_ox[i]>0:
            plt.plot(wnpar[i,0],cn_par_optimal[i],'k.') #optimal N||
        else:
            plt.plot(wnpar[i,0],-cn_par_optimal[i],'k.') #optimal N||
            # cn_par_optimal is always as |cn_par_optimal|
            # so we choose the proper sign from cnpar_ox
                
    if istart==1:    
        plt.subplot(325) #-------------------------
        plt.hold(True)
        plt.grid(True)
        plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
        plt.title('$alphast$ $[degree]$',horizontalalignment='center',y=0.8)
        y_mn= np.amin(alphast[:])-1.0
        y_mx= np.amax(alphast[:])+1.0
        dy= 0.35*(y_mx-y_mn)
        plt.ylim((y_mn-dy,y_mx+dy))
        for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
            if remainder(i,6)==0: col='b'
            if remainder(i,6)==1: col='g'
            if remainder(i,6)==2: col='r'
            if remainder(i,6)==3: col='c'    
            if remainder(i,6)==4: col='m' 
            if remainder(i,6)==5: col='k'  
            plt.plot(wnpar[i,0],alphast[i],'o',color=col)
        
        plt.subplot(326) #-------------------------
        plt.hold(True)
        plt.grid(True)
        plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
        plt.title('$betast$ $[degree]$',horizontalalignment='center',y=0.8)
        y_mn= np.amin(betast[:])-1.0
        y_mx= np.amax(betast[:])+1.0
        dy= 0.35*(y_mx-y_mn)
        plt.ylim((y_mn-dy,y_mx+dy))
        for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
            if remainder(i,6)==0: col='b'
            if remainder(i,6)==1: col='g'
            if remainder(i,6)==2: col='r'
            if remainder(i,6)==3: col='c'    
            if remainder(i,6)==4: col='m' 
            if remainder(i,6)==5: col='k'  
            plt.plot(wnpar[i,0],betast[i],'o',color=col)
          
    plt.savefig('genray_rays_OX_transmission.png') 
    show()
#stop

#--------------------------------------------------------------------------
fig80=plt.figure()  # E-wave-field along RAYS vs pol.distance(t)
plt.subplot(231) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$|E_X/E|$ $along$ $ray$',y=1.03)
plt.xlim( (ws_min, ws_max) )
plt.ylim( (0., 1.05) )  
xticks(np.arange(ws_min, ws_max, step=ws_step) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        ereal=cwexde[0,i,0:Nm].copy()
        eimag=cwexde[1,i,0:Nm].copy()
        ea= sqrt(ereal**2 + eimag**2)
        plt.plot(ws[i,0:Nm],ea[0:Nm],color=col,linewidth=linw)  
plt.subplot(232) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$|E_Y/E|$ $along$ $ray$',y=1.03)
plt.ylim( (0., 1.05) ) 
plt.xlim( (ws_min, ws_max) )
xticks(np.arange(ws_min, ws_max, step=ws_step) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        ereal=cweyde[0,i,0:Nm].copy()
        eimag=cweyde[1,i,0:Nm].copy()
        ea= sqrt(ereal**2 + eimag**2)
        plt.plot(ws[i,0:Nm],ea[0:Nm],color=col,linewidth=linw)  
plt.subplot(233) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$|E_Z/E|$ $along$ $ray$',y=1.03)
plt.ylim( (0., 1.05) ) 
plt.xlim( (ws_min, ws_max) )
xticks(np.arange(ws_min, ws_max, step=ws_step) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        ereal=cwezde[0,i,0:Nm].copy()
        eimag=cwezde[1,i,0:Nm].copy()
        ea= sqrt(ereal**2 + eimag**2)
        plt.plot(ws[i,0:Nm],ea[0:Nm],color=col,linewidth=linw)
plt.subplot(234) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.ylabel('$Power$ $Flux$ $(erg/s/cm^2)$')
plt.xlim( (ws_min, ws_max) )
xticks(np.arange(ws_min, ws_max, step=ws_step) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],fluxn[i,0:Nm],color=col,linewidth=linw)
plt.xlabel('$poloidal$ $distance$ $along$ $ray$  $(cm)$')
plt.subplot(236) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.ylabel('$Lin.Damping$ $k_i$ $(cm^{-1})$')
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],salphal[i,0:Nm],color=col,linewidth=linw)
plt.xlabel('$poloidal$ $distance$ $along$ $ray$  $(cm)$')
plt.savefig('genray_rays_Ewave_p.png') 
show()


fig8=plt.figure()  # E+ and E- wave-field along RAYS vs distance(t)
plt.subplot(611) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$|E_{+}/E|$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    one_ov_root2=1./sqrt(2.)
    if Nm>0:
        Ex_real=cwexde[0,i,0:Nm].copy()
        Ex_imag=cwexde[1,i,0:Nm].copy()
        Ey_real=cweyde[0,i,0:Nm].copy()
        Ey_imag=cweyde[1,i,0:Nm].copy()
        # Form Eplus polarization as (Ex+i*Ey)/sqrt(2)  
        # WARNING: are Ex and Ey actually perp to Beq ?  BH: Yes. Stix Frame.
        # This normalization gives that (Eplus_abs**2+Emin**2+E_Z**2) =1.)
        Eplus_real= one_ov_root2*(Ex_real - Ey_imag)  #
        Eplus_imag= one_ov_root2*(Ex_imag + Ey_real)
        Eplus_abs= sqrt(Eplus_real**2 + Eplus_imag**2)
        plt.plot(ws[i,0:Nm],Eplus_abs,color=col,linewidth=linw)  
plt.subplot(612) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$|E_{-}/E|$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        Ex_real=cwexde[0,i,0:Nm].copy()
        Ex_imag=cwexde[1,i,0:Nm].copy()
        Ey_real=cweyde[0,i,0:Nm].copy()
        Ey_imag=cweyde[1,i,0:Nm].copy()
        # Form Eminus as (Ex-i*Ey)/sqrt(2.)
        #   WARNING: are Ex and Ey actually perp to Beq ?
        Eminus_real= one_ov_root2*(Ex_real + Ey_imag)  #
        Eminus_imag= one_ov_root2*(Ex_imag - Ey_real)
        Eminus_abs= sqrt(Eminus_real**2 + Eminus_imag**2)
        plt.plot(ws[i,0:Nm],Eminus_abs,color=col,linewidth=linw)  
plt.subplot(613) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$|E_Z/E|$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        ereal=cwezde[0,i,0:Nm].copy()
        eimag=cwezde[1,i,0:Nm].copy()
        ea= sqrt(ereal**2 + eimag**2)
        plt.plot(ws[i,0:Nm],ea,color=col,linewidth=linw)
plt.subplot(614) #-------------------------
#To check sum =1. (Yes):
#print 'sum polorizations**2=', Eplus_abs**2+Eminus_abs**2+ea**2
plt.hold(True)
isp=0 # electrons, by default
if Nsp>0:
    isp=isp_wc-1 # ion (isp_wc is the genray index; isp is the python index)
    msme= mass[isp]/mass[0]        # m_s/m_e ratio
    Z_s= charge[isp]/charge[0]    # Z_s/e
    print 'm_s/m_e=', msme, '  Z_s=',Z_s
    print 'isp_name=',isp_name
    plt.subplot(614) #------------------------- for 1st ion species
    plt.hold(True)
    plt.grid(True)
    plt.ylabel('$\omega/\omega_{c}$['+isp_name[isp]+']')
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
        if Nm>0:
            plt.plot(ws[i,0:Nm],f/(28e5*sbtot[i,0:Nm]*Z_s/msme),color=col,linewidth=linw)
plt.subplot(615) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$K_i$'+' $(cm^{-1})$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],salphal[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(616) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel(r'$N_{\perp}$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],wnper[i,0:Nm],color=col,linewidth=linw)
      
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
plt.savefig('genray_rays_Ewave_plus-min_p.png') 
show()    



#--------------------------------------------------------------------------
fig81=plt.figure()  # Vgroup/c along RAYS vs pol.distance(t)
plt.subplot(231) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$V_{grp,R}/c$  $along$ $ray$',y=1.03)
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],vgr_r[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(232) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$V_{grp,\phi}/c$  $along$ $ray$',y=1.03)
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],vgr_phi[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(233) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$V_{grp,Z}/c$  $along$ $ray$',y=1.03)
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],vgr_z[i,0:Nm],color=col,linewidth=linw)
plt.subplot(235) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.ylabel('$|V_{grp}/c|$  $along$ $ray$')
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        vgrr=vgr_r[i,0:Nm]
        vgrphi=vgr_phi[i,0:Nm]
        vgrz=vgr_z[i,0:Nm]
        vgr=sqrt(vgrr**2+vgrphi**2+vgrz**2)
        plt.plot(ws[i,0:Nm],vgr[0:Nm],color=col,linewidth=linw)
plt.xlabel('$poloidal$ $distance$ $along$ $ray$  $(cm)$')
plt.savefig('genray_rays_Vgroup_p.png') 
show()



#--------------------------------------------------------------------------
fig5=plt.figure()  # wce/w and (wpe/w)^2 along RAYS vs distance
qm=  charge[0]/mass[0]  # First species is electrons
q2m= charge[0]**2/mass[0]

plt.subplot(231) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$|\omega_{ce}/\omega|$')
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        wce_w= abs(28e5*qm*sbtot[i,0:Nm])/f
        plt.plot(ws[i,0:Nm],wce_w,color=col,linewidth=linw)
        
plt.subplot(234) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$|\omega/\omega_{ce}|$')
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        w_wce= f/abs(28e5*qm*sbtot[i,0:Nm])
        plt.plot(ws[i,0:Nm],w_wce,color=col,linewidth=linw)
plt.xlabel('$poloidal$ $distance$ $along$ $ray$  $(cm)$')

plt.subplot(232) #------------------------- value of RF frequency
txt=r"$f (MHz) =$"+r"$%1.3f$" %(f*1e-6)
plt.title(txt,y=1.03)   # 
plt.axis('off')

plt.subplot(233) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$(\omega_{pe}/\omega)^2$',y=1.03)
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        #print 806.2e5*sene[0,0]*q2m/f**2, sene[0,0], q2m
        plt.plot(ws[i,0:Nm],806.2e5*sene[i,0:Nm]*q2m/f**2,color=col,linewidth=linw)    
#plt.xlabel('$pol.distance$ $along$ $ray$  $(cm)$')

if Nsp>1: # plot w_c/w along rays
    plt.subplot(235) #-------------------------
    plt.grid(True)
    #isp=1  # First ion species
    isp=isp_wc-1 # ion (isp_wc is the genray index; isp is the python index)
    qm=  charge[isp]/mass[isp]  # First ion species
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
        if Nm>0:
            plt.plot(ws[i,0:Nm],abs(28e5*qm*sbtot[i,0:Nm]/f),color=col,linewidth=linw)
    plt.xlabel('$pol.distance$ $along$ $ray$  $(cm)$')
    txt1=r"  $m/m_e=$"+r"$%1.0f$" %(mass[isp])
    txt2=r"  $q/e=$"+r"$%1.0f$" %(charge[isp])
    plt.ylabel('$For$'+txt1+txt2)
    plt.title('$|\omega_{ci}/\omega|$',y=1.03)
    xticks(np.arange(ws_min, ws_max, step=ws_step) )
    plt.xlim( (ws_min, ws_max) )
    
#--------------------------------------------
if Nsp>1:
    isp=isp_wc-1 # ion (isp_wc is the genray index; isp is the python index)
    mime= mass[isp]/mass[0]        # m_s/m_e ratio
    Z_s= charge[isp]/charge[0]    # Z_s/e
    print 'm_s/m_e=', mime, '  Z_s=',Z_s
    plt.subplot(236) #------------------------- for 1st ion species
    plt.hold(True)
    plt.grid(True)
    plt.title('$\omega/\omega_{c}$ ['+isp_name[isp]+']')
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
        linwi=2*linw-0.5*i
        linwi=max(0.5,linwi)
        plt.plot(ws[i,0:Nm],f/(28e5*sbtot[i,0:Nm]*Z_s/mime),color=col,linewidth=linwi)
    plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
else:
    plt.subplot(236) #------------------------- w_UH for electrons
    plt.hold(True)
    plt.grid(True)
    plt.title('$\omega_{UHe}/\omega$' )
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
        wuh_w = sqrt( (28e5*sbtot[i,0:Nm]/f)**2  +  806.2e5*sene[i,0:Nm]*q2m/f**2)
        plt.plot(ws[i,0:Nm],wuh_w, color=col,linewidth=linw)
    plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
    
plt.savefig('genray_rays_wc_wp_p.png') 
show()



#--------------------------------------------------------------------------
fig5=plt.figure()  # Refractive indices along RAYS vs R(t)
ax=plt.subplot(231) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$n_R$',y=1.03)
locs, labels = plt.xticks() 
plt.setp(labels, rotation=90)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(wr[i,0:Nm-1],wn_r[i,0:Nm-1],color=col,linewidth=linw) 
        plt.plot(wr[i,0],wn_r[i,0],'k.') # dot at t=0        
plt.subplot(232) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$n_Z$',y=1.03)
locs, labels = plt.xticks() 
plt.setp(labels, rotation=90)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(wr[i,0:Nm],wn_z[i,0:Nm],color=col,linewidth=linw)  
        plt.plot(wr[i,0],wn_z[i,0],'k.') # dot at t=0
plt.subplot(233) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title(r'$n_\phi$',y=1.03)
locs, labels = plt.xticks() 
plt.setp(labels, rotation=90)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(wr[i,0:Nm],wn_phi[i,0:Nm],color=col,linewidth=linw)
        plt.plot(wr[i,0],wn_phi[i,0],'k.') # dot at t=0
plt.subplot(234) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.ylabel(r'$n_{\perp}$'+' $along$ $ray$')
plt.xlabel('$R$ $along$ $ray$  $(cm)$')
locs, labels = plt.xticks() 
plt.setp(labels, rotation=90)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(wr[i,0:Nm],wnper[i,0:Nm],color=col,linewidth=linw)
        plt.plot(wr[i,0],wnper[i,0],'k.') # dot at t=0
plt.subplot(235) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title(r'$n_{||}$',verticalalignment='center',y=0.8)
plt.xlabel('$R$ $along$ $ray$  $(cm)$')
locs, labels = plt.xticks() 
plt.setp(labels, rotation=90)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(wr[i,0:Nm],wnpar[i,0:Nm],color=col,linewidth=linw)
        plt.plot(wr[i,0],wnpar[i,0],'k.') # dot at t=0
plt.subplot(236) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$lin.damp.$ $k_i$$(cm^{-1})$',y=0.8)
plt.xlabel('$R$ $along$ $ray$  $(cm)$')
locs, labels = plt.xticks() 
plt.setp(labels, rotation=90)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(wr[i,0:Nm],salphal[i,0:Nm],color=col,linewidth=linw)
        plt.plot(wr[i,0],salphal[i,0],'k.') # dot at t=0
plt.savefig('genray_rays_refr-index_R.png') 
show()



#--------------------------------------------------------------------------
fig6=plt.figure()  # Refractive indices along RAYS vs pol.distance(t)
plt.subplot(231) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$n_R$',y=1.03)
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],wn_r[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(232) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$n_Z$',y=1.03)
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],wn_z[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(233) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title(r'$n_\phi$',y=1.03)
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],wn_phi[i,0:Nm],color=col,linewidth=linw)
plt.subplot(234) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.ylabel(r'$n_{\perp}$')
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],wnper[i,0:Nm],color=col,linewidth=linw)
plt.subplot(235) #-------------------------
# plt.hold(True)
plt.grid(True)
#plt.title(r'$n_{||}$'+' $along$ $ray$',verticalalignment='center',y=1.03)
plt.title(r'$n_{||}$',verticalalignment='center',y=1.03)
plt.xlabel('$poloidal$ $distance$ $along$ $ray$  $(cm)$')
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],wnpar[i,0:Nm],color=col,linewidth=linw)
plt.subplot(236) #-------------------------
# plt.hold(True)
plt.grid(True)
#plt.title(r'$|n|$'+' $along$ $ray$',verticalalignment='center',y=1.03)
plt.title(r'$|n|$',verticalalignment='center',y=1.03)
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-nskip # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],sqrt(wnpar[i,0:Nm]**2+wnper[i,0:Nm]**2),color=col,linewidth=linw)
plt.savefig('genray_rays_refr-index_p.png') 
show()


# Renormalize power to the value used in experiment:
pwrscale=1.0 #1.9 # multiplication factor for input power
# (powers p and currents j will be renormalized by pwrscale)
power_total= power_total*pwrscale
power_inj_total= power_inj_total*pwrscale
powtot_e= powtot_e*pwrscale
powtot_i= powtot_i*pwrscale
if n_powden_cl>0:
    powtot_cl= powtot_cl*pwrscale

powden_e= np.asarray(powden_e)*pwrscale
powden_i= np.asarray(powden_i)*pwrscale
if n_powden_cl>0:
    powden_cl=np.asarray(powden_cl)*pwrscale
powden=   np.asarray(powden)*pwrscale
s_cur_den_parallel=np.asarray(s_cur_den_parallel)*pwrscale
s_cur_den_onetwo=  np.asarray(s_cur_den_onetwo)*pwrscale
s_cur_den_toroidal=np.asarray(s_cur_den_toroidal)*pwrscale
s_cur_den_poloidal=np.asarray(s_cur_den_poloidal)*pwrscale
delpwr= np.asarray(delpwr)*pwrscale

# factor 1e-4 is to convert erg/s/cm^3 -> kW/m^3
powden=  np.asarray(powden)/1.0e4 # kW/m^3
powden_e=np.asarray(powden_e)/1.0e4 # kW/m^3
powden_i=np.asarray(powden_i)/1.0e4 # kW/m^3
if n_powden_cl>0:
    powden_cl=np.asarray(powden_cl)/1.0e4 # kW/m^3
rhomax=1.05*np.amax(rho_bin_center)

#------------


#--------------------------------------------------------------------------
fig7=plt.figure()   # Power profiles
plt.subplot(221) #-------------------------
plt.ylabel('$p_e$ $(kW/m^3)$') 
# plt.hold(True)
plt.grid(True)
txt="$P_{e}$"+"$=$%3.3f" %(powtot_e) +" $kW$"
plt.title(txt,y=1.02)
plt.plot(rho_bin_center,(powden_e),'r',linewidth=linw)
axis([0.,rhomax,0.,1.05*np.amax(powden_e)+1.e-2])

plt.subplot(222) #-------------------------
# plt.hold(True)
plt.grid(True)
txt="$P_{i}$"+"$=$%3.3f" %(powtot_i) +" $kW$"
plt.title(txt,y=1.02)
plt.plot(rho_bin_center,(powden_i),'b',linewidth=linw*2)
#plt.title('$p_i$ $(Maxw.$ $ions)$ $and$ $p_e$ $(black)$')
#plt.plot(rho_bin_center,np.asarray(powden_e)*1.e-7,'k',linewidth=linw)
axis([0.,rhomax,0.,1.05*np.amax(powden_i)+1.e-2])
#axis([0.,1., 0., 0.35])

plt.subplot(223) #-------------------------
# plt.hold(True)
plt.grid(True)
txt="$P_{cl}$"+"$=$%3.3f" %(powtot_cl) +" $kW$"
plt.title(txt,verticalalignment='center',y=1.02)
plt.xlabel(r'$\rho$')
plt.ylabel('$Collisional$ $p_{cl}$')
if n_powden_cl>0:
    plt.plot(rho_bin_center,(powden_cl),'r',linewidth=linw)
    axis([0.,rhomax,0.,1.05*np.amax(powden_cl)+1.e-2])

plt.subplot(224) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.xlabel(r'$\rho$')
txt="$P_{total}$"+"$=$%3.3f" %(power_total) +" $kW$"
plt.title(txt,verticalalignment='center',y=1.02)
plt.plot(rho_bin_center,(powden),'r',linewidth=linw)
axis([0.,rhomax,0.,1.05*np.amax(powden)+1.e-2])

savefig('genray_profiles_power.png')
show() 




#--------------------------------------------------------------------------
fig1=plt.figure()   # Current profiles
plt.subplot(221) #-------------------------
title('$Current$ $Profiles$  $(A/cm^2)$',y=1.03)
# plt.hold(True)
plt.grid(True)
plt.ylabel('$<j_{||}>$')
plt.plot(rho_bin_center,s_cur_den_parallel,'r',linewidth=linw)
xlim( (0,rhomax) ) 
plt.subplot(222) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$<j.B>/B_0$',y=1.03)
plt.plot(rho_bin_center,s_cur_den_onetwo,'r',linewidth=linw)
xlim( (0,rhomax) ) 
plt.subplot(223) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.xlabel(r'$\rho$')
plt.ylabel('$<j_{||}>RB<R^{-2}>/(<B><1/R>)$')
#plt.title('$<j_{||}>RB<1/R^2>/(<B><1/R>)$',verticalalignment='center',y=1.03)
plt.plot(rho_bin_center,s_cur_den_toroidal,'r',linewidth=linw)
xlim( (0,rhomax) ) 
plt.subplot(224) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.xlabel(r'$\rho$')
plt.title('$<j_{||}>B_{pol}($'+r'$\theta$'+'$=0)/<B>$',verticalalignment='center',y=0.8)
plt.plot(rho_bin_center,s_cur_den_poloidal,'r',linewidth=linw)
xlim( (0,rhomax) ) 
savefig('genray_profiles_J.png')
show() 


fig7=plt.figure()   # Pe and J|| profiles
#----------------
plt.subplot(221) #------------------------- Pe
plt.ylabel('$p_e$ $(kW/m^3)$') 
# plt.hold(True)
plt.grid(True)
txt="$P_{e}$"+"$=$%3.3f" %(powtot_e/1e3) +" $MW$"
plt.title(txt,y=1.01)
plt.plot(rho_bin_center,(powden_e),'r',linewidth=linw*2)
axis([0.,rhomax,0.,1.05*np.amax(powden_e)+1.e-2])
#----------------
#'parallel_cur_total' [A]
parallel_cur_total=dat.variables['parallel_cur_total']
print(parallel_cur_total.long_name,  parallel_cur_total.getValue(),  parallel_cur_total.units)
parallel_cur_total=np.asscalar(parallel_cur_total[:])/1e3  # kA
print('parallel_cur_total[kA]= ', parallel_cur_total)
#'toroidal_cur_total' [A]
toroidal_cur_total=dat.variables['toroidal_cur_total']
print(toroidal_cur_total.long_name,  toroidal_cur_total.getValue(),  toroidal_cur_total.units)
toroidal_cur_total=np.asscalar(toroidal_cur_total[:])/1e3  # kA
print('toroidal_cur_total[kA]= ', toroidal_cur_total)
#'poloidal_cur_total' [A]
poloidal_cur_total=dat.variables['poloidal_cur_total']
print(poloidal_cur_total.long_name,  poloidal_cur_total.getValue(),  poloidal_cur_total.units)
poloidal_cur_total=np.asscalar(poloidal_cur_total[:])/1e3  # kA
print('poloidal_cur_total[kA]= ', poloidal_cur_total)
#'GA_tor_cur_total'   [A]
GA_tor_cur_total=dat.variables['GA_tor_cur_total']
print(GA_tor_cur_total.long_name,  GA_tor_cur_total.getValue(),  GA_tor_cur_total.units)
GA_tor_cur_total=np.asscalar(GA_tor_cur_total[:])/1e3  # kA
print('GA_tor_cur_total[kA]= ', GA_tor_cur_total)
plt.subplot(223) #------------------------- J||
txt="$I_{CD}$"+"$=$%3.6f" %(parallel_cur_total/1e3) +" $MA$"
title(txt,x=1.5)
# plt.hold(True)
plt.grid(True)
plt.ylabel('$<j_{||}>/P_e$  $(A/m^2/W)$')
j_P_min= 10*np.min(s_cur_den_parallel[:])/powtot_e
j_P_max= 10*np.max(s_cur_den_parallel[:])/powtot_e
j_P_min=min(j_P_min,0)
j_P_max=max(j_P_max,0)
dj= 0.1*(j_P_max-j_P_min) +1e-3
axis([0.,rhomax, j_P_min-dj, j_P_max+dj])
#if  (j_P_max>0 and j_P_max<0.01):
#    plt.yticks([0.0, 0.001, 0.002, 0.003, 0.004]) #, 0.005, 0.006, 0.007, 0.008])
plt.plot(rho_bin_center,10*s_cur_den_parallel/powtot_e,'r',linewidth=linw*2)
# s_cur_den_parallel is in A/cm^2, powtot_e is in kW,
# (A/cm^2)/kW = 10000*(A/m^2)/kW = 10*(A/m^2)/W
plt.xlabel(r'$\rho$')
savefig('genray_Pe_JII.png')
show()



fig2=plt.figure()
wpe_w_2= 806.2e5*sene[:,:]*q2m/f**2

plt.subplot(211) #-------------------------
# plt.hold(True)
plt.grid(True)
#plt.title(delpwr.long_name)
#plt.xlabel('$poloidal$ $distance$ $along$ $ray$  $ws$  $(cm)$')
plt.title('$power$ $in$ $ray$ $(MW)$',y=1.03)
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
# Plot with different thickness of lines, depending on ray
linw_mx= 6*linw # largest thickness
linw_mn= 0.5    # smallest
if Nrays>1:
    dlinw= (linw_mx-linw_mn)/(Nrays-1)
else:
    dlinw=0
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= nrayelt[i] # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],delpwr[i,0:Nm]/1e13,linewidth=linw_mx-i*dlinw)

plt.subplot(212) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$(\omega_{pe}/\omega)^2$  $along$ $ray$',y=0.8)
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
#print 'shape of wpe_w_2 ', wpe_w_2.shape
#print 'shape of sene ', sene.shape
y_mn= np.min(wpe_w_2) # second arg: Nm=points_along_ray (can be 0)
y_mx= np.max(wpe_w_2)
dy= 0.35*(y_mx-y_mn)
plt.ylim((y_mn,y_mx+dy))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    if Nm>0:
        #print 806.2e5*sene[0,0]*q2m/f**2, sene[0,0], q2m
        plt.plot(ws[i,0:Nm],806.2e5*sene[i,0:Nm]*q2m/f**2,color=col,linewidth=linw_mx-i*dlinw)
    
plt.xlabel('$pol.distance$ $along$ $ray$  $(cm)$')
plt.savefig('genray_rays_delpwr_lin.png') 
#plt.savefig('genray_rays_delpwr_lin.eps') 
show()


fig2=plt.figure() # Same as above, but log10 scale for delpwr
plt.subplot(211) #-------------------------
# plt.hold(True)
plt.grid(True)
#plt.title(delpwr.long_name)
#plt.xlabel('$poloidal$ $distance$ $along$ $ray$  $ws$  $(cm)$')
plt.title('$power$ $in$ $ray$ $(MW)$',y=1.03)
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
# Plot with different thickness of lines, depending on ray
linw_mx= 6*linw # largest thickness
linw_mn= 0.5    # smallest
p_tiny=1e-10
if Nrays>1:
    dlinw= (linw_mx-linw_mn)/(Nrays-1)
else:
    dlinw=0
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= nrayelt[i] # max number of points along a ray
    if Nm>0:
        semilogy(ws[i,0:Nm],(p_tiny+delpwr[i,0:Nm])/1e13,linewidth=linw_mx-i*dlinw)

plt.subplot(212) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$(\omega_{pe}/\omega)^2$  $along$ $ray$',y=0.8)
xticks(np.arange(ws_min, ws_max, step=ws_step) )
plt.xlim( (ws_min, ws_max) )
dy= 0.35*(y_mx-y_mn)
plt.ylim((y_mn,y_mx+dy))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    if Nm>0:
        #print 806.2e5*sene[0,0]*q2m/f**2, sene[0,0], q2m
        plt.plot(ws[i,0:Nm],806.2e5*sene[i,0:Nm]*q2m/f**2,color=col,linewidth=linw_mx-i*dlinw)
    
plt.xlabel('$pol.distance$ $along$ $ray$  $(cm)$')
plt.savefig('genray_rays_delpwr_log.png') 
#plt.savefig('genray_rays_delpwr_log.eps') 
show()


fig2=plt.figure() # Same as above, but delpwr vs rho(along ray)
plt.subplot(211) #-------------------------
# plt.hold(True)
plt.grid(True)
rho_min=np.min(spsi)
rho_min=min(0.,rho_min) # extend lower range to 0.
rho_max=np.max(spsi)
rho_max=max(1.,rho_max) # extend upper range to 1.
rho_step= ((rho_max*1.1-rho_min)/4)  # step for xticks() labeling ;
                     #    ()/4 means 4 labels only
delpwr_max=np.max(delpwr/1e13)
plt.title('$power$ $in$ $ray$ $(MW)$',y=1.03)
#xticks(np.arange(rho_min, rho_max, step=rho_step) )
plt.xlim( (rho_min, rho_max) )
plt.ylim( (delpwr_max/1e5, delpwr_max*2) )
# Plot with different thickness of lines, depending on ray
linw_mx= 6*linw # largest thickness
linw_mn= 0.5    # smallest
if Nrays>1:
    dlinw= (linw_mx-linw_mn)/(Nrays-1)
else:
    dlinw=0
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= nrayelt[i] # max number of points along a ray
    if Nm>0:
        semilogy(spsi[i,0:Nm],delpwr[i,0:Nm]/1e13,linewidth=linw_mx-i*dlinw)
        semilogy(spsi[i,0],delpwr[i,0]/1e13,'ko') # mark starting point

plt.subplot(212) #-------------------------
# plt.hold(True)
plt.grid(True)
plt.title('$(\omega_{pe}/\omega)^2$  $along$ $ray$',y=0.8)
#xticks(np.arange(rho_min, rho_max, step=rho_step) )
plt.xlim( (rho_min, rho_max) )
dy= 0.35*(y_mx-y_mn)
plt.ylim((y_mn,y_mx+dy))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i]) # max number of points along a ray
    if Nm>0:
        #print 806.2e5*sene[0,0]*q2m/f**2, sene[0,0], q2m
        plt.plot(spsi[i,0:Nm],806.2e5*sene[i,0:Nm]*q2m/f**2,color=col,linewidth=linw_mx-i*dlinw)
        plt.plot(spsi[i,0],806.2e5*sene[i,0]*q2m/f**2,'ko') # starting point
    
plt.xlabel('$rho$ $along$ $ray$')
plt.savefig('genray_rays_delpwr_vs_rho.png') 
show()




dat.close() # close genray.nc
         
elapsed_time = time.time() - e0
cpu_time = time.clock() - c0
print('elapsed and cpu time since start (sec.) =', elapsed_time, cpu_time)
print('FINISHED')
