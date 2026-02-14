# genray_c_plot.py'
# Plots genray.nc (output data file produced by GENRAY-C)
# Yuri Petrov   CompX   2011

from numpy import *
from mpl_toolkits.mplot3d import Axes3D

from pylab import *
from matplotlib import rc 
from matplotlib.pyplot import cm, figure, axes, plot, xlabel, ylabel,  \
     title, savefig, show

import netCDF4
print 'netCDF4:', netCDF4.__version__

from netCDF4 import Dataset

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import time
import pylab as pylab
#matplotlib.interactive(True) # no plots on screen
matplotlib.interactive(False) # with plots on screen
print 'matplotlib version:', matplotlib.__version__

# Detect which OS you are using; By default, it is 'windows'
ios='windows' #Will be changed to 'linux', if Python detects it, see below. 
# Some of functionality (such as using Latex) depends on OS
if sys.platform.startswith('win32'):
    print 'OS is windows'
    ios='windows'    
elif sys.platform.startswith('linux'):
    print 'OS is linux'
    ios='linux'
    #Render with externally installed LateX:
    matplotlib.rc('text', usetex = True) #ONLY FOR LINUX !does not work for windows
# Other possible values for sys.platform.startswith(): 
#     'darwin' for MacOS, 'cygwin' for Windows/Cygwin, 'aix' for AIX



# Specify filename with data (normally it is genray.nc, unless it was renamed):
filenm='genray.nc'
#filenm='WHAM2_modelb2_B0.853T_Rlim0.0964m_2.1.nc'
#filenm='WHAM2.6.nc'
#filenm='STEP_SPR14_helicon400.nc'
#filenm='STEP_SPR08_helicon400.nc'

#-----------------------------------------
iplot_cql=0 #1 #To add plots from CQL3D run, set iplot_cql=1, specify file_cql_krf
#file_cql_krf='gdt_hhfw_run8.3_krf001.nc'
#file_cql_krf='gdt_hhfw_run8.3_krf002.nc'
#-----------------------------------------

# Specify for plots:
rhomax=1.2 #2.0 #0. #3.0  # In the code, rhomax is defined as max of rho 
            # over ALL equilibrium grid - it could be very large.
            # If the power deposition happened within plasma (rho<0)
            # the power profiles may look "compressed".
            # For clarity, set the limit of rho for plots like Pe(rho).
            # If the whole rho range is desired, set rhomax to 0.
            # it will be determined automatically.

Zmin=-80. #-75. #-190. #[cm] Z Limits for plots. Set to 0 if you want it done automatically
Zmax=+80. #+75. #+190. #[cm] Z Limits for plots. Set to 0 if you want it done automatically          
            
fnt  = 10 #20 #18 #10     # font size for axis numbers (see 'param=' below) 
linw = 1.0    # LineWidth for contour plots
Ncont= 10     # Number of contour levels for PSI (pol.flux) or density
#add labels in contour lines of omega/omega_c over (R,Z) plane:
contour_labels=2 # 0 - do not add labels; 1 - add manually;  2 - automatically
isp_wc=1  #  Species number for which you want to plot omega/omega_c res.layers
          #  1- electrons, 2 or higher - ions
# levels for omega/omega_c  to be plotted:
Dmin=1
Dmax=3 #6 #15
nres=1 #2 # Specify resonance harmonic number such as omega ~ nres*omega_ci
#(for plots genray_rays_VIIres_p.png)
levels_wwc=np.arange(Dmin,Dmax,1) # from Dmin to Dmax

arr_len=3.  # Specify arrow length for plots of (Nr,Nz) at start point

i_cps=0 #Added [2018-10], for new feature "cross_polarization_scat"
        #To skip related plots, set i_cps to 0
        
i_plot_powerRZ= -1 # Choices: 0, -1, +1, 2
# The 2D plots of absorbed power in (Z,R) plane are based on data accumulated 
# over R=sqrt(X^2 +Y^2) radial grid. If ray is propagating 
# at negative X and/or negative Y, the value of R is still positive.
# In general, the path of rays in (Z,X) plane may not coincide 
# with power pattern in (Z,R) grid because of such mapping.
# For convenience, the power deposition P(Z,R) can be also plotted
# as mirror reflection over the Z axis, i.e. as P(Z,-R),
# or as the combination of the two, i.e. P(Z,R) & P(Z,-R).
# i_plot_powerRZ= 0 ==> Do not make such plots at all.
# i_plot_powerRZ=+1 ==> Plot P(Z, R) (where R is .ge.0).
# i_plot_powerRZ=-1 ==> Plot mirror image P(Z,-R) 
# i_plot_powerRZ= 2 ==> Plot the combination P(Z,R)&P(Z,-R) 

# For plots of Dispersion Function vs x, showing Hot plasma roots
# and also plots of cold Nperp vs x:
plot_Nperp_prof=1   # 0- do not plot; 1-plot (takes some time)
dispers_log_scale=0 # 1 for plotting as log10(Nperp), 0 - linear scale
hot_roots=1 # 1 - plot dispersion from id=6;  0 - do not plot
            # Plotting hot roots is time-consuming, so if not needed, set to 0.
Nper_cold_ylim=300. #1.e5 # Set upper limit for plot of cold roots 
                    #(if set to 0, the upper limit will be 1e5)
Nper_hot_ylim=300. #1.e5  # Set upper limit for plot of hot roots
                    #(if set to 0, the upper limit will be 1e5)
Nper_im_ylim=0. #40.0   # Set upper limit for plot of imaginary N_im refr. index
                    #(if set to 0, the upper limit will be found automatically)
#-------------------


# Constants -----------------------------------
pi=3.14159265358979
clight= 2.99792458e10   # speed of light [cm/s]
e     = 4.8032e-10      # e-charge [cgs]
e_mass= 9.10938291e-28  # e-mass   [gram]
p_mass= 1.67262158e-24  # proton mass    [gram]
ergtkev=1.6022e-09
mp=p_mass/e_mass

#===================  for PLOTS =============================================
#set fonts and line thicknesses
params = {
    'axes.linewidth': linw,
    'lines.linewidth': linw,
    'axes.labelsize': fnt+4,
    'font.size': fnt+4,
    'legend.fontsize': fnt,
    'xtick.labelsize':fnt,
    'ytick.labelsize':fnt,
    'font.weight'  : 'regular'
}
mpl.rcParams.update(params)
#rc.defaults() #to restore defaults
mpl.rcParams['font.size']=fnt+2  # set font size for text in mesh-plots

e0 = time.time()  # elapsed time since the epoch
c0 = time.clock() # total cpu time spent in the script so far



# Open the genray netcdf file, read only
dat= Dataset(filenm, 'r', format='NETCDF4')
#==============================================================================
print dat.file_format # print which format was used for the genray.nc file
print '========================================'
print 'The genray file, ',filenm,', contains:'
print '========================================'
print "The global attributes: ",dat.dimensions.keys()        
print "File contains variables: ",dat.variables.keys()
print '========================================'

print '---------------------------------------GENERAL PROFILES'
# some eqdsk data:
eqdsk_x=dat.variables['eqdsk_x']
print 'eqdsk_x is: ', eqdsk_x.long_name, eqdsk_x.shape
eqdsk_y=dat.variables['eqdsk_y']
print 'eqdsk_y is: ', eqdsk_y.long_name, eqdsk_y.shape
eqdsk_r=dat.variables['eqdsk_r']
print 'eqdsk_r is: ', eqdsk_r.long_name, eqdsk_r.shape
eqdsk_z=dat.variables['eqdsk_z']
print 'eqdsk_z is: ', eqdsk_z.long_name, eqdsk_z.shape
eqdsk_psi=dat.variables['eqdsk_psi']  
# THIS IS peqd array in Genray-C, always ascending so that psilim>psimag.
# The direction of Bz can be any.
print 'eqdsk_psi is: ', eqdsk_psi.long_name, eqdsk_psi.shape
print 'min(eqdsk_psi),max(eqdsk_psi)=',np.min(eqdsk_psi),np.max(eqdsk_psi)
# These two values (psimag,psilim) are added on 07-28-2016 into genray.nc:
psimag=dat.variables['psimag'].getValue()  #getValue() for scalar
psimag=np.asscalar(psimag)
psilim=dat.variables['psilim'].getValue()  #getValue() for scalar
psilim=np.asscalar(psilim)
print 'psimag,psilim=',psimag,psilim

eqdsk_x=np.asarray(eqdsk_x)
eqdsk_y=np.asarray(eqdsk_y)
eqdsk_z=np.asarray(eqdsk_z)
eqdsk_r=np.asarray(eqdsk_r)
eqdsk_psi=np.asarray(eqdsk_psi)
#print 'eqdsk_r=',eqdsk_r
#print 'eqdsk_x=',eqdsk_x
#print 'eqdsk_y=',eqdsk_y
nxeqd=np.size(eqdsk_x)
nzeqd=np.size(eqdsk_z)
nreqd=np.size(eqdsk_r)
print 'nxeqd,nzeqd,nreqd=',nxeqd,nzeqd,nreqd
print 'eqdsk_r[0], eqdsk_r[nreqd-1], eqdsk_z[0],eqdsk_z[nzeqd-1]=',\
eqdsk_r[0], eqdsk_r[nreqd-1], eqdsk_z[0],eqdsk_z[nzeqd-1]

model_b=dat.variables['model_b'].getValue()  #getValue() for scalar
model_b=np.asscalar(model_b)
print 'model_b=',model_b
model_rho_dens=dat.variables['model_rho_dens'].getValue()  #getValue() for scalar
model_rho_dens=np.asscalar(model_rho_dens)
if model_b==4: # FRC
    rs_frc=dat.variables['rs_frc'].getValue()  #getValue() for scalar
    rs_frc= np.asscalar(rs_frc)*100 # cm FRC separatrix
    r0_frc= rs_frc/sqrt(2)  #   cm
        
rscan=dat.variables['xscan'] # it was 'rscan' before 09-30-2014
rscan=rscan[:]*100  # m->cm
rhoscan=dat.variables['rhoscan']  #  rho(xscan)
wcw=dat.variables['wcwscan']  # wce/w as a func of rscan
wpw=dat.variables['wpwscan']  # wpe/w as a func of rscan
wuw=dat.variables['wuwscan']  # wuh/w as a func of rscan (Upper_hybrid/omega)
wlw=dat.variables['wlwscan']  # wlh/w as a func of rscan (Lower_hybrid/omega)[2022]
wcw=np.asarray(wcw)
wpw=np.asarray(wpw)
wuw=np.asarray(wuw)
wlw=np.asarray(wlw)

# Values of Y and Z at which the scan is done:
yscan=dat.variables['y_save_disp'] # a single value in the code
zscan=dat.variables['z_save_disp'] # a single value [m]
yscan=np.asscalar(yscan[:])*100 # cm
zscan=np.asscalar(zscan[:])*100 # cm
rscan_size=np.size(rscan)


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
        n_wall=r_wall.size
        print 'r_wall is: ', r_wall.long_name, r_wall.shape
        print 'z_wall is: ', z_wall.long_name, z_wall.shape  
        print 'n_wall=',n_wall
        r_wall=np.asarray(r_wall)*100
        z_wall=np.asarray(z_wall)*100
finally:
    print '----------------------------------------'

if plot_Nperp_prof==1:
    try:
        try:
            Nperp2m=dat.variables['Nperp2m']
        except:
            print('No data on Nperp2(r) in genray.nc')
            plot_Nperp_prof=0 # reset to 0, meaning: no plots of Disp.Func.
        else:
            Npara=dat.variables['Npara']
            Npera=dat.variables['Npera']
            Nperp2m=dat.variables['Nperp2m']
            Nperp2p=dat.variables['Nperp2p']
            Nperp_im= dat.variables['Nperp_im']
            Nperp2hot=dat.variables['Nperp2hot']
            ddd=dat.variables['ddd']
            inpar=np.size(Nperp2m,0)
            inper=np.size(ddd,1)
            print 'inpar=',inpar
            print 'rscan_size=',rscan_size
            print 'inper=',inper
    finally:
        print '----------------------------------------'
        
rscan_max=np.amax(rscan[:])
rscan_min=np.amin(rscan[:])
print 'rscan_mx=', rscan_max


# Added [2018-10], for new feature "cross_polarization_scat" :
if i_cps==1:
    i_cps=dat.variables['i_cps'].getValue()  
    i_cps= np.asscalar(i_cps) # =1 if CPS is enabled
    nray_ant=dat.variables['nray_ant'].getValue()  
    nray_ant= np.asscalar(nray_ant) # Primary rays launched from antenna
    nray_add=dat.variables['nray_add'].getValue()  
    nray_add= np.asscalar(nray_add) # Offspring rays scattered from primary ray(s)

    # N. of start pts (for offspring rays) on each primary ray:
    nray_start_cps=dat.variables['nray_start_cps'].getValue()  
    nray_start_cps= np.asscalar(nray_start_cps) 

    print 'i_cps=',i_cps
    print 'nray_ant, nray_add=',nray_ant,nray_add
    print 'nray_start_cps=',nray_start_cps

    detector_cps_r= dat.variables['detector_cps_r']
    print 'detector_cps_r',detector_cps_r.long_name,detector_cps_r[:],\
    detector_cps_r.units
    detector_cps_r= detector_cps_r.getValue()  #getValue() for scalar
    detector_cps_r= np.asscalar(detector_cps_r) # m

    detector_cps_z= dat.variables['detector_cps_z']
    print 'detector_cps_z',detector_cps_z.long_name,detector_cps_z[:],\
    detector_cps_z.units
    detector_cps_z= detector_cps_z.getValue()  #getValue() for scalar
    detector_cps_z= np.asscalar(detector_cps_z) # m

    detector_cps_phi= dat.variables['detector_cps_phi']
    print 'detector_cps_phi',detector_cps_phi.long_name,detector_cps_phi[:],\
    detector_cps_phi.units
    detector_cps_phi= detector_cps_phi.getValue()  #getValue() for scalar
    detector_cps_phi= np.asscalar(detector_cps_phi) # degrees

    detector_cps_diam= dat.variables['detector_cps_diam']
    print 'detector_cps_diam',detector_cps_diam.long_name,detector_cps_diam[:],\
    detector_cps_diam.units
    detector_cps_diam= detector_cps_diam.getValue()  #getValue() for scalar
    detector_cps_diam= np.asscalar(detector_cps_diam) # m

    # Evaluate (X,Y,Z) coords of detector position:
    detector_cps_x= detector_cps_r*cos(detector_cps_phi*pi/180) # m
    detector_cps_y= detector_cps_r*sin(detector_cps_phi*pi/180) # m
    print 'detector_cps_x=',detector_cps_x, ' m'
    print 'detector_cps_y=',detector_cps_y, ' m'
    print 'detector_cps_z=',detector_cps_z, ' m'

    i_hit_detector_cps= dat.variables['i_hit_detector_cps']
    print 'i_hit_detector_cps', i_hit_detector_cps.long_name
    i_hit_detector_cps= np.asarray(i_hit_detector_cps)
    print 'i_hit_detector_cps= ', i_hit_detector_cps
    # Note: the indexing of i_hit_detector_cps[] array goes over ALL rays,
    # starting with primary rays (i_hit_detector_cps[]=0 for primary rays)

    i_bad_initial_conditions= dat.variables['i_bad_initial_conditions']
    print 'i_bad_initial_conditions', i_bad_initial_conditions.long_name
    i_bad_initial_conditions= np.asarray(i_bad_initial_conditions)
    print 'i_bad_initial_conditions= ', i_bad_initial_conditions
    # Note: the indexing of i_bad_initial_conditions[] array goes over ALL rays,
    # starting with primary rays 
    # (i_bad_initial_conditions[]=0 for rays that were started, 1 otherwise)
# Done [2018-10] new feature "cross_polarization_scat" 


freqcy=dat.variables['freqcy']
print 'freqcy =',freqcy.long_name, freqcy[:], freqcy.units
freqcy=freqcy.getValue()  #getValue() for scalar
f=freqcy  #Hz
omega=f*2*pi

mass=dat.variables['dmas']   #  dmas=Mass(species)/Mass_electron
print 'mass=', mass.long_name, mass[:], mass.units
# To find the thermal speed:
#	 vi=sqrt(2T/m),    vi in (cm/sec),T(keV)
#	 vi=1.87d9*dsqrt(T(i)/dmas(i))


charge=dat.variables['charge']
print 'charge=', charge.long_name, charge[:], charge.units

rho_bin_center=dat.variables['rho_bin_center']
print 'rho_bin_center is: ', rho_bin_center.long_name, rho_bin_center.shape
#print rho_bin_center[:]

rho_bin=dat.variables['rho_bin']
print 'rho_bin is: ', rho_bin.long_name, rho_bin.shape


try:
    try:
        binvol=dat.variables['binvol']
    except:
        print('No data on binvol in genray.nc')
        i_binvol=0
    else:
        i_binvol=1
        binvol=dat.variables['binvol'] #
        print 'binvol is: ', binvol.long_name, binvol.shape
finally:
    print '---'
    
densprof=dat.variables['densprof']
print 'densprof is: ', densprof.long_name, densprof.shape

temprof=dat.variables['temprof']
print 'temprof is: ', temprof.long_name, temprof.shape

Nsp=temprof[:,0].size  # number of species
print 'Number of species: Nsp=',Nsp
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
    if (np.absolute(Zs-4.)<0.1) & (np.absolute(msme-12*mp)<100):
        isp_name.append('C+4')
    if (np.absolute(Zs-5.)<0.1) & (np.absolute(msme-12*mp)<100):
        isp_name.append('C+5')
    if (np.absolute(Zs-6.)<0.1) & (np.absolute(msme-12*mp)<100):
        isp_name.append('C+6')
    if (np.absolute(Zs-54.)<0.1) & (np.absolute(msme-130*mp)<1000):
        isp_name.append('Xe+54')
    # OTHER NAMES CAN BE ADDED HERE.
    print 'isp=',isp,' m=',mass[isp],' Z=',charge[isp],' isp_name=',isp_name[isp]
    #txt1=r"$m/m_e =$"+r"$%1.0f$" %(mass[i])
    #txt2=r"$q/e =$"+r"$%1.0f$" %(charge[i])
print isp_name
    

w_dens_vs_x_nc=dat.variables['w_dens_vs_x_nc']
print 'w_dens_vs_x_nc is: ', w_dens_vs_x_nc.long_name, w_dens_vs_x_nc.shape
w_temp_vs_x_nc=dat.variables['w_temp_vs_x_nc']
print 'w_temp_vs_x_nc is: ', w_temp_vs_x_nc.long_name, w_temp_vs_x_nc.shape
w_x_densprof_nc=dat.variables['w_x_densprof_nc']
print 'w_x_densprof_nc is: ', w_x_densprof_nc.long_name, w_x_densprof_nc.shape
w_bmod_vs_x_nc=dat.variables['w_bmod_vs_x_nc']
print 'w_bmod_vs_x_nc is: ', w_bmod_vs_x_nc.long_name, w_bmod_vs_x_nc.shape

w_dens_vs_y_nc=dat.variables['w_dens_vs_y_nc']
print 'w_dens_vs_y_nc is: ', w_dens_vs_y_nc.long_name, w_dens_vs_y_nc.shape
w_temp_vs_y_nc=dat.variables['w_temp_vs_y_nc']
print 'w_temp_vs_y_nc is: ', w_temp_vs_y_nc.long_name, w_temp_vs_y_nc.shape
w_y_densprof_nc=dat.variables['w_y_densprof_nc']
print 'w_y_densprof_nc is: ', w_y_densprof_nc.long_name, w_y_densprof_nc.shape
w_bmod_vs_y_nc=dat.variables['w_bmod_vs_y_nc']
print 'w_bmod_vs_y_nc is: ', w_bmod_vs_y_nc.long_name, w_bmod_vs_y_nc.shape

# Added 01-12-2016
bmodprofxz=dat.variables['bmodprofxz']
D=bmodprofxz # [Tesla]  B(X,Z)
Dmin=np.min(D)
Dmax=np.max(D) #*0.01
print 'Bmin[T]=',Dmin, '  Bmax[T]=',Dmax
isp=0 # electrons, by default
if Nsp>1:
    isp=isp_wc-1 # ion (isp_wc is the genray-c index; isp is the python index)
if isp>-1: # if isp<0, no plots for w/wc
    mime= mass[isp]/mass[0]      # m_s/m_e ratio
    Z_s= charge[isp]/charge[0]    # Z_s/e
    w_wc= f/(28e9*D[:]*Z_s/mime)  # omega/omega_c
    print 'For species#', isp_wc, '(', isp_name[isp], ') in GENRAY-C, we have:'
    print 'w/wc min/max =', f/(28e9*Dmax*Z_s/mime), f/(28e9*Dmin*Z_s/mime)
    print 'But plotted only these levels:', levels_wwc


# Density profiles on xy, xz and yz-grids
densprofxy=dat.variables['densprofxy']
print 'densprofxy is: ', densprofxy.long_name, densprofxy.shape
densprofxz=dat.variables['densprofxz']
print 'densprofxz is: ', densprofxz.long_name, densprofxz.shape
densprofyz=dat.variables['densprofyz']
print 'densprofyz is: ', densprofyz.long_name, densprofyz.shape

# O-X modes Transmission data: Only present for i_ox=2:
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
finally:
    print '----------------------------------------'
    

# Define min/max for major radius range:
if n_wall==0:  # wall is not defined; get limits from eqdsk:
    Rmax=np.amax(eqdsk_r)*100
    Rmin=np.amin(eqdsk_r)*100 # meters -> cm
    Rmax=np.asscalar(Rmax)
    Rmin=np.asscalar(Rmin)
else:
    Rmax=np.amax(r_wall)
    Rmin=np.amin(r_wall)
    Rmax=np.asscalar(Rmax)
    Rmin=np.asscalar(Rmin)
    if Zmin==0: 
        Zmin=np.amin(z_wall)
        Zmin=np.asscalar(Zmin)
    if Zmax==0: 
        Zmax=np.amax(z_wall)
        Zmax=np.asscalar(Zmax)
        
Rmax=max(Rmax,abs(Rmin))
Rmin=-Rmax
    
print 'Rmin,Rmax=', Rmin,Rmax , 'Zmin,Zmax=', Zmin,Zmax 


#----------------- For processing data from CQL3D run ------------------------
if iplot_cql==1:
    # Open file_cql_krf netcdf file, read only
    dat_cql= Dataset(file_cql_krf, 'r', format='NETCDF4')
    print dat_cql.file_format # print which format was used
    print '========================================'
    print 'The CQL3D file, ',file_cql_krf,', contains:'
    print '========================================'
    print "The global attributes: ",dat_cql.dimensions.keys()        
    print "File contains variables: ",dat_cql.variables.keys()
    print '========================================'
    delpwr_cql=dat_cql.variables['delpwr']
    print 'delpwr_cql is: ', delpwr_cql.long_name, delpwr_cql.shape
    delpwr_cql= delpwr_cql[:,:]/1.e10  # Convert to kW
    #--- Linear damping:     
    urfpwrl_cql=dat_cql.variables['urfpwrl']
    print 'urfpwrl_cql is: ', urfpwrl_cql.long_name, urfpwrl_cql.shape
    urfpwrl_cql= urfpwrl_cql[:,:] 
    # Note: dprfl= urfpwrl_cql[i,irayelt]*delpwr_cql[i,irayelt] is the 
    #              Lin.absorped power at given ray element.
    #--- Collisional damping:     
    urfpwrc_cql=dat_cql.variables['urfpwrc']
    print 'urfpwrc_cql is: ', urfpwrc_cql.long_name, urfpwrc_cql.shape
    urfpwrc_cql= urfpwrc_cql[:,:] 
    # Note: dprfc= urfpwrc_cql[i,irayelt]*delpwr_cql[i,irayelt] is the 
    #              Coll.absorped power at given ray element.
    # E-wave-field along rays
    cwexde_cql=dat_cql.variables['cwexde']
    print 'cwexde_cql is: ', cwexde_cql.long_name, cwexde_cql.shape
    cweyde_cql=dat_cql.variables['cweyde']
    print 'cweyde_cql is: ', cweyde_cql.long_name, cweyde_cql.shape
    cwezde_cql=dat_cql.variables['cwezde']
    print 'cwezde_cql is: ', cwezde_cql.long_name, cwezde_cql.shape
    # RAY trajectories
    wr_cql=dat_cql.variables['wr']
    print 'wr_cql is: ', wr_cql.long_name, wr_cql.shape
    wz_cql=dat_cql.variables['wz']
    print 'wz_cql is: ', wz_cql.long_name, wz_cql.shape
    ws_cql=dat_cql.variables['ws']
    print 'ws_cql is: ', ws_cql.long_name, ws_cql.shape
    # Number of rays
    Nrays_cql=dat_cql.variables['nray'].getValue() 
    Nrays_cql=np.asscalar(Nrays_cql)
    print 'Number of rays: Nrays_cql=',Nrays_cql
    # Number of elements for each ray:
    nrayelt_cql=dat_cql.variables['nrayelt']
    nrayelt_cql=np.asarray(nrayelt_cql)
    print 'nrayelt_cql= ', nrayelt_cql
    #------------------------------------
    # Setup Z-grid for accumulation of absorbed power, etc.
    z_cql_min=np.min(eqdsk_z)*100 # meters -> cm
    z_cql_max=np.max(eqdsk_z)*100 # meters -> cm
    # Or maybe simply use Zmin;Zmax
    z_cql_min=Zmin
    z_cql_max=Zmax
    nzgrid_cql=200  # Hardwired, for now
    print 'z_cql_min, z_cql_max[cm]=', z_cql_min, z_cql_max
    dzgrid_cql= (z_cql_max-z_cql_min)/(nzgrid_cql-1.)
    zgrid_cql=np.linspace(z_cql_min, z_cql_max, nzgrid_cql, endpoint=True)
    dzgrid_cql= zgrid_cql[1]-zgrid_cql[0]
    #print zgrid_cql.size, dzgrid_cql, dzgrid_cql*(nzgrid_cql-1)
    #-----------------------------------
    prf_accum_cql=zeros((nzgrid_cql,Nrays_cql)) # initialize
    # The prf_accum_cql array will hold the RF absorbed power.
    prf_accum_allrays=zeros((nzgrid_cql)) #summed over rays for given krf wave type
    # Also, for Linear absorption:
    prfl_accum_cql=zeros((nzgrid_cql,Nrays_cql)) # initialize
    prfl_accum_allrays=zeros((nzgrid_cql))
    tot_absorbed=0.0 # For verification only: total absorbed power
    # Eplus
    Ep2=zeros((nzgrid_cql)) # |E+/E|^2 to be accumulated; from all rays
    Em2=zeros((nzgrid_cql)) # |E_/E|^2 to be accumulated; from all rays
    for i in range(0,Nrays_cql,1):  # i goes from 0 to Nrays_cql-1
        Nm= np.asscalar(nrayelt_cql[i])-1 # max number of points along a ray
        # To remove the last point: set Nm= np.asscalar(nrayelt_cql[i])-2
        tot_absorbed= tot_absorbed + (delpwr_cql[i,0]-delpwr_cql[i,Nm])
        # The above includes Linear and Coll damping (if used in run)
        for irayelt in range(1,Nm,1): # step along ray
            # Note that we start from second element (1) , because we form 
            # the diffrence of delpwr between irayelt and irayelt-1 elements
            zrayelt= wz_cql[i,irayelt] #[cm] Z at given ray element
            #--- Linear damping:
            dprfl= urfpwrl_cql[i,irayelt]*delpwr_cql[i,irayelt] # kW absorbed
            #--- Coll. damping:
            dprfc= urfpwrc_cql[i,irayelt]*delpwr_cql[i,irayelt] # kW absorbed
            #--- Power absorbed at ray element (including Lin.power, if any):
            dprf= -(delpwr_cql[i,irayelt]-delpwr_cql[i,irayelt-1]) # kW absorbed
            # dprf is power absorbed at ray element (with accuracy of one step).
            #--- Subtract Lin. and Coll. damped power; To be plotted separately:
            dprf= dprf-dprfl-dprfc
            # Find where ray element positioned w.r.to zgrid_cql :
            iz= np.floor( (zrayelt-z_cql_min)/dzgrid_cql ) #Based on uniform grid
            iz=max(iz,0) #In case ray is out of zgrid_cql, set lower iz=0
            iz=min(iz,nzgrid_cql-2) # set upper iz=nzgrid_cql-2
            #Note: the largest element of zgrid_cql is zgrid_cql[nzgrid_cql-1]
            # With the above procedure, 
            #  zrayelt is in  [zgrid_cql[iz] ; zgrid_cql[iz+1])  range
            if zrayelt<=zgrid_cql[iz]:
                print zgrid_cql[iz],zgrid_cql[iz+1],zrayelt
            if zrayelt>=zgrid_cql[iz+1]:
                print zgrid_cql[iz],zgrid_cql[iz+1],zrayelt
            # Weight factors for linear partition between iz and iz+1 nodes:
            wght_l= (zgrid_cql[iz+1]-zrayelt)/dzgrid_cql
            wght_r= 1.0-wght_l
            prf_accum_cql[iz,i]=  prf_accum_cql[iz,i]  +wght_l*dprf
            prf_accum_cql[iz+1,i]=prf_accum_cql[iz+1,i]+wght_r*dprf
            #The power from ray "i", absorbed at element "irayelt" is added 
            # to grid nodes iz and iz+1
            prfl_accum_cql[iz,i]=  prfl_accum_cql[iz,i]  +wght_l*dprfl
            prfl_accum_cql[iz+1,i]=prfl_accum_cql[iz+1,i]+wght_r*dprfl
            #----- Now Erf field at given ray element
            Ex_real=cwexde_cql[0,i,irayelt] # Real(Ex)/|E|
            Ex_imag=cwexde_cql[1,i,irayelt] # Imag(Ex)/|E|
            Ey_real=cweyde_cql[0,i,irayelt] # Real(Ey)/|E|
            Ey_imag=cweyde_cql[1,i,irayelt] # Imag(Ey)/|E|
            # Form Eplus as (Ex+i*Ey)/sqrt(2.)   
            Eplus_real=  (Ex_real-Ey_imag)/sqrt(2.)
            Eplus_imag=  (Ex_imag+Ey_real)/sqrt(2.)
            Eplus2= Eplus_real**2 + Eplus_imag**2 #== |E+/E|^2
            # Form Eminus as (Ex-i*Ey)/sqrt(2.)  
            Eminus_real= (Ex_real+Ey_imag)/sqrt(2.)
            Eminus_imag= (Ex_imag-Ey_real)/sqrt(2.)
            Eminus2= Eminus_real**2 + Eminus_imag**2 #== |E_/E|^2
            Ep2[iz]=   Ep2[iz]   +wght_l*Eplus2 # Added
            Ep2[iz+1]= Ep2[iz+1] +wght_r*Eplus2
            Em2[iz]=   Em2[iz]   +wght_l*Eminus2 # Added
            Em2[iz+1]= Em2[iz+1] +wght_r*Eminus2

    prf_accum_allrays= np.sum(prf_accum_cql,axis=1)
    print np.sum(prf_accum_allrays)
    print 'np.sum(prf_accum_cql) [kW]=', np.sum(prf_accum_cql)
    prfl_accum_allrays= np.sum(prfl_accum_cql,axis=1)
    print np.sum(prfl_accum_allrays)
    print 'np.sum(prfl_accum_cql) [kW]=', np.sum(prfl_accum_cql)
    print 'tot_absorbed, from delpwr [kW]=', tot_absorbed # Same as sum(prf_accum_allrays), to 5th digit
    
    #--------------------------------------------------------------------------
    fig1=plt.figure() # RF power deposited to general species. 
    #   if Linear and/or Coll. damping was used in CQL3D run, 
    #   they are NOT included in this plot (Lin. damping is plotted separately)
    ax = plt.subplot(1, 1, 1)
    plt.hold(True)
    txt1= "$[kW]=$"+r"$%1.3f$" %(np.sum(prf_accum_allrays)) 
    plt.title(file_cql_krf+':'+' $Absorbed$ $Power$'+txt1,fontsize=18,y=1.03)
    plt.xlabel('$Z$ $(cm)$  $[uniform$ $grid$ $in$ $Z]$' , fontsize=18)
    plt.ylabel('$Color:$ $Each$ $ray,$ $Black:Sum$ $ (kW) $' ,fontsize=18)
    plt.grid(True)
    for i in range(0,Nrays_cql,1):  # i goes from 0 to Nrays_cql-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        #Nm= np.asscalar(nrayelt_cql[i])-1 # max number of points along a ray
        linwi=4*linw-0.6*i
        linwi=max(1,linwi)
        plt.plot(zgrid_cql,prf_accum_cql[:,i],color=col,linewidth=linwi)
    #From All rays:
    plt.plot(zgrid_cql,prf_accum_allrays,color='k',linewidth=4)
    plt.axis([z_cql_min, z_cql_max, 0., 1.05*np.max(prf_accum_allrays)] )
    plt.savefig('cql_rays_pwr_vs_z.png') 
    show()
    #.............................................................
    fig2=plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.hold(True)
    txt1= "$[kW]=$"+r"$%1.3f$" %(np.sum(prfl_accum_allrays)) 
    plt.title(file_cql_krf+':'+' $Lin$ $Power$'+txt1,fontsize=18,y=1.03)
    plt.xlabel('$Z$ $(cm)$  $[uniform$ $grid$ $in$ $Z]$' , fontsize=18)
    plt.ylabel('$Color:$ $Each$ $ray,$ $Black:Sum$ $ (kW) $' ,fontsize=18)
    plt.grid(True)
    for i in range(0,Nrays_cql,1):  # i goes from 0 to Nrays_cql-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        #Nm= np.asscalar(nrayelt_cql[i])-1 # max number of points along a ray
        linwi=4*linw-0.6*i
        linwi=max(1,linwi)
        plt.plot(zgrid_cql,prfl_accum_cql[:,i],color=col,linewidth=linwi)
    #From All rays:
    plt.plot(zgrid_cql,prfl_accum_allrays,color='k',linewidth=4)
    plt.axis([z_cql_min, z_cql_max, 0., 1.05*np.max(prfl_accum_allrays)] )
    plt.savefig('cql_rays_pwrl_vs_z.png') 
    show()
    #...................................
    fig3=plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.hold(True)
    plt.title(file_cql_krf+':'+' $Power$ $in$ $ray$ $channel$',fontsize=18,y=1.03)
    plt.xlabel('$distance$ $along$ $ray$  $(cm)$' , fontsize=18)
    plt.ylabel('$ delpwr $  $ (kW) $' ,fontsize=18)
    ws_max=0.
    plt.grid(True)
    for i in range(0,Nrays_cql,1):  # i goes from 0 to Nrays_cql-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        Nm= np.asscalar(nrayelt_cql[i])-1 # max number of points along a ray
        ws_max1=np.amax(ws_cql[:,:])
        ws_max=np.amax([ws_max1,ws_max])
        linwi=4*linw-0.6*i
        linwi=max(1,linwi)
        plt.plot(ws_cql[i,0:Nm],delpwr_cql[i,0:Nm],color=col,linewidth=linwi)
    plt.axis([0.,1.05*ws_max,0.,1.05*np.amax(delpwr_cql[:,0])] )
    plt.savefig('cql_rays_delpwr.png') 
    show()

    #--------------------------------------------------------------------------
    fig4=plt.figure() # Eplus accumulated over Z-grid. 
    ax = plt.subplot(1, 1, 1)
    plt.hold(True)
    txt1= "$[kW]=$"+r"$%1.3f$" %(np.sum(prf_accum_allrays)) 
    plt.title(file_cql_krf+':'+' $|E_{+/-}/E|^2$',fontsize=18,y=1.03)
    plt.xlabel('$Z$ $(cm)$  $[uniform$ $grid$ $in$ $Z]$' , fontsize=18)
    plt.ylabel('$(a.u.)$  $Red:|E_{+}/E|^2;$ $Blue:|E_{-}/E|^2$' ,fontsize=18)
    plt.grid(True)
    #From All rays:
    plt.plot(zgrid_cql,Ep2,color='r',linewidth=1) # |E_{+}/E|^2: red
    plt.plot(zgrid_cql,Em2,color='b',linewidth=1) # |E_{-}/E|^2: blue
    Ep2_max= np.max(Ep2)
    Em2_max= np.max(Em2)
    plt.axis([z_cql_min, z_cql_max, 0., 1.05*max(Ep2_max,Em2_max)] )
    plt.savefig('cql_rays_Erf_vs_z.png') 
    show()
    
    stop #Stops here in case of iplot_cql==1; Comment, if you wish to continue




#================= genray.nc PLOTS ======================================
fig0=plt.figure()  # B(r=0,z) vs Z coord.
ix0= int(nxeqd/2) #530
x00= eqdsk_x[ix0]
Baxial=bmodprofxz[:,ix0]
f=freqcy
plt.subplot(211) #------------------------- B(Z)
plt.hold(True)
plt.grid(True)
plt.minorticks_on() # To add minor ticks
plt.tick_params(which='both',  width=1)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4, color='k')
plt.title('$z-scan$  $at$ $y=0,$ $x=$' +"%5.3f" %(x00) +'$cm$',y=1.03)
plt.ylabel('$B(z)$  $(T)$' )
#plt.xlabel('$ Z $   $ (cm) $')
plt.plot(eqdsk_z[:]*100,Baxial[:],'b',linewidth=1.5)
#--------------
plt.subplot(212) #------------------------- omega/omega_c vs Z
plt.hold(True)
plt.grid(True)
plt.minorticks_on() # To add minor ticks
plt.tick_params(which='both',  width=1)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4, color='k')
if Nsp>1:
    isp=isp_wc-1 # ion (isp_wc is the genray-c index; isp is the python index)
    mime= mass[isp]/mass[0]      # m_s/m_e ratio
    Z_s= charge[isp]/charge[0]    # Z_s/e
    plt.ylabel('$\omega/\omega_{c}$ ['+isp_name[isp]+']')
else:  # electrons only: 
    isp=0 # electrons
    mime=1.0
    Z_s=1.0
    plt.ylabel('$\omega/\omega_{ce}$')
plt.xlabel('$z$   $ (cm) $')
plt.plot(eqdsk_z[:]*100,f/(28e9*Baxial[:]*Z_s/mime),'b',linewidth=1.5)       
#--------------
savefig('genray_B_vs_Z.png',format='png')
show() 


#-----------------------------------------------------------------------------
fig0=plt.figure()  # B vs X coord. at y=0 z=0 plane
dzz= (eqdsk_z[nzeqd-1]-eqdsk_z[0])/(nzeqd-1)  # [m]  grid spacing
iz0=  int( (0.01*zscan-eqdsk_z[0])/dzz ) +1
print 0.01*zscan,   eqdsk_z[iz0-1], eqdsk_z[iz0],  eqdsk_z[iz0+1]
#iz0= int(nzeqd/2) #This is the index for z=0 (approximately)  
z00= eqdsk_z[iz0] # The exact value of Z (will be printed in plot)  [m]
BvsX=bmodprofxz[iz0,:]
BvsX_min= np.amin(BvsX[:])
BvsX_max= np.amax(BvsX[:])
BvsX_lim= max(BvsX_min,BvsX_max)
f=freqcy
plt.subplot(221) #------------------------- B(X)
plt.hold(True)
plt.grid(True)
plt.title('$x-scan$  $at$ $y=0,$ $z=$' +"%5.3f" %(z00*100) +'$cm$' ,y=1.03)
plt.ylabel('$B(x)$  $(T)$' )
#plt.xlabel('$ X $   $ (cm) $')
if model_b !=4:  # not equal to 4
    #Not FRC: Plot the magnitude of B (as saved in bmodprofxz)
    plt.plot(eqdsk_x[:]*100,BvsX[:],'b',linewidth=1.5)
if model_b==4: # FRC
    # Reverse sign of B at x=r0_frc and plot the position of separatrix and B=0
    # set the vertical limits for this plot:
    #plt.ylim( (-BvsX_lim, BvsX_lim) ) 
    # Find indexes corresponding to |X|<r0_frc (null point):
    ir0= np.where( (eqdsk_x<r0_frc/100.) & (eqdsk_x>-r0_frc/100.) )
    #print 'ir0=',ir0
    # Find indexes corresponding to |X|>r0_frc :
    irn= np.where( (eqdsk_x<-r0_frc/100.) ) # where X<-r0
    irp= np.where( (eqdsk_x> r0_frc/100.) ) # where X>+r0
    #print 'irs=',irs
    # Plot B inside null points as negative:
    plt.plot(eqdsk_x[ir0]*100,-BvsX[ir0],'b',linewidth=1.5)
    # and outside of null points as positive:
    plt.plot(eqdsk_x[irn]*100,+BvsX[irn],'b',linewidth=1.5)
    plt.plot(eqdsk_x[irp]*100,+BvsX[irp],'b',linewidth=1.5)  
    ## Plot B inside null points as positive:
    # plt.plot(eqdsk_x[ir0]*100,+BvsX[ir0],'b',linewidth=1.5)
    ## and outside of null points as negative:
    # plt.plot(eqdsk_x[irn]*100,-BvsX[irn],'b',linewidth=1.5)
    # plt.plot(eqdsk_x[irp]*100,-BvsX[irp],'b',linewidth=1.5)    
    # Designate null point and separatrix point by dashed lines:
    plt.plot([+rs_frc, +rs_frc],[-BvsX_lim, BvsX_lim],'r--')
    plt.plot([+r0_frc, +r0_frc],[-BvsX_lim, BvsX_lim],'r--')
    plt.plot([-rs_frc, -rs_frc],[-BvsX_lim, BvsX_lim],'r--')
    plt.plot([-r0_frc, -r0_frc],[-BvsX_lim, BvsX_lim],'r--')
    # Also plot the horizontal line for B=0 level:
    plt.plot(eqdsk_x[:]*100,BvsX[:]*0,'k',linewidth=1.)
    #--------------
plt.subplot(223) #------------------------- omega/omega_c vs X
plt.hold(True)
plt.grid(True)
if Nsp>1:
    isp=isp_wc-1 # ion (isp_wc is the genray-c index; isp is the python index)
    mime= mass[isp]/mass[0]      # m_s/m_e ratio
    Z_s= charge[isp]/charge[0]    # Z_s/e
    plt.ylabel('$\omega/\omega_{c}$ ['+isp_name[isp]+']')
else:  # electrons only: 
    isp=0 # electrons
    mime=1.0
    Z_s=1.0
    plt.ylabel('$\omega/\omega_{ce}$')
plt.xlabel('$x$   $ (cm) $')
plt.yscale('log') # log10 scale !!! Comment it, if you want a linear scale
plt.plot(eqdsk_x[:]*100,f/(28e9*BvsX[:]*Z_s/mime),'b',linewidth=1.5)       
#--------------
savefig('genray_B_vs_X.png',format='png')
show() 


#--------------------------------------------------------------------------

fig0=plt.figure()  # omega_ce/omega profiles vs x coord.
f=freqcy  #[0]
if Nsp>1:  # add ion species to calculate omega_LH
    isp=1
    msme= np.asscalar(mass[isp]/mass[0])      # ms/me ratio
    Zs= np.asscalar(charge[isp]/charge[0])    # Zs/e
    #print msme,Zs, Zs/msme
    wpi_w= wpw*(Zs/sqrt(msme)) #  == omega_pi/omega    BASED on ASSUMPTION ni=ne
    #print wpw/wpi_w
    wci_w= (Zs/msme)*wcw    #  == omega_ci/omega
    # Define omega_LH/omega. When wLH_w=1, eps_xx becomes 0 (resonance)
    wLH_w=sqrt( ((wcw*wci_w)**2 +(wci_w*wpw)**2 +(wpi_w*wcw)**2)/(wcw**2+wpw**2) )
    #From GENRAY-C coding [2022-04]
    #  Lower Hybrid resonance frequency :
    # (WLH/W)^2 = SUM[(Wpi/W)^2]*[Wce^2/(Wce^2 + Wpe^2)]
    # [ Assuming Wci^2 << W^2 << Wce^2 ]
    # wlwscan(i)= sqrt(xi*ye2/(ye2+xe)) ! (Lower hybrid)/omega

plt.subplot(221) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title( r"$f(GHz)=$"+r"$%1.4f$" %(f*1e-9) + \
   '  $Scan$ $at$ $(y,z)=($' +"%5.0f" %(yscan) +',' +\
   "%5.0f" %(zscan) +'$)cm$',y=1.15)
plt.ylabel('$\omega_{ce}/\omega$ $,$ $\omega_{pe}/\omega$ $,$ $\omega_{UH}/\omega$')
#plt.xlabel('$ x $   $ (cm) $')
plt.grid(True)
plt.plot(rscan[:],wcw[:],'b',linewidth=1.5) # fce/f
plt.plot(rscan[:],wpw[:],'r',linewidth=1.5) # fpe/f
plt.plot(rscan[:],wuw[:],'k',linewidth=1.5) # fuh/f
plt.plot(rscan[:],wlw[:],'g',linewidth=1.5) # flh/f [2022]
#if Nsp>1:  # ion species is present: plot omega/omega_LH
#    plt.plot(rscan[:],1./wLH_w[:],'g',linewidth=1.5) # w/wLH

wcw_max=np.amax(wcw[:])
wcw_min=np.amin(wcw[:])
wpw_max=np.amax(wpw[:])
wpw_min=np.amin(wpw[:])
wlw_max=np.amax(wlw[:])
wlw_min=np.amin(wlw[:])
wuw_max=np.amax(wuw[:])
wuw_min=np.amin(wuw[:])
print 'min/max of wpe/w =',wpw_min,wpw_max, '  (along Z[cm]=',zscan,')'
print 'min/max of wce/w =',wcw_min,wcw_max, '  (along Z[cm]=',zscan,')'
print 'min/max of wUH/w =',wuw_min,wuw_max, '  (along Z[cm]=',zscan,')'
print 'min/max of wLH/w =',wlw_min,wlw_max, '  (along Z[cm]=',zscan,')'
www_max=max(3*wcw_max,2*wpw_max)
plt.plot([0, Rmax],[0, 0],'k')
if model_b==4: # FRC
    plt.plot([rs_frc, rs_frc],[0, www_max],'r--')
    plt.plot([r0_frc, r0_frc],[0, www_max],'r--')
plt.xlim( (0., Rmax) )  
# set the upper limit for this plot:
plt.ylim( (0., www_max) ) 
#--------------
plt.subplot(222) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$\omega/\omega_{ce}$'+\
     ' $,$ $\omega/\omega_{pe}$'+\
     ' $,$ $\omega/\omega_{UH}$'+\
     ' $,$ $\omega/\omega_{LH}$' , y=1.03)
#plt.xlabel('$ x $   $ (cm) $')
plt.grid(True)
#plt.title('$Scan$  $at$ $(y,z)=($' +"%5.0f" %(yscan) +',' +\
#"%5.0f" %(zscan) +'$)cm$' ,y=1.03)
n_scan=rscan_size-2
plt.plot(rscan[:],1./wcw[:],'b',linewidth=1.5) # w/wce
plt.text(Rmax*0.75, (1./wcw[n_scan]),'$ce$', va='top',fontsize=fnt+4,color='b')
plt.plot(rscan[:],1./wpw[:],'r',linewidth=1.5) # w/wpe
plt.text(Rmax*0.75, (1./wpw[n_scan]),'$pe$', va='top',fontsize=fnt+4,color='r')
plt.plot(rscan[:],1./wuw[:],'k',linewidth=1.5) # w/wuh
plt.text(Rmax*0.85, (1./wuw[n_scan]),'$UH$', va='top',fontsize=fnt+4,color='k')
if np.min(wlw)>0:
    plt.plot(rscan[:],1./wlw[:],'g',linewidth=1.5) # w/wlh [2022]
    plt.text(Rmax*0.85, (1./wlw[n_scan]),'$LH$', va='top',fontsize=fnt+4,color='g')
plt.plot([0, Rmax],[0, 0],'k')
owcw_maxa=np.amax(1./wcw[n_scan]) 
owcw_mina=np.amin(1./wcw[n_scan])   
nwwc_max= ceil(owcw_maxa) 
nwwc_min= ceil(owcw_mina) 
nj= floor(1./wcw[0])
for i in range(1,rscan_size,1):
    if (1./wcw[i-1] -nj)==0:
        #print rscan[i-1], nj,  1./wcw[i-1], 1./wcw[i]
        plt.plot(rscan[i-1],1./wcw[i-1],'k.') # w/wce
        nj=nj-1        
    elif (1./wcw[i-1] -nj)*(1./wcw[i] -nj) <0:
        #print rscan[i-1], nj,  1./wcw[i-1], 1./wcw[i]
        plt.plot(rscan[i-1],1./wcw[i-1],'k.') # w/wce
        nj=nj-1

#if Nsp>1:  # ion species is present: plot omega/omega_LH
#    plt.plot(rscan[:],1./wLH_w[:],'g',linewidth=1.5) # w/wLH
#    plt.text(Rmax*0.8,1./wLH_w[n_scan],'$LH$',va='top',fontsize=fnt+4,color='g')

if model_b==4: # FRC
    plt.plot([rs_frc, rs_frc],[0, nwwc_max],'r--')
    plt.plot([r0_frc, r0_frc],[0, nwwc_max],'r--')
    plt.ylim( (0., nwwc_max*5) ) # set upper limit to 5*(omega/omega_c)_max
plt.xlim( (0., Rmax) )  
plt.yscale('log') # log10 scale !!! Comment it, if you want a linear scale
# To avoid large peaks where B->0,  set the upper limit for this plot:
#plt.ylim( (0., nwwc_max*5) ) # set upper limit to 5*(omega/omega_c)_max
#--------------
plt.subplot(223) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel(r'$\rho$ $coordinate$')
plt.xlabel('$ x $   $ (cm) $')
plt.grid(True)
plt.plot(rscan[:],rhoscan[:],'b',linewidth=1.5)  
plt.plot([0, Rmax],[0, 0],'k')
rho_max= rhoscan[n_scan]
if model_b==4: # FRC
    plt.plot([rs_frc, rs_frc],[0, rho_max],'r--')
    plt.plot([r0_frc, r0_frc],[0, rho_max],'r--')
plt.xlim( (0., Rmax) ) 
#--------------
plt.subplot(224) #-------------------------
plt.hold(True)
plt.grid(True)
if i_binvol==1:
    rho_mx= np.max(rho_bin_center[:])
    binvol_mx= np.max(binvol[:])
    text(0.5,binvol_mx,'$ binvol $   $ (cm^3) $')
    plt.xlabel(r'$\rho$ $coordinate$')
    plt.grid(True)
    plt.plot(rho_bin_center[:],binvol[:],'b',linewidth=1.5)  
    plt.plot([0, rho_mx],[0, 0],'k')
    plt.xlim( (0., rho_mx) )
    plt.ylim( (0., binvol_mx*1.1) ) 

savefig('genray_wcw_vs_R.png',format='png')
show() 
#stop


#--------------------------------------------------------------------------
fig0=plt.figure()  # Te and ne profiles vs x or y coord.
plt.subplot(221)
plt.hold(True)
plt.ylabel('$ T $   $ (keV) $')
plt.grid(True)
T_max= 1.05*np.amax(w_temp_vs_x_nc[:,:])+0.001  # upper limit for plots
x_max= np.max(w_x_densprof_nc[:]*100)
for i in range(0,Nsp,1):
    isp=i
    if remainder(i,4)==0: col='k'    #black
    if remainder(i,4)==1: col='b'    #blue
    if remainder(i,4)==2: col='r'
    if remainder(i,4)==3: col='g'    
    #plt.plot(rho_bin[:],temprof[i,:],color=col,linewidth=linw*(Nsp-i) )
    plt.plot(w_x_densprof_nc[:]*100,w_temp_vs_x_nc[i,:],color=col,linewidth=linw*(Nsp-i) )
    text((1-isp)*x_max*0.07, np.max(w_temp_vs_x_nc[i,:]), isp_name[isp],color=col)
if model_b==4: # FRC
    plt.plot([+rs_frc, +rs_frc],[0, T_max],'r--')
    plt.plot([+r0_frc, +r0_frc],[0, T_max],'r--')
    plt.plot([-rs_frc, -rs_frc],[0, T_max],'r--')
    plt.plot([-r0_frc, -r0_frc],[0, T_max],'r--')
    plt.title('$Dashed:$ $B=0$ $and$ $Separatrix$',y=1.03)
plt.axis([-Rmax,Rmax,0.,T_max])  

plt.subplot(223)
plt.hold(True)
plt.xlabel('$ x $   $ (cm) $')
#plt.xlabel(r'$\rho$')
plt.ylabel('$ n $   $ (cm^{-3}) $')
plt.grid(True)
n_max= 1.05*np.amax(w_dens_vs_x_nc[:,:])+1 # upper limit for plots
for i in range(0,Nsp,1):
    isp=i
    if remainder(i,4)==0: col='k'
    if remainder(i,4)==1: col='b'
    if remainder(i,4)==2: col='r'
    if remainder(i,4)==3: col='g'    
    #plt.plot(rho_bin[:],densprof[i,:],color=col,linewidth=linw*(Nsp-i) )
    #print i
    #print w_x_densprof_nc[:]*100
    #print w_dens_vs_x_nc[i,:]
    plt.plot(w_x_densprof_nc[:]*100,w_dens_vs_x_nc[i,:],color=col,linewidth=linw*(Nsp-i) )
    text((1-isp)*x_max*0.07,np.max(w_dens_vs_x_nc[i,:]), isp_name[isp],color=col)
if model_b==4: # FRC
    plt.plot([+rs_frc, +rs_frc],[0, n_max],'r--')
    plt.plot([+r0_frc, +r0_frc],[0, n_max],'r--')
    plt.plot([-rs_frc, -rs_frc],[0, n_max],'r--')
    plt.plot([-r0_frc, -r0_frc],[0, n_max],'r--')
plt.axis([-Rmax,Rmax,0.,n_max])  

plt.subplot(222)
plt.hold(True)
#plt.ylabel('$ T $   $ (keV) $')
plt.grid(True)
for i in range(0,Nsp,1):
    isp=i
    if remainder(i,4)==0: col='k'
    if remainder(i,4)==1: col='b'
    if remainder(i,4)==2: col='r'
    if remainder(i,4)==3: col='g'    
    plt.plot(w_y_densprof_nc[:]*100,w_temp_vs_y_nc[i,:],color=col,linewidth=linw*(Nsp-i) )
    text((1-isp)*x_max*0.07,np.max(w_temp_vs_y_nc[i,:]), isp_name[isp],color=col)
plt.axis([-Rmax,Rmax,0.,1.05*np.max(w_temp_vs_y_nc[:,:])+0.001])  

plt.subplot(224)
plt.hold(True)
plt.xlabel('$ y $   $ (cm) $')
#plt.ylabel('$ n $   $ (cm^{-3}) $')
plt.grid(True)
for i in range(0,Nsp,1):
    isp=i
    if remainder(i,4)==0: col='k'
    if remainder(i,4)==1: col='b'
    if remainder(i,4)==2: col='r'
    if remainder(i,4)==3: col='g'    
    plt.plot(w_y_densprof_nc[:]*100,w_dens_vs_y_nc[i,:],color=col,linewidth=linw*(Nsp-i))  
    text((1-isp)*x_max*0.07,np.max(w_dens_vs_y_nc[i,:]), isp_name[isp],color=col)    
plt.axis([-Rmax,Rmax,0.,1.05*np.amax(w_dens_vs_y_nc[0,:])+1])    
savefig('genray_profiles_T-n.png',format='png') # try pdf,eps,ps
show() 
#stop

#--------------------------------------------------------------------------

fig0=plt.figure()  # Characteristic frequencies vs R coord.
f9=freqcy/1e9  #[0] GHz units
plt.subplot(221) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title( '$-f_{pe}$    $--f_{UH}$    $-.-f_{ce}$' ,y=1.03)
plt.ylabel('$f/GHz$')
#plt.xlabel('$R$   $(cm)$')
plt.yscale('log') 
plt.grid(True)
plt.plot(rscan[:],wpw[:]*f9,'r',linewidth=1.5) # fpe
plt.plot(rscan[:],wuw[:]*f9,'k--',linewidth=1.5) # fuh
plt.plot(rscan[:],wcw[:]*f9,'b-.',linewidth=1.5) # fce  GHz
fce_max=np.amax(wcw[:]*f9)
fce_min=np.amin(wcw[:]*f9)
fpe_max=np.amax(wpw[:]*f9)
fpe_min=np.amin(wpw[:]*f9)
print 'max of fpe[GHz] =',fpe_max, '  (along Z[cm]=',zscan,')'
print 'max of fce[GHz] =',fce_max, '  (along Z[cm]=',zscan,')'
fff_max=max(3*fce_max,3*fpe_min)
plt.xlim( (0., Rmax) )  # limits in R axis in cm
# set the upper limit for this plot:
plt.ylim( (0.1, 100) )    # in GHz
#--------------
if Nsp>1:
    plt.subplot(223) #-------------------------
    plt.hold(True)
    plt.grid(True)
    plt.title( '$--f_{cH}/f$     $-.-f_{cD}/f$' ,y=1.03)
    plt.ylabel('$fc/f$')
    plt.xlabel('$R$   $(cm)$')
    #plt.yscale('log') 
    plt.grid(True)
    plt.plot(rscan[:],wcw[:]/1837,'g--',linewidth=1.5) # fcH/f  
    plt.plot(rscan[:],wcw[:]/3674,'b-.',linewidth=1.5) # fcD/f  
    plt.xlim( (0., Rmax) )  # limits in R axis in cm

if Nsp>1:
    plt.subplot(224) #-------------------------
    plt.hold(True)
    plt.grid(True)
    plt.title( '$--f/f_{cH}$     $-.-f/f_{cD}$' ,y=1.03)
    plt.ylabel('$f/fc$')
    plt.xlabel('$R$   $(cm)$')
    #plt.yscale('log') 
    plt.grid(True)
    plt.plot(rscan[:],1837/wcw[:],'g--',linewidth=1.5) # f/fcH
    plt.plot(rscan[:],3674/wcw[:],'b-.',linewidth=1.5) # f/fcD
    plt.xlim( (0., Rmax) )  # limits in R axis in cm

plt.subplot(222) #-------------------------
txt1= ' $Scan$  $at$ $(y,z)=($' +"%5.0f" %(yscan) +',' +\
"%5.0f" %(zscan) +'$)cm$' 
txt2=  r"$f (GHz) =$"+r"$%1.6f$" %(f*1e-9)
dx=0.7/3.0
plt.text(0.03, 0.7-dx*0, txt1 , va='center',fontsize=fnt+4)
plt.text(0.03, 0.7-dx*1, txt2 , va='center',fontsize=fnt+4) 
plt.axis([0., 1., 0., 1.])
plt.axis('off')
#--------------
savefig('genray_fce_vs_R.png',format='png')
show() 


if plot_Nperp_prof==1 & hot_roots==1: # ddd(r,Nper,Npar) profiles vs x,Nperp coord. (Hot roots)
#--------------------------------------------------------------------------
    X,Y = np.meshgrid(rscan, Npera)
    #D=ddd[0,:,:]
    #print shape(D)
    #print shape(X)
    for i in range(0,inpar,1):
        fig0=plt.figure()  # ddd(r,Nper,Npar) profiles vs x coord.
        plt.hold(True)
        plt.title('$Dispersion$ $function$ $D(x,N_{\perp})=0$ $levels$  $For$ $N_{||}=$'\
        +"%4.3f" %(Npara[i]),y=1.03)
        plt.ylabel('$ N_{\perp}$     $Red:$ $Hot$ $roots$ $(id=6)$')
        plt.xlabel('$x$ $(cm)$   $Scan$ $at$ $(y,z)=($' +"%5.0f" %(yscan) +',' +\
"%5.0f" %(zscan) +'$)cm$' )
        plt.grid(True)
        #plt.plot([0, Rmax],[0, 0],'k')
        D=ddd[i,:,:]
        #D_min=np.min(D)
        #D_max=np.max(D)
        #print 'D_min, D_max=', D_min, D_max # Usually -1e-24; +1e24
        #levels= np.arange(-100.,100.,2.0)
        #CS=plt.contour(X,Y,D[:],Ncont*5,linewidth=linw,cmap=plt.cm.jet)
        #CB=plt.colorbar(orientation='vertical', shrink=0.5, format='%.2e')
        if dispers_log_scale==1:
            plt.yscale('log')
        levels=[0] # level=0.0 of Real(D)
        CS=plt.contour(X,Y,D[:],levels,linewidth=linw,colors='r')
        Nper_im=Nperp_im[i,:]  # just to initialize
        for iscan in range(0,rscan_size,1):
            if Nperp_im[i,iscan]>1.0e-10:
                Nper_im[iscan]= Nperp_im[i,iscan]
            else:
                Nper_im[iscan]= 1.0e-5  # to avoid negative Nperp
        plt.plot(rscan[:],Nper_im[:],'k',linewidth=1.0)  # Imaginary/damping
        #axis([0,rscan_max,0.,np.max(Npera)]) 
        plt.xlim( (rscan_min, rscan_max) )
        if Nper_hot_ylim>0:
            plt.ylim( (1.e-3, Nper_hot_ylim) )
        else:
            plt.ylim( (1.e-3, 1.e4) )
        savefig('genray_hotD_id6_vs_R-Nperp_'+"%3.3f" %(Npara[i])+'.png',format='png')
        print 'file with hot roots is saved for Npara=',Npara[i]
        #show() 
    #stop

if plot_Nperp_prof==1:     # Nperp2 profiles vs x coord.   Cold roots
#--------------------------------------------------------------------------
    for i in range(0,inpar,1):
        fig0=plt.figure()  # Nperp2 profiles vs x coord.
        plt.hold(True)
        plt.title('$ioxm=-1 (blue)$  $and$  $ioxm=+1(green)$     $For$ $N_{||}=$'\
        +"%4.3f" %(Npara[i]))
        if dispers_log_scale==1:
            plt.ylabel('$N_{\perp} $     $Blue/Green:$ $Cold$ $plasma$ $roots$ $(id=2)$')
        else:
            plt.ylabel('$N_{\perp}^2 $     $Blue/Green:$ $Cold$ $plasma$ $roots$ $(id=2)$')
        plt.xlabel('$x$ $(cm)$   $Scan$ $at$ $(y,z)=($' +"%5.0f" %(yscan) +',' +\
"%5.0f" %(zscan) +'$)cm$' )
        plt.grid(True)
        if dispers_log_scale==1:
            plt.yscale('log') # !!!!!!!!!!!!!!!! log scale !!!!!!!!!!!!
            Nperp_m=np.zeros(rscan_size) #0*Nperp2m[i,:]  # just to initialize
            Nperp_p=np.zeros(rscan_size)  # just to initialize
            Nper_im=np.zeros(rscan_size)  # just to initialize
            #print 'min/max Nperp_p', np.min(Nperp_p),np.min(Nperp_p)
            for iscan in range(0,rscan_size,1):
                if Nperp2m[i,iscan]>1.0e-10:
                    Nperp_m[iscan]= np.sqrt(Nperp2m[i,iscan])
                else:
                    Nperp_m[iscan]= 1.0e-5  # to avoid negative Nperp^2
                if Nperp2p[i,iscan]>1.0e-10:
                    Nperp_p[iscan]= np.sqrt(Nperp2p[i,iscan])
                else:
                    Nperp_p[iscan]= 1.0e-5  # to avoid negative Nperp^2
                if Nperp_im[i,iscan]>1.0e-10:
                    Nper_im[iscan]= Nperp_im[i,iscan]
                else:
                    Nper_im[iscan]= 1.0e-5  # to avoid negative Nperp
            # Log10 scale for Nperp, with lower limit:
            plt.plot(rscan[:],Nper_im[:],'k',linewidth=1.0)  # Imaginary/damping
            plt.plot(rscan[:],Nperp_m[:],'b.',linewidth=2.5)  
            plt.plot(rscan[:],Nperp_p[:],'g.',linewidth=1.0)
            if Nper_cold_ylim>0:
                plt.ylim( (1.e-3, Nper_cold_ylim) )
            else:
                plt.ylim( (1.e-3, 1.e4) ) 
        else:
            # Linear scale for Nperp^2 (can be negative), no lower limit
            #plt.plot(rscan[:],Nper_im[:]^2,'k',linewidth=1.0)  # Imaginary/damping
            plt.plot(rscan[:],Nperp2m[i,:],'b',linewidth=1.0)  
            plt.plot(rscan[:],Nperp2p[i,:],'g',linewidth=0.5)
            plt.plot(rscan[:],rscan[:]*0,'k') # horizontal line Nperp=0
            #plt.ylim( (-2, 80**2) )
            plt.ylim((-1.e3, 5.e4))
            #plt.lim((-Nper_cold_ylim**2 , Nper_cold_ylim**2))
            # Vertical-axis Limits are set in above, for now.
            
        plt.xlim( (rscan_min, rscan_max) )
        savefig('genray_coldD_id2_vs_R_'+"%3.3f" %(Npara[i])+'.png',format='png')
        print 'file with cold roots is saved for Npara=',Npara[i]
        #show() 
    stop
    

if plot_Nperp_prof==1:     # Nper_im vs x coord.  
#--------------------------------------------------------------------------
    for i in range(0,inpar,1):
        fig0=plt.figure()  # Nperp2 profiles vs x coord.
        plt.hold(True)
        plt.title('        $Im(N_{\perp})$ $(damping)$   $For$ $N_{||}=$'\
        +"%4.3f" %(Npara[i]))
        #plt.ylabel('$N_{\perp} $     $Cold$ $plasma$ $roots$ $(id=2)$')
        plt.xlabel('$x$ $(cm)$   $Scan$ $at$ $(y,z)=($' +"%5.0f" %(yscan) +',' +\
"%5.0f" %(zscan) +'$)cm$' )
        plt.grid(True)
        #plt.yscale('log') # !!!!!!!!!!!!!!!! log scale !!!!!!!!!!!!
        #Nperp_m=Nperp2m[i,:]  # just to initialize
        #Nperp_p=Nperp2m[i,:]  # just to initialize
        Nper_im=Nperp_im[i,:]  # just to initialize
        Nper_im_min=np.amin(Nperp_im[i,:])
        Nper_im_max=np.amax(Nperp_im[i,:])
        print 'Nper_im  min/max:',Nper_im_min,Nper_im_max
        if Nper_im_max>Nper_im_min:
            #plt.plot(rscan[:],Nperp_m[:],'b',linewidth=2.5)  
            #plt.plot(rscan[:],Nperp_p[:],'g',linewidth=1.0)
            plt.plot(rscan[:],Nper_im[:],'k',linewidth=1.0)  # Imaginary/damping
            #plt.plot([0, rscan_max],[0, 0],'k')
            plt.xlim( (rscan_min, rscan_max) )
            if Nper_im_ylim>0:
                plt.ylim( (0, Nper_im_ylim) )
            #Nperp2m_max=np.max(Nperp2m[i,:])
            # Comment next line if you want automatic limits in Nperp (cold)
            #plt.ylim((-100.,Nperp2m_max))  # Un-comment if you want specific limits in Nperp
            savefig('genray_Nper_im_vs_R_'+"%3.3f" %(Npara[i])+'.png',format='png')
            #print 'file with cold roots is saved for Npara=',Npara[i]
            #show() 
    #stop    

#[Nov-2014] Added: power deposition profiles over (R,Z) rectangular grid.
#  [erg/sec]
Rgrid=dat.variables['Rgrid']
print 'Rgrid is: ', Rgrid.long_name, Rgrid.shape
Zgrid=dat.variables['Zgrid']
print 'Rgrid is: ', Zgrid.long_name, Zgrid.shape
NRgrid=Rgrid[:].size  
NZgrid=Zgrid[:].size  
spwr_rz_e=dat.variables['spwr_rz_e']
print 'spwr_rz_e is: ', spwr_rz_e.long_name #, spwr_rz_e.shape
spwr_rz_i=dat.variables['spwr_rz_i']
print 'spwr_rz_i is: ', spwr_rz_i.long_name #, spwr_rz_i.shape
spwr_rz_cl=dat.variables['spwr_rz_cl']
print 'spwr_rz_cl is: ', spwr_rz_cl.long_name #, spwr_rz_cl.shape

Rgrid= Rgrid[:]*100  # m->cm
Zgrid= Zgrid[:]*100  # m->cm
spwr_rz_e=  spwr_rz_e[:]/1.0e10  # kW
spwr_rz_i=  spwr_rz_i[:]/1.0e10  # kW
spwr_rz_cl= spwr_rz_cl[:]/1.0e10  # kW
sum_spwr_rz_e= sum(spwr_rz_e)
sum_spwr_rz_i= sum(spwr_rz_i)
sum_spwr_rz_cl=sum(spwr_rz_cl)
print 'sum(spwr_rz_e)=    ',sum_spwr_rz_e,   ' kW'
print 'sum(spwr_rz_i)=    ',sum_spwr_rz_i,   ' kW'
print 'sum(spwr_rz_cl)=   ',sum_spwr_rz_cl,  ' kW'

# Define min/max for major radius range:
if n_wall==0:  # wall is not defined; get limits from eqdsk:
    Rmax=np.amax(eqdsk_r)*100
    Rmin=np.amin(eqdsk_r)*100 # meters -> cm
else:
    Rmax=np.amax(r_wall)
    Rmin=np.amin(r_wall)
    if Zmin==0: Zmin=np.amin(z_wall)
    if Zmax==0: Zmax=np.amax(z_wall)
    Rmax=max(Rmax,abs(Rmin))
    Rmin=-Rmax
    if Zmax==0: Zmax=min(Zmax,Zgrid[NZgrid-1])
    if Zmin==0: Zmin=max(Zmin,Zgrid[0])
    
print 'Rmin,Rmax=', Rmin,Rmax , 'Zmin,Zmax=', Zmin,Zmax    


    
# Define boundary for top-view (R-phi) plots:
phi=arange(301.)*2*pi/300.  # tor.angle in rad.
bndry_in_X= Rmin*cos(phi)
bndry_in_Y= Rmin*sin(phi)
bndry_out_X= Rmax*cos(phi)
bndry_out_Y= Rmax*sin(phi)







#--------------------------------------------------------------------------
if Nsp>1:   # For e or ions, if present
    fig0=plt.figure()  # Valfven, nperp_Alfven, (Kperp*rho) profiles vs x or y coord., at z=0
    plt.subplot(221)
    plt.hold(True)
    plt.ylabel('$ V_{Alfven} $   $ (cm/s) $')
    plt.title('$Profiles$ $at$ $Z=0$',y=1.02)
    plt.grid(True)
    
    # find sum(density[i]*mass[i]/mass_proton)
    i=0 # electrons
    nm= w_dens_vs_x_nc[i,:]*mass[i]/mp
    for i in range(1,Nsp,1):
        nm= nm+ w_dens_vs_x_nc[i,:]*mass[i]/mp  #summed over all species; func of x
    # Alfven speed:
    VA=2.18e15*w_bmod_vs_x_nc[:]/sqrt(nm[:])
    print 'nm=',np.max(nm), '  max(w_bmod_vs_x_nc)=', np.max(w_bmod_vs_x_nc)
    for i in range(1,Nsp,1):
        if remainder(i,4)==0: col='k'
        if remainder(i,4)==1: col='b'
        if remainder(i,4)==2: col='r'
        if remainder(i,4)==3: col='g'    
        plt.plot(w_x_densprof_nc[:]*100,VA,color=col,linewidth=linw*(Nsp-i) )
    #plt.axis([-Rmax,Rmax,0.,T_max])  #
    
    plt.subplot(223)
    plt.hold(True)
    plt.xlabel('$ x $   $ (cm) $')
    plt.grid(True)
    for i in range(1,Nsp,1):
        if remainder(i,4)==0: col='k'
        if remainder(i,4)==1: col='b'
        if remainder(i,4)==2: col='r'
        if remainder(i,4)==3: col='g'    
        isp=isp_wc-1 # ion (isp_wc is the genray-c index; isp is the python index) 
        mime= mass[isp]/mass[0]      # m_s/m_e ratio
        Z_s= charge[isp]/charge[0]    # Z_s/e
        wc_w_x= (28e9*w_bmod_vs_x_nc[:]*Z_s/mime)/f  # omega_ci/omega
        w_wc_x= 1/wc_w_x
        plt.ylabel('$\omega_{c}/\omega$ ['+isp_name[isp]+']')
        plt.plot(w_x_densprof_nc[:]*100,wc_w_x,color=col,linewidth=linw*(Nsp-i) )
    #plt.axis([-Rmax,Rmax,0.,n_max])  
    
    plt.subplot(222)
    plt.hold(True)
    plt.grid(True)
    plt.title('$N_{\perp}$ $=$ $c/V_{Alfven}$',y=1.02)
    for i in range(1,Nsp,1):
        if remainder(i,4)==0: col='k'
        if remainder(i,4)==1: col='b'
        if remainder(i,4)==2: col='r'
        if remainder(i,4)==3: col='g'    
        #isp=isp_wc-1 # ion (isp_wc is the genray-c index; isp is the python index)
        #mime= mass[isp]/mass[0]      # m_s/m_e ratio
        #Z_s= charge[isp]/charge[0]    # Z_s/e
        #wci_w= (28e9*w_bmod_vs_x_nc[:]*Z_s/mime)/f  # omega/omega_ci
        nper_A= clight/VA
        plt.plot(w_y_densprof_nc[:]*100,nper_A,color=col,linewidth=linw*(Nsp-i) )
    #plt.axis([-Rmax,Rmax,0.,1.05*np.amax(w_temp_vs_y_nc[0,:])+0.001])
    #plt.yscale('log') # log10 scale !!! Comment it, if you want a linear scale    
    
    plt.subplot(224)
    plt.hold(True)
    plt.xlabel('$ x $   $ (cm) $')
    plt.grid(True)
    for i in range(1,Nsp,1):
        if remainder(i,4)==0: col='k'
        if remainder(i,4)==1: col='b'
        if remainder(i,4)==2: col='r'
        if remainder(i,4)==3: col='g'  
        Vthi= 1.87e9*sqrt(w_temp_vs_x_nc[i,:]/mass[i]) # cm/s   sqrt(2T/m)
        isp=i # ion species 
        mime= mass[isp]/mass[0]      # m_s/m_e ratio
        Z_s= charge[isp]/charge[0]    # Z_s/e
        wci_w= (28e9*w_bmod_vs_x_nc[:]*Z_s/mime)/f  # omega_ci/omega
        Kper_rho= (Vthi/VA)/wci_w
        # Kperp based on Alfven speed:  omega/VA
        # Larmor radius for thermal ion species: Vthi/omega_ci
        # Then, Kperp*rhoi= (Vthi/VA)*(omega/omega_ci)
        plt.title(r'$estimate$ $for$ $K_{\perp,A}{\rho_{Ti}}$  ['+isp_name[isp]+']',y=1.02)
        plt.plot(w_x_densprof_nc[:]*100,Kper_rho,color=col,linewidth=linw*(Nsp-i))   
    #plt.axis([-Rmax,Rmax,0.,1.05*np.amax(w_dens_vs_y_nc[0,:])+1])    
    #plt.yscale('log') # log10 scale !!! Comment it, if you want a linear scale    
    #plt.ylim((0,100.)) #upper limit: Kperp*rho=100 (can be too large when B~0)
    savefig('genray_profiles_VA.png',format='png') # try pdf,eps,ps
    show() 
#--------------------------------------------------------------------------
#stop




#============== MORE DATA, FROM RAY-TRACING ====================================
print '---------------------------------------- RAYS DATA and PLOTS'
# Averaged current densities
s_cur_den_parallel=dat.variables['s_cur_den_parallel']
print 's_cur_den_parallel is: ', s_cur_den_parallel.long_name, s_cur_den_parallel.shape
s_cur_den_onetwo=dat.variables['s_cur_den_onetwo']
print 's_cur_den_onetwo is: ', s_cur_den_onetwo.long_name, s_cur_den_onetwo.shape
s_cur_den_toroidal=dat.variables['s_cur_den_toroidal']
print 's_cur_den_toroidal is: ', s_cur_den_toroidal.long_name, s_cur_den_toroidal.shape
s_cur_den_poloidal=dat.variables['s_cur_den_poloidal']
print 's_cur_den_poloidal is: ', s_cur_den_poloidal.long_name, s_cur_den_poloidal.shape

# Power profiles [erg/sec/cm^3]
powden=dat.variables['powden']
print 'powden is: ', powden.long_name, powden.shape
powden_e=dat.variables['powden_e']
print 'powden_e is: ', powden_e.long_name, powden_e.shape
powden_i=dat.variables['powden_i']
print 'powden_i is: ', powden_i.long_name, powden_i.shape

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
        print 'powden_cl is: ', powden_cl.long_name, powden_cl.shape
finally:
    print '----------------------------------------'
    
    

# Total power [erg/sec]
print '---------------------------------------- [erg/sec]'
power_total=dat.variables['power_total'].getValue()  #getValue() for scalar
print 'power_total is: ', power_total
power_inj_total=dat.variables['power_inj_total'].getValue()  #getValue() for scalar
print 'power_inj_total is: ', power_inj_total
powtot_e=dat.variables['powtot_e'].getValue()  #getValue() for scalar
print 'powtot_e is: ', powtot_e
powtot_i=dat.variables['powtot_i'].getValue()  #getValue() for scalar
print 'powtot_i is: ', powtot_i
powtot_cl=dat.variables['powtot_cl'].getValue()  #getValue() for scalar
print 'powtot_cl is: ', powtot_cl
# For iabsorp.eq.3 only:
#powtot_s=dat.variables['powtot_s'].getValue()  #getValue() for scalar
#print 'powtot_s is: ', powtot_s
print '---------------------------------------- [erg/sec]'

delpwr=dat.variables['delpwr']
print 'delpwr is: ', delpwr.long_name, delpwr.shape

# CONVERT to kW ===================================================
powden_e= powden_e[:]/1.0e7  # W/cm^3 (unless binvol=1.0, then W)
powden_i= powden_i[:]/1.0e7  # W/cm^3 (unless binvol=1.0, then W)
powden_cl=powden_cl[:]/1.0e7 # W/cm^3 (unless binvol=1.0, then W)
powden=   powden[:]/1.0e7    # W/cm^3 (unless binvol=1.0, then W)
#------------
powtot_e=    powtot_e/1.0e10  # kW
powtot_i=    powtot_i/1.0e10  # kW
powtot_cl=   powtot_cl/1.0e10 # kW
power_total= power_total/1.0e10 # kW
power_inj_total= power_inj_total/1.0e10 # kW
delpwr= delpwr[:,:]/1.0e10 # kw


ws=dat.variables['ws']
print 'ws is: ', ws.long_name, ws.shape


#sdpwr=dat.variables['sdpwr']
#print 'sdpwr is: ', sdpwr.long_name, sdpwr.shape
spower=dat.variables['spower']
print 'spower is: ', spower.long_name, spower.shape


# RAY trajectories
wx=dat.variables['wx']
print 'wx is: ', wx.long_name, wx.shape
wy=dat.variables['wy']
print 'wy is: ', wy.long_name, wy.shape
wr=dat.variables['wr']
print 'wr is: ', wr.long_name, wr.shape
wz=dat.variables['wz']
print 'wz is: ', wz.long_name, wz.shape
wphi=dat.variables['wphi']
print 'wphi is: ', wphi.long_name, wphi.shape
# Number of elements for each ray:
nrayelt=dat.variables['nrayelt']
print 'nrayelt is: ', nrayelt.long_name, nrayelt.shape
# Refractive indices along rays:
wnper=dat.variables['wnper']
print 'wnper is: ', wnper.long_name, wnper.shape
wnpar=dat.variables['wnpar']
print 'wnpar is: ', wnpar.long_name, wnpar.shape
wn_x=dat.variables['wn_x']
print 'wn_x is: ', wn_x.long_name, wn_x.shape
wn_y=dat.variables['wn_y']
print 'wn_y is: ', wn_y.long_name, wn_y.shape
wn_r=dat.variables['wn_r']
print 'wn_r is: ', wn_r.long_name, wn_r.shape
wn_z=dat.variables['wn_z']
print 'wn_z is: ', wn_z.long_name, wn_z.shape
wn_phi=dat.variables['wn_phi']
print 'wn_phi is: ', wn_phi.long_name, wn_phi.shape
# E-wave-field along rays
cwexde=dat.variables['cwexde']
print 'cwexde is: ', cwexde.long_name, cwexde.shape
cweyde=dat.variables['cweyde']
print 'cweyde is: ', cweyde.long_name, cweyde.shape
cwezde=dat.variables['cwezde']
print 'cwezde is: ', cwezde.long_name, cwezde.shape
# fluxn = Power Flux along rays:
fluxn=dat.variables['fluxn']
print 'fluxn is: ', fluxn.long_name, fluxn.shape
# Total magnetic field along rays:
sbtot=dat.variables['sbtot']
print 'sbtot is: ',sbtot.long_name, sbtot.shape, sbtot.units
# Bz magnetic field along rays [G]:
sb_z=dat.variables['sb_z']
print 'sb_z is: ',sb_z.long_name, sb_z.shape, sb_z.units
# Br magnetic field along rays [G]:
sb_r=dat.variables['sb_r']
print 'sb_r is: ',sb_r.long_name, sb_r.shape, sb_r.units
# --- Also availabe in genray.nc file: sb_x, sb_y, and sb_phi (toroidal)
# density along rays [1/cm^3]
sene=dat.variables['sene']
print 'sene is: ', sene.long_name, sene.shape, sene.units

# Number of rays
Nrays=wz[:,0].size  
print 'Number of rays: Nrays=',Nrays
Nrayelt=wz[0,:].size
print 'Max Number of ray elements: Nrayelt=',Nrayelt

# Kim [1/cm] along rays, for electrons+ions (dep. on i_salphal setting):
salphal=dat.variables['salphal']
print 'salphal is: ', salphal.long_name, salphal.shape

# Kim for ion species:
if (Nsp>1) & (isp_wc>1):
    print 'If salphas is unavailable, change isp_wc to 1 and rerun this script'
    salphas=dat.variables['salphas'] # [ksp=0:nbulk,nrays,nrelta]
    print 'salphas is: ', salphas.long_name, salphas.shape, ' (Nsp-1,Nrays-1,Nrayelt-1)'
    salphas=np.asarray(salphas)
    print 'salphas is: ', salphas.shape #, salphas[Nsp-1,Nrays-1,Nrayelt-1]
    salphas_all_i=np.zeros((Nrays,Nrayelt))
    print 'salphas_all_i is: ', salphas_all_i.shape
    salphas_all_i= np.sum(salphas, axis=0)
    print 'salphas_all_i is: ', salphas_all_i.shape


# Vgroup normalized to c:
vgr_x=dat.variables['vgr_x']
print 'vgr_x is: ', vgr_x.long_name, vgr_x.shape
vgr_y=dat.variables['vgr_y']
print 'vgr_y is: ', vgr_y.long_name, vgr_y.shape
vgr_z=dat.variables['vgr_z']
print 'vgr_z is: ', vgr_z.long_name, vgr_z.shape

# Added Jan-2016: VTi along rays, for each species:
vthermal=dat.variables['vthermal']
print 'vthermal [cm/s] is: ', vthermal.long_name, vthermal.shape 
# Based on wvthermal(1,k)=1.87d9*dsqrt(tempe_xyz(x,y,z,k)/dmas(k)), see genray.f
# vth=sqrt(2T/m),    vth in (cm/sec),T(keV)
# vth= 1.87e9*sqrt(T/dmas)


print 'Bz at ray starting point:',sb_z[0,0], ' Gauss'

print '----------------------------------------'



#--------------------------------------------------------------------------
fig2=plt.figure()
ax = plt.subplot(1, 1, 1)
plt.hold(True)
plt.title('$Power$ $in$ $ray$ $channel$  $(kW)$' , fontsize=26,y=1.03)
plt.xlabel('$distance$ $along$ $ray$  $(cm)$' , fontsize=26)
plt.ylabel('$ delpwr $  $ (kW) $' ,fontsize=26)
ws_max=0.
plt.grid(True)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    ws_max1=np.amax(ws[:,:])
    ws_max=np.amax([ws_max1,ws_max])
    linwi=6*linw-0.6*i
    linwi=max(1,linwi)
    plt.plot(ws[i,0:Nm],delpwr[i,0:Nm],color=col,linewidth=linwi)
plt.axis([0.,1.05*ws_max,0.,1.05*np.amax(delpwr[:,0])] )
plt.savefig('genray_rays_delpwr.png') 
show()



#--------------------------------------------------------------------------
if i_cps==1:
    fig3=plt.figure()  # starting values over (R,Npar) plane
    plt.hold(True)
    plt.grid(True)
    plt.title('$Init.$  $values$  $for$  $each$ $ray.$   $Large$ $circle$ $-$ $got$ $to$ $detector$',y=1.03)
    plt.xlabel('$R$  $(cm)$')
    plt.ylabel('$Refr.$ $index$  $N_{||}$ ')

    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
        print i,Nm
        Nparstart=wnpar[i,0] #  refractive index at starting point
        Rstart=   wr[i,0]    #  R at starting point [cm]
        if i_bad_initial_conditions[i]==0: # only show those that were started
            if i_hit_detector_cps[i]==1:
                plt.plot(Rstart,Nparstart,'o',color=col)
            else:
                plt.plot(Rstart,Nparstart,'.',color=col)
    
    plt.savefig('genray_start_in_R-Npar.png') 
    show()
#--------------------------------------------------------------------------



#--------------------------------------------------------------------------
fig4=plt.figure()  # RAYS in top-view (R-phi or X-Y)
#ax = plt.subplot(1, 1, 1)
#ax.set_aspect(1.0)
plt.hold(True)
plt.grid(True)
plt.title('$Rays$  $x(t),y(t)$  $and$  $Electron$ $Density$ $at$ $z=0$  $(cm^{-3})$',y=1.03)
plt.xlabel('$X$  $(cm)$')
plt.ylabel('$Y$  $(cm)$')
X,Y = np.meshgrid(eqdsk_x, eqdsk_y)
X=X*100
Y=Y*100
D=densprofxy
#CS=plt.contour(X,Y,D[:],Ncont,linewidth=linw,cmap=cm.jet)
#CB=plt.colorbar(orientation='vertical',shrink=0.5,format='%.2e') 
#CB.add_lines
#if ios=='windows': 
#    CB.lines.set_linewidth(5)
#CB.set_label("$(cm^{-3})$", size=14)

for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    print i,Nm
    X=wx[i,0:Nm]
    Y=wy[i,0:Nm]
    plt.plot(X,Y,color=col,linewidth=linw)  
    # plot arrow for refractive vector (Nx,Nz) at starting point: 
    Nx0=wn_x[i,0]
    Ny0=wn_y[i,0]
    N0=sqrt(Nx0*Nx0+Ny0*Ny0)
    #plt.arrow(wx[i,0],wy[i,0],(Nx0/N0)*arr_len,(Ny0/N0)*arr_len,'->')
    # plot small dot at the launching point:
    plt.plot(wx[i,0],wy[i,0],'k.')
    
plt.plot(wx[0,0],wy[0,0],'ko') # small circle at launching point
plt.plot(bndry_in_X,  bndry_in_Y,  'k',linewidth=linw*2)
plt.plot(bndry_out_X, bndry_out_Y, 'k',linewidth=linw*2)
if i_cps==1:
    #------------- To designate CPS detector, draw a circle around its position:
    cps_phi=arange(31.)*2*pi/30.  # simply an arc angle [0;2pi] in rad.
    detector_cps_circle_x= detector_cps_x +0.5*detector_cps_diam*cos(cps_phi)
    detector_cps_circle_y= detector_cps_y +0.5*detector_cps_diam*sin(cps_phi)
    # [m] to [cm]
    detector_cps_circle_x= detector_cps_circle_x*100
    detector_cps_circle_y= detector_cps_circle_y*100
    plt.plot(detector_cps_circle_x, detector_cps_circle_y, 'k',linewidth=linw)
    #-------------
plt.axis('equal')    
plt.savefig('genray_rays_in_R-phi.png') 
show()
#stop


#--------------------------------------------------------------------------
fig0=plt.figure() # RAYS and Density profile in cross-sectional view X-Z
ax = plt.subplot(111)
ax.set_aspect(1.0)
plt.axis([Zmin,Zmax,Rmin*1.1,Rmax*1.1])
plt.hold(True)
xmin=amin(eqdsk_x)*100*1.1  # Limits for plots
xmax=amax(eqdsk_x)*100*1.1  # m->cm, and add abit
ymin=amin(eqdsk_y)*100*1.1  # Limits for plots
ymax=amax(eqdsk_y)*100*1.1  # m->cm, and add abit
zmin=amin(eqdsk_z)*100*1.1  # Limits for plots
zmax=amax(eqdsk_z)*100*1.1  # m->cm, and add abit
#plt.title(r'$Rays,$ $n_e$ $(cm^{-3})$ $and$ $levels$ $\omega / \omega_c$  ['+isp_name[isp_wc-1]+']',y=1.03)
plt.title(r'$Rays$ $and$ $n_e$ $(cm^{-3})$',y=1.03)
plt.ylabel('$X$  $(cm)$')
plt.xlabel('$Z$  $(cm)$')
X,Zx = np.meshgrid(eqdsk_x, eqdsk_z)
X=X*100
Zx=Zx*100
D=densprofxz # during Genray-C run: saved for 1st species (k=1)
Dmin=np.min(D)
Dmax=np.max(D)
print 'np.min(D), np.max(D), Dmax-Dmin=', Dmin, Dmax, Dmax-Dmin
if (Dmax-Dmin)/(Dmin+Dmax)>1e-8:
    CS=plt.contour(Zx,X,D[:],Ncont,linewidth=linw,cmap=cm.jet)
    CB=plt.colorbar(orientation='horizontal',shrink=1.0,format='%1.0e')
    if ios=='windows': 
        CB.lines.set_linewidth(10) # for color lines in colorbar
    CB.set_label("$(cm^{-3})$", size=14)


#plt.plot(eqdsk_z, densprofxz[32,:])
#plt.grid(True)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    # To remove the last point: set Nm= np.asscalar(nrayelt[i])-2
    print i,Nm
    plt.plot(wz[i,0:Nm],wx[i,0:Nm],color=col,linewidth=linw)
    # plot arrow for refractive vector (Nx,Nz) at starting point: 
    #Nx0=wn_x[i,0]
    #Nz0=wn_z[i,0]
    #N0=sqrt(Nx0*Nx0+Nz0*Nz0)
    #plt.arrow(wx[i,0],wz[i,0],(Nx0/N0)*arr_len,(Nz0/N0)*arr_len,'->') 
    # plot small dot at the launching point:
    plt.plot(wz[i,0],wx[i,0],'k.')
    
CS_wwc=plt.contour(Zx,X,w_wc[:],levels_wwc,linewidth=linw,cmap=cm.jet)
if contour_labels==1:
    plt.clabel(CS_wwc,inline=1,orientation='horizontal',fmt='%1.0f',fontsize=11,colors='k',manual=True)
if contour_labels==2:
    plt.clabel(CS_wwc,inline=1,orientation='horizontal',fmt='%1.0f',fontsize=11,colors='k',manual=False)

# plot small circle at the launching point:
plt.plot(wz[0,0],wx[0,0],'ko')
# plot walls, if any:
if n_wall>0:
    plt.plot(z_wall,r_wall,'k',linewidth=linw*2)
    
if i_cps==1:
    #------------- To designate CPS detector, draw a circle around its position:
    cps_phi=arange(31.)*2*pi/30.  # simply a [0;2pi] arc angle in rad.
    detector_cps_circle_x= detector_cps_x +0.5*detector_cps_diam*cos(cps_phi)
    detector_cps_circle_z= detector_cps_z +0.5*detector_cps_diam*sin(cps_phi)
    # [m] to [cm]
    detector_cps_circle_x= detector_cps_circle_x*100
    detector_cps_circle_z= detector_cps_circle_z*100
    plt.plot(detector_cps_circle_z, detector_cps_circle_x, 'k',linewidth=linw)
    # Here we have [Z,X] plotting coords !
    #-------------    
plt.savefig('genray_rays_Ne_inXZ.png') 
plt.show()


#--------------------------------------------------------------------------
fig0=plt.figure() # RAYS and Density profile in cross-sectional view Y-Z
ax = plt.subplot(121)
ax.set_aspect(1.0)
plt.axis([Rmin*1.1,Rmax*1.1,Zmin,Zmax])
plt.hold(True)
plt.title('$Rays$ $and$ $Density$ $(cm^{-3})$',y=1.03)
plt.xlabel('$Y$  $(cm)$')
plt.ylabel('$Z$  $(cm)$')
Y,Zx = np.meshgrid(eqdsk_y, eqdsk_z)
Y=Y*100
Zx=Zx*100
D=densprofyz # during Genray-C run: saved for last species (k=nbulk)
Dmin=np.min(D)
Dmax=np.max(D)
print 'np.min(D), np.max(D), Dmax-Dmin=', Dmin, Dmax, Dmax-Dmin
if (Dmax-Dmin)/(Dmin+Dmax)>1e-8:
    CS=plt.contour(Y,Zx,D[:],Ncont,linewidth=linw,cmap=cm.jet)
    CB=plt.colorbar(orientation='vertical', shrink=0.5, format='%.2e') 
    if ios=='windows': 
        CB.lines.set_linewidth(5)
    CB.set_label("$(cm^{-3})$", size=14)
#plt.plot(eqdsk_y, densprofyz[32,:])
plt.grid(True)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(wy[i,0:Nm],wz[i,0:Nm],color=col,linewidth=linw)
    # plot arrow for refractive vector (Ny,Nz) at starting point: 
    Ny0=wn_y[i,0]
    Nz0=wn_z[i,0]
    N0=sqrt(Ny0*Ny0+Nz0*Nz0)
    #plt.arrow(wy[i,0],wz[i,0],(Ny0/N0)*arr_len,(Nz0/N0)*arr_len,'->')
    # plot small dot at the launching point:
    plt.plot(wy[i,0],wz[i,0],'k.')
    
# plot small circle at the launching point:
plt.plot(wy[0,0],wz[0,0],'ko')
# plot walls, if any:
if n_wall>0:
    plt.plot(r_wall,z_wall,'k',linewidth=linw*2)
#-----------------------------------------------------    

isp=Nsp-1 # last ion species, or electrons if Nsp=1 (one species only)
mime= mass[isp]/mass[0]        # m_s/m_e ratio
Z_s= charge[isp]/charge[0]    # Z_s/e
ax = plt.subplot(122)
ax.set_aspect(1.0)
plt.axis([Rmin*1.1,Rmax*1.1,Zmin,Zmax])
plt.hold(True)
plt.title('$n$ $(cm^{-3})$ $for$ '+isp_name[isp],y=1.03)
plt.xlabel('$Y$  $(cm)$')
#plt.ylabel('$Z$  $(cm)$')
if (Dmax-Dmin)/(Dmin+Dmax)>1e-8:
    CS=plt.contour(Y,Zx,D[:],Ncont,linewidth=linw,cmap=cm.jet)
    CB=plt.colorbar(orientation='vertical', shrink=0.5, format='%.2e') 
    if ios=='windows': 
        CB.lines.set_linewidth(5)
    CB.set_label("$(cm^{-3})$", size=14)
# plot walls, if any:
if n_wall>0:
    plt.plot(r_wall,z_wall,'k',linewidth=linw*2)
    
plt.savefig('genray_rays_n_inYZ.png') 
plt.show()



#--------------------------------------------------------------------------
fig3=plt.figure()  # RAYS and Pol.Flux in cross-section view X-Z
ax = plt.subplot(1, 1, 1)
ax.set_aspect(1.0)
plt.axis([Zmin,Zmax,Rmin*1.1,Rmax*1.1])
plt.hold(True)
plt.ylabel('$X$  $(cm)$')
plt.xlabel('$Z$  $(cm)$')
#plt.title(r'$Rays,$ $\Psi_p/2\pi$ $and$ $\omega / \omega_c$  ['+isp_name[isp_wc-1]+']',y=1.03)
plt.title(r'$Rays$ $and$ $\Psi_p/2\pi$',y=1.03)
#R,Zr = np.meshgrid(eqdsk_x[ix0:ix0*2], eqdsk_z[:])
R,Zr = np.meshgrid(eqdsk_r[:], eqdsk_z[:])
R=R*100
Zr=Zr*100
PSI = eqdsk_psi
# For setting the min/max levels in contour plots of PSI:
PSImin=psimag 
PSImax=psilim +0.30*(psilim-psimag) # consider to plot few lines outside of LCFS
                                    # that's why 30% is added
print 'psimag,psilim=', psimag,psilim, '  PSImin,PSImax=', PSImin,PSImax


levels=np.arange(PSImin,PSImax,(PSImax-PSImin)/(Ncont-1)) #all Ncont levels here
print 'min(PSI),max(PSI)=',np.min(PSI),np.max(PSI)
print 'levels of PSI:',levels

CS=plt.contour(Zr, R,PSI[:],levels,linewidth=linw,cmap=cm.jet)
CS=plt.contour(Zr,-R,PSI[:],levels,linewidth=linw,cmap=cm.jet)
CB=plt.colorbar(orientation='horizontal', shrink=1.0, format='%.1e') 
if ios=='windows': 
    CB.lines.set_linewidth(10)
CB.set_label("$(Tesla*m^2)$", size=14)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(wz[i,0:Nm],wx[i,0:Nm],color=col,linewidth=linw)
    # plot arrow for refractive vector (Nx,Nz) at starting point: 
    Nx0=wn_x[i,0]
    Nz0=wn_z[i,0]
    N0=sqrt(Nx0*Nx0+Nz0*Nz0)
    #plt.arrow(wx[i,0],wz[i,0],(Nx0/N0)*arr_len,(Nz0/N0)*arr_len,'->',color=col) 
    # plot small dot at the launching point:
    plt.plot(wz[i,0],wx[i,0],'k.')

CS_wwc=plt.contour(Zx,X,w_wc[:],levels_wwc,linewidth=linw,cmap=cm.jet)
if contour_labels==1:
    plt.clabel(CS_wwc,inline=1,orientation='horizontal',fmt='%1.0f',fontsize=11,colors='k',manual=True)
if contour_labels==2:
    plt.clabel(CS_wwc,inline=1,orientation='horizontal',fmt='%1.0f',fontsize=11,colors='k',manual=False)
    
# plot small circle at the launching point:
plt.plot(wz[0,0],wx[0,0],'ko')

# plot walls, if any:
if n_wall>0:
    plt.plot(z_wall,r_wall,'k',linewidth=linw*2)
plt.savefig('genray_rays_PolFlux_inXZ.png') 
plt.show()


#-------------------------------------------------------------------------------
fig3=plt.figure()  # RAYS and omega/omega_c in cross-section view X-Z
ax = plt.subplot(1, 1, 1)  # rays over |B|
ax.set_aspect(1.0)
plt.axis([Zmin,Zmax,Rmin*1.1,Rmax*1.1])
plt.hold(True)
plt.ylabel('$X$  $(cm)$',fontsize=fnt+4)
plt.xlabel('$Z$  $(cm)$',fontsize=fnt+4)
isp=0 # electrons, by default
if Nsp>1:
    isp=isp_wc-1 # ion (isp_wc is the genray-c index; isp is the python index)
plt.title(r'$Rays$  $and$  $\omega / \omega_c$  ['+isp_name[isp]+']',fontsize=fnt+4,y=1.03)
X,Zx = np.meshgrid(eqdsk_x, eqdsk_z)
X=X*100
Zx=Zx*100    
#if ((Dmax>1000*Dmin) & (Dmin>1e-4)): 
#    # if Bmin is at least 1 Gauss, check range in |B|, 
#    # reduce to avoid peaks near wires (of magnetic coils, if present)
#    Dmin=BvsX_min
#    Dmax=max(BvsX_max,50*BvsX_min)
#    print 'For plots of rays over |B| levels, Bmin and Bmax are adjusted to:'
#    print 'Bmin[T]=',Dmin, '  Bmax[T]=',Dmax
#print shape(D), shape(Z), Dmin, Dmax
#levels=np.arange(Dmin,Dmax,(Dmax-Dmin)/(Ncont-1))

# Three levels of PSI (to mark psimag and psilim)
levels=np.arange(PSImin,PSImax,(PSImax-PSImin)/5)  # 5 lines only
CS=plt.contour( Zr, R,PSI[:],levels,linewidth=linw*3,colors='k')
# Reflection:
CS=plt.contour( Zr,-R,PSI[:],levels,linewidth=linw*3,colors='k')
# Levels of omega/omega_c :
CS_wwc=plt.contour(Zx,X,w_wc[:],levels_wwc,linewidth=linw,cmap=cm.jet)
if contour_labels==1:
    plt.clabel(CS_wwc,inline=1,fmt='%1.0f',fontsize=18,colors='k',manual=True)
if contour_labels==2:
    plt.clabel(CS_wwc,inline=1,fmt='%1.0f',fontsize=18,colors='k',manual=False)
CB=plt.colorbar(orientation='horizontal', shrink=1.0, format='%1.0f') 
if ios=='windows': 
    CB.lines.set_linewidth(5)
CB.set_label(r'$\omega / \omega_c$', size=14)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(wz[i,0:Nm],wx[i,0:Nm],color=col,linewidth=linw)
    # plot arrow for refractive vector (Nx,Nz) at starting point: 
    Nx0=wn_x[i,0]
    Nz0=wn_z[i,0]
    N0=sqrt(Nx0*Nx0+Nz0*Nz0)
    #plt.arrow(wx[i,0],wz[i,0],(Nx0/N0)*arr_len,(Nz0/N0)*arr_len,'->',color=col) 
    # plot small dot at the launching point:
    plt.plot(wz[i,0],wx[i,0],'k.')
    
# plot small circle at the launching point:
plt.plot(wz[0,0],wx[0,0],'ko')
# plot walls, if any:
if n_wall>0:
    plt.plot(z_wall,r_wall,'k',linewidth=linw*2)
plt.savefig('genray_rays_w_wc_inXZ.png') 
plt.show()


#--------------------------------------------------------------------------
fig7=plt.figure()  # B-equilibrium-field along RAYS vs distance(t)
plt.subplot(311) #-------------------------
plt.hold(True)
plt.grid(True)
plt.minorticks_on() # To add minor ticks
plt.tick_params(which='both',  width=1)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4, color='k')
plt.ylabel('$B_Z$  $(G)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],sb_z[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(312) #-------------------------
plt.hold(True)
plt.grid(True)
plt.minorticks_on() # To add minor ticks
plt.tick_params(which='both',  width=1)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4, color='k')
plt.ylabel('$B_R$  $(G)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],sb_r[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(313) #-------------------------
plt.hold(True)
plt.grid(True)
plt.minorticks_on() # To add minor ticks
plt.tick_params(which='both',  width=1)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4, color='k')
plt.ylabel('$|B|$  $(G)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],sbtot[i,0:Nm],color=col,linewidth=linw)
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
plt.savefig('genray_rays_BzBr_p.png') 
show()


#--------------------------------------------------------------------------
fig8=plt.figure()  # E-wave-field along RAYS vs distance(t)
plt.subplot(231) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$|E_X/E|$ $along$ $ray$',y=1.03)
plt.ylim((0.0,1.05))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    Ex_real=cwexde[0,i,0:Nm].copy()
    Ex_imag=cwexde[1,i,0:Nm].copy()
    Ex_abs= sqrt(Ex_real**2 + Ex_imag**2)  #
    plt.plot(ws[i,0:Nm],Ex_abs,color=col,linewidth=linw)  
plt.subplot(232) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$|E_Y/E|$ $along$ $ray$',y=1.03)
plt.ylim((0.0,1.05))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    Ey_real=copy(cweyde[0,i,0:Nm])
    Ey_imag=cweyde[1,i,0:Nm].copy()
    Ey_abs= sqrt(Ey_real**2 + Ey_imag**2)
    plt.plot(ws[i,0:Nm],Ey_abs,color=col,linewidth=linw)  
plt.subplot(233) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$|E_Z/E|$ $along$ $ray$',y=1.03)
plt.ylim((0.0,1.05))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    ereal=cwezde[0,i,0:Nm].copy()
    eimag=cwezde[1,i,0:Nm].copy()
    ea= sqrt(ereal**2 + eimag**2)
    plt.plot(ws[i,0:Nm],ea,color=col,linewidth=linw)
plt.subplot(234) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$Power Flux$ $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],fluxn[i,0:Nm],color=col,linewidth=linw)
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
plt.subplot(236) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title(' $Lin.Damp.$'+' $K_{im}$'+'$(cm^{-1})$') #may include ions, dep.on i_salphal
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    if Nsp>1:  # if ions present, subtract damping from all IONS
        #plt.plot(ws[i,0:Nm],salphal[i,0:Nm]-salphas_all_i[i,0:Nm],color=col,linewidth=linw)
        #NEED TO CHECK: #salphal may include ions, but it dep.on i_salphal values !!!
        #For now, plot salphal, without subtraction of ion (possible) contribution
        plt.plot(ws[i,0:Nm],salphal[i,0:Nm],color=col,linewidth=linw)
    else:        
        plt.plot(ws[i,0:Nm],salphal[i,0:Nm],color=col,linewidth=linw) # electrons only
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')

if (Nsp>1) & (isp_wc>1):  # if ions present, plot Kim for ion species#1
    plt.subplot(235) #-------------------------
    plt.hold(True)
    plt.grid(True)
    plt.title(' $Lin.Damp.(ions)$'+' $K_{im}$')
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
        plt.plot(ws[i,0:Nm],salphas_all_i[i,0:Nm],color=col,linewidth=linw) 
        #                        # ALL IONS, summed up above in this script
                                 
plt.savefig('genray_rays_Ewave_p.png') 
show()


fig8=plt.figure()  # E+ and E- wave-field along RAYS vs distance(t)
plt.subplot(611) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$|E_{+}/E|$')
plt.ylim((0.0,1.05))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    Ex_real=cwexde[0,i,0:Nm].copy()
    Ex_imag=cwexde[1,i,0:Nm].copy()
    Ey_real=cweyde[0,i,0:Nm].copy()
    Ey_imag=cweyde[1,i,0:Nm].copy()
    # Form Eplus as (Ex+i*Ey)/sqrt(2.)  
    Eplus_real= (Ex_real - Ey_imag)/sqrt(2.)
    Eplus_imag= (Ex_imag + Ey_real)/sqrt(2.)
    Eplus_abs= sqrt(Eplus_real**2 + Eplus_imag**2)
    plt.plot(ws[i,0:Nm],Eplus_abs,color=col,linewidth=linw)  
plt.subplot(612) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$|E_{-}/E|$')
plt.ylim((0.0,1.05))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    Ex_real=cwexde[0,i,0:Nm].copy()
    Ex_imag=cwexde[1,i,0:Nm].copy()
    Ey_real=cweyde[0,i,0:Nm].copy()
    Ey_imag=cweyde[1,i,0:Nm].copy()
    # Form Eminus as (Ex-i*Ey)/sqrt(2.)
    Eminus_real= (Ex_real + Ey_imag)/sqrt(2.)
    Eminus_imag= (Ex_imag - Ey_real)/sqrt(2.)
    Eminus_abs= sqrt(Eminus_real**2 + Eminus_imag**2)
    plt.plot(ws[i,0:Nm],Eminus_abs,color=col,linewidth=linw)  
plt.subplot(613) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$|E_{||}/E|$') # WARNING: is Ez/E==cwezde actually parallel to Beq ?
plt.ylim((0.0,1.05))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    ereal=cwezde[0,i,0:Nm].copy()
    eimag=cwezde[1,i,0:Nm].copy()
    ea= sqrt(ereal**2 + eimag**2)
    plt.plot(ws[i,0:Nm],ea,color=col,linewidth=linw)
plt.subplot(614) #-------------------------
plt.hold(True)
isp=0 # electrons, by default
if Nsp>1:
    isp=isp_wc-1 # ion (isp_wc is the genray-c index; isp is the python index)
    mime= mass[isp]/mass[0]        # m_s/m_e ratio
    Z_s= charge[isp]/charge[0]    # Z_s/e
    print 'm_s/m_e=', mime, '  Z_s=',Z_s
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
        Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
        plt.plot(ws[i,0:Nm],f/(28e5*sbtot[i,0:Nm]*Z_s/mime),color=col,linewidth=linw)
plt.subplot(615) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$K_{im}$'+' $(cm^{-1})$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
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
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],wnper[i,0:Nm],color=col,linewidth=linw)
      
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
plt.savefig('genray_rays_Ewave_plus-min_p.png') 
show()    
    
#stop   
   
fig8=plt.figure()  # Z,R, Kperp*rho_Larm, V//res, Ki, along RAYS vs distance(t)
ax=plt.subplot(611) #-------------------------   
plt.hold(True)
plt.grid(True)
plt.ylabel('$Z$'+' $(cm)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],wz[i,0:Nm],color=col,linewidth=linw)  
ax=plt.subplot(612) #-------------------------   
plt.hold(True)
plt.grid(True)
plt.ylabel('$|R|$'+' $(cm)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],wr[i,0:Nm],color=col,linewidth=linw)  
ax=plt.subplot(613) #-------------------------   
plt.hold(True)
plt.grid(True)
plt.ylabel(r'$ K_{\perp} \rho_{Larm} $')
if Nsp>1:
    isp=isp_wc-1 # ion (isp_wc is the genray-c index; isp is the python index)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    ###plt.plot(ws[i,0:Nm],wnpar[i,0:Nm],color=col,linewidth=linw)  
    mime= mass[isp]/mass[0]      # m_s/m_e ratio
    Z_s= charge[isp]/charge[0]    # Z_s/e
    fc= 28e5*sbtot[i,0:Nm]*Z_s/mime  # Hz
    wc_w= fc/f   # omega_c/omega for given species [isp]
    #NII=  wnpar[i,0:Nm]
    Vti=  vthermal[isp,i,0:Nm]  # cm/s
    rho_L= Vti/(fc*2*pi)  # Vthermal/omega_c  # cm
    if Nm>0:
        Kperp_rho= (wnper[i,0:Nm]/wc_w)*(Vti/clight)
        #The above is Nperp*(omega/omega_c)*(Vt/c) == Kperp*rho_thermal
        plt.plot(ws[i,0:Nm],Kperp_rho,color=col,linewidth=linw)
        #plt.plot(ws[i,0:Nm],rho_L,color=col,linewidth=linw)
        print 'mi/me=', mime,  'max of Vti   along ray=', np.max(Vti)
        print 'mi/me=', mime,  'max of rho_L along ray=', np.max(rho_L)
        print 'mi/me=', mime,  'max of wc/w along ray=', np.max(wc_w)
        print 'mi/me=', mime,  'max of Nperp along ray=', np.max(wnper[i,0:Nm])

ax=plt.subplot(614) #-------------------------   
isp=0  # (species 0 is for electrons)
#nres=7 #0
if Nsp>1:
    isp=isp_wc-1 # ion (isp_wc is the genray-c index; isp is the python index)
    #nres=7 # Specify resonance harmonic number such as omega ~ nres*omega_ci
plt.hold(True)
plt.grid(True)
plt.ylabel(r"${V_{||n}}/{V_{T}}$"+'\n'+'$(n=$'+r"$%i$"%(nres) +'$)$')
# To make a new line, use +'\n'+ in the above.
#if nres==0:
#    plt.ylabel(r"${%omega}/{K_{||}V_{T}}$" )
#if nres==1:
#    plt.ylabel(r"${V_{||n=1}}/{V_{T}}$" )
#if nres==2:
#    plt.ylabel(r"${V_{||n=2}}/{V_{T}}$" )
#if nres==3:
#    plt.ylabel(r"${V_{||n=3}}/{V_{T}}$" )
#if nres==4:
#    plt.ylabel(r"${V_{||n=4}}/{V_{T}}$" )
#if nres==5:
#    plt.ylabel(r"${V_{||n=5}}/{V_{T}}$" )
plt.ylim((-10,+10))
#ax.set_yticks([-6,-4,-2,0,2,4,6])
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    if Nm>0:
        mime= mass[isp]/mass[0]      # m_s/m_e ratio
        Z_s= charge[isp]/charge[0]    # Z_s/e
        wc_w= (28e5*sbtot[i,0:Nm]*Z_s/mime)/f   # omega_ci/omega
        NII=  wnpar[i,0:Nm]
        Vti=  vthermal[isp,i,0:Nm]
        VV= (1-nres*wc_w)/(NII*Vti/clight)
        #The above is V||res/Vti  = (omega-nres*omega_c)/(K||*Vti)
        plt.plot(ws[i,0:Nm],VV,color=col,linewidth=linw) 
        
ax=plt.subplot(615) #-------------------------   
nres=nres+1 # Specify resonance harmonic number such as omega ~ nres*omega_ci
plt.hold(True)
plt.grid(True)
plt.ylabel(r"${V_{||n}}/{V_{T}}$"+'\n'+'$(n=$'+r"$%i$"%(nres) +'$)$')
plt.ylim((-10,+10))
#ax.set_yticks([-6,-4,-2,0,2,4,6])
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    if Nm>0:
        mime= mass[isp]/mass[0]      # m_s/m_e ratio
        Z_s= charge[isp]/charge[0]    # Z_s/e
        wc_w= (28e5*sbtot[i,0:Nm]*Z_s/mime)/f   # omega_ci/omega
        NII=  wnpar[i,0:Nm]
        Vti=  vthermal[isp,i,0:Nm]
        VV= (1-nres*wc_w)/(NII*Vti/clight)
        #The above is V||res/Vti  = (omega-nres*omega_c)/(K||*Vti)
        plt.plot(ws[i,0:Nm],VV,color=col,linewidth=linw) 

ax=plt.subplot(616) #-------------------------
plt.hold(True)
plt.grid(True)
plt.ylabel('$K_{im}$'+'$(cm^{-1})$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],salphal[i,0:Nm],color=col,linewidth=linw)  # e+i (all)
    
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
plt.savefig('genray_rays_VIIres_p.png') 
show()
#stop


fig8=plt.figure()  #  N// along RAYS vs distance(t)
ax=plt.subplot(312) #-------------------------   
plt.hold(True)
plt.grid(True)
plt.title('$Parallel$ $refraction$ $index$ $n_{||}$',fontsize=fnt+4,y=1.03)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    if Nm>0:
        plt.plot(ws[i,0:Nm],wnpar[i,0:Nm],color=col,linewidth=linw)  
plt.xlabel('$distance$ $along$ $ray$  $(cm)$',fontsize=fnt+4)
plt.savefig('genray_rays_NII_p.png') 
show()   

if Nsp>1:
    isp=isp_wc-1 # ion (isp_wc is the genray-c index; isp is the python index)
    fig8=plt.figure()  #  (omega/Kper)/Vti = (1/nper)/(Vti/c)      along RAYS
    ax=plt.subplot(312) #-------------------------   
    plt.hold(True)
    plt.grid(True)
    plt.title(r'$({\omega}/K_{\perp})/V_{ti}$'+'$ ['+isp_name[isp]+']$',fontsize=fnt+4,y=1.03)
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
        if Nm>0:
            Vti=  vthermal[isp,i,0:Nm]  # cm/s
            #  (omega/Kper)/Vti = (1/nper)/(Vti/c)
            wKperVti= (1/wnper[i,0:Nm])/(Vti/clight)
            plt.plot(ws[i,0:Nm],wKperVti,color=col,linewidth=linw)  
    plt.xlabel('$distance$ $along$ $ray$  $(cm)$',fontsize=fnt+4)
    plt.savefig('genray_rays_wKperVti_p.png') 
    show()   
#stop

#--------------------------------------------------------------------------
fig81=plt.figure()  # Vgroup/c along RAYS vs distance(t)
plt.subplot(231) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$V_{group,x}/c$  $along$ $ray$')
Vgr_mn= np.amin(vgr_x[:])
Vgr_mx= np.amax(vgr_x[:])
Vgr_mn=max(Vgr_mn,-1.0)
Vgr_mx=min(Vgr_mx, 1.0)
plt.ylim((Vgr_mn-1e-4,Vgr_mx+1e-4))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],vgr_x[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(234) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$V_{group,y}/c$  $along$ $ray$')
Vgr_mn= np.amin(vgr_y[:])
Vgr_mx= np.amax(vgr_y[:])
Vgr_mn=max(Vgr_mn,-1.0)
Vgr_mx=min(Vgr_mx, 1.0)
plt.ylim((Vgr_mn-1e-4,Vgr_mx+1e-4))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],vgr_y[i,0:Nm],color=col,linewidth=linw)  
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
plt.subplot(233) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$V_{group,z}/c$  $along$ $ray$')
Vgr_mn= np.amin(vgr_z[:])
Vgr_mx= np.amax(vgr_z[:])
Vgr_mn=max(Vgr_mn,-1.0)
Vgr_mx=min(Vgr_mx, 1.0)
plt.ylim((Vgr_mn-1e-4,Vgr_mx+1e-4))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],vgr_z[i,0:Nm],color=col,linewidth=linw)
plt.subplot(236) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$|V_{group}/c|$')
vgr=sqrt(vgr_x[:]**2 +vgr_y[:]**2 +vgr_z[:]**2)
#print 'vgr shape:', vgr.shape
Vgr_mn= np.amin(vgr[:])
Vgr_mx= np.amax(vgr[:])
Vgr_mn=max(Vgr_mn, 0.0)
Vgr_mx=min(Vgr_mx, 1.0)
plt.ylim((Vgr_mn,Vgr_mx))
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],vgr[i,0:Nm],color=col,linewidth=linw)
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
plt.savefig('genray_rays_Vgroup_p.png') 
show()

#stop

#--------------------------------------------------------------------------
fig5=plt.figure()  # wce/w and (wpe/w)^2 along RAYS vs distance
f=freqcy
q2m= charge[0]**2/mass[0]
plt.subplot(231) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$\omega/\omega_{LH}$  $along$ $rays$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    #---------------------
    wpe_w_2= 806.2e5*sene[i,0:Nm]*q2m/f**2  # (omega_{pe}/omega)^2
    wpe_w= sqrt(wpe_w_2)       # omega_{pe}/omega  along ray
    wce_w=28e5*sbtot[i,0:Nm]/f # omega_{ce}/omega  along ray
    if Nsp>1:  # add ion species to calculate omega_LH
        isp=1
        msme= np.asscalar(mass[isp]/mass[0])      # ms/me ratio
        Zs= np.asscalar(charge[isp]/charge[0])    # Zs/e
        #print msme,Zs, Zs/msme
        wpi_w= wpe_w*(Zs/sqrt(msme)) #  == omega_pi/omega    BASED on ASSUMPTION ni=ne
        #print wpe_w/wpi_w
        wci_w= (Zs/msme)*wce_w    #  == omega_ci/omega
        # Define omega_LH/omega. When wLH_w=1, eps_xx becomes 0 (resonance)
        wLH_w=sqrt( ((wce_w*wci_w)**2 +(wci_w*wpe_w)**2 +(wpi_w*wce_w)**2)/(wce_w**2+wpe_w_2) )    
        #---------------------
        plt.plot(ws[i,0:Nm],1/wLH_w,color=col,linewidth=linw)
#plt.xlabel('$distance$ $along$ $ray$  $(cm)$')

plt.subplot(232) #-------------------------
txt=r"$f_{MHz}=$"+r"$%1.1f$" %(f*1e-6)
plt.title(txt)   # 
plt.axis('off')

plt.subplot(233) #-------------------------


plt.hold(True)
plt.ylabel('$Power$ $in$ $ray$ $channel$  $(kW)$' )
#plt.ylabel('$ delpwr $  $ (kW) $' )
ws_max=0.
plt.grid(True)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    ws_max1=np.amax(ws[:,Nm])
    ws_max=np.amax([ws_max1,ws_max])
    plt.plot(ws[i,0:Nm],delpwr[i,0:Nm],color=col,linewidth=linw)
#plt.axis([0.,1.05*ws_max,0.,1.05*np.amax(delpwr[:,0])] )
plt.ylim((0.,1.05*np.amax(delpwr[:,0])))

#--------------------------------------------
plt.subplot(235) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$(\omega_{pe}/\omega)^2$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    #print 806.2e5*sene[0,0]*q2m/f**2, sene[0,0], q2m
    wpe_w_2= 806.2e5*sene[i,0:Nm]*q2m/f**2  # (omega_{pe}/omega)^2
    wpe_w= sqrt(wpe_w_2)  # omega_{pe}/omega along ray
    plt.plot(ws[i,0:Nm],wpe_w_2,color=col,linewidth=linw)
#plt.xlabel('$distance$ $along$ $ray$  $(cm)$')


#--------------------------------------------
plt.subplot(234) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$\omega/\omega_{ce}$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    wce_w= f/(28e5*sbtot[i,0:Nm]) # omega/omega_ce
    plt.plot(ws[i,0:Nm],wce_w,color=col,linewidth=linw)
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
#--------------------------------------------
if Nsp>1:
    isp=isp_wc-1 # ion (isp_wc is the genray-c index; isp is the python index)
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

#stop

#--------------------------------------------------------------------------
fig5=plt.figure()  # Refractive indices along RAYS vs R(t)
ax=plt.subplot(231) #-------------------------
plt.hold(True)
plt.grid(True)
#plt.title('$N_R$ $vs$ $R$') # N_R ws R
plt.title('$N_X$ $vs$ $X$')  # N_X ws x
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    #plt.plot(wr[i,0:Nm],wn_r[i,0:Nm],color=col,linewidth=linw)
    plt.plot(wx[i,0:Nm],wn_x[i,0:Nm],color=col,linewidth=linw)    
plt.subplot(232) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$N_Z$ $vs$ $R$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(wr[i,0:Nm],wn_z[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(233) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title(r'      $N_\phi$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(wr[i,0:Nm],wn_phi[i,0:Nm],color=col,linewidth=linw)
plt.subplot(234) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title(r'$N_{\perp}$'+' $along$ $ray$')
plt.xlabel('$R$ $along$ $ray$  $(cm)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(wr[i,0:Nm],wnper[i,0:Nm],color=col,linewidth=linw)
plt.subplot(235) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title(r'$N_{||}$'+' $along$ $ray$')
plt.xlabel('$R$ $along$ $ray$  $(cm)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(wr[i,0:Nm],wnpar[i,0:Nm],color=col,linewidth=linw)
plt.subplot(236) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$   Lin.Damping$'+' $N_{im}$')
plt.xlabel('$R$ $along$ $ray$  $(cm)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(wr[i,0:Nm],salphal[i,0:Nm]*clight/omega,color=col,linewidth=linw)
    # salphal is K_im [1/cm]; Includes electrons and ions (if present)
plt.savefig('genray_rays_refr-index_R.png') 
show()

#stop




#--------------------------------------------------------------------------
fig6=plt.figure()  # Refractive indices along RAYS vs distance(t)
plt.subplot(231) #-------------------------
plt.hold(True)
plt.grid(True)
#plt.title('$N_R$ $vs$ $R$') # N_R ws R
plt.title('$N_X$ $vs$ $X$')  # N_X ws x
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    #plt.plot(ws[i,0:Nm],wn_r[i,0:Nm],color=col,linewidth=linw)
    plt.plot(ws[i,0:Nm],wn_x[i,0:Nm],color=col,linewidth=linw)    
plt.subplot(232) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title('$N_Z$ $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],wn_z[i,0:Nm],color=col,linewidth=linw)  
plt.subplot(233) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title(r'     $N_\phi$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],wn_phi[i,0:Nm],color=col,linewidth=linw)
plt.subplot(234) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title(r'$N_{\perp}$'+' $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],wnper[i,0:Nm],color=col,linewidth=linw)
plt.subplot(235) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title(r'$N_{||}$'+' $along$ $ray$')
plt.xlabel('$distance$ $along$ $ray$  $(cm)$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],wnpar[i,0:Nm],color=col,linewidth=linw)
plt.subplot(236) #-------------------------
plt.hold(True)
plt.grid(True)
plt.title(r'$|N|$'+' $along$ $ray$')
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    plt.plot(ws[i,0:Nm],sqrt(wnpar[i,0:Nm]**2+wnper[i,0:Nm]**2),color=col,linewidth=linw)
plt.savefig('genray_rays_refr-index_p.png') 
show()





#--------------------------------------------------------------------------
if i_ox==2:  # only for O-X transmission case
    fig8=plt.figure()
    plt.subplot(221) #-------------------------
    plt.hold(True)
    plt.grid(True)
    #plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
    plt.ylabel('$1=converted;$   $0-not$')
    plt.title('Indicator of OX conversion')
    plt.ylim((0.,1.1))
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        plt.plot(wnpar[i,0],i_ox_conversion[i],'o',color=col)
    plt.subplot(222) #-------------------------
    plt.hold(True)
    plt.grid(True)
    plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
    #text(wnpar[0,0],1.02,'   $pwr(Xmode)/pwr(Omode)$')
    plt.title('Transm.coef=P(X)/P(O)')
    
    plt.ylim((0.,1.1))
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        plt.plot(wnpar[i,0],transm_ox[i],'o',color=col)
    plt.subplot(223) #-------------------------
    plt.hold(True)
    plt.grid(True)
    plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
    plt.ylabel('$N_{||}$ $at$ $OX$ $transm.$')
    y_mn= min(0.,1.1*np.amin(cnpar_ox[:]))
    y_mx= max(0.,1.1*np.amax(cnpar_ox[:]))
    plt.ylim((y_mn,y_mx))
    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
        if remainder(i,6)==0: col='b'
        if remainder(i,6)==1: col='g'
        if remainder(i,6)==2: col='r'
        if remainder(i,6)==3: col='c'    
        if remainder(i,6)==4: col='m' 
        if remainder(i,6)==5: col='k'  
        plt.plot(wnpar[i,0],cnpar_ox[i],'o',color=col)
#    plt.subplot(224) #-------------------------
#    plt.hold(True)
#    plt.xlabel('$N_{||}$ $at$ $starting$ $point$')
#    plt.ylabel('$N$ $along$ $B$x$grad(\Psi)$')
#    for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
#        if remainder(i,6)==0: col='b'
#        if remainder(i,6)==1: col='g'
#        if remainder(i,6)==2: col='r'
#        if remainder(i,6)==3: col='c'    
#        if remainder(i,6)==4: col='m' 
#        if remainder(i,6)==5: col='k'  
#        plt.plot(wnpar[i,0],cn_b_gradpsi[i],'o',color=col)
    
    plt.savefig('genray_rays_OX_transmission.png') 
    show()



#--------------------------------------------------------------------------
fig7=plt.figure()   # Power profiles

NR=rho_bin_center[:].size
# find radial index corresponding to rho~1 (approximately)
for ir in range(1,NR,1):
    if np.abs(rho_bin_center[ir]-1.0)<0.1 :
        NR1= ir
    else:
        NR1=NR

if rhomax==0:
    rhomax=1.05*np.amax(rho_bin_center)
rhomax=min(rhomax, rho_bin_center[NR-1])
print 'rho_bin_center: min/max',np.amin(rho_bin_center),np.amax(rho_bin_center)
print 'rho_bin: min/max', np.amin(rho_bin),np.amax(rho_bin)
drho= rho_bin_center[1]-rho_bin_center[0]  # step in rho grid
irmx= int(rhomax/drho) # index corresponding to rhomax
irmx=min(irmx,NR-1)
print 'rhomax for plots:', rhomax
print 'irmx, rho_bin_center[irmx]', irmx, rho_bin_center[irmx]

plt.subplot(221) #-------------------------
txt="$p_{e}$"+"$=$%3.3f" %(powtot_e) +" $kW$"
plt.title(txt,y=1.03)
plt.hold(True)
plt.grid(True)
plt.xlim((0.,rhomax))
plt.ylim((0.,1.05*np.amax(powden_e[0:irmx])+1.e-8))
#plt.axis([0.,rhomax,0.,1.05*np.amax(powden_e[0:irmx])+1.e-8])
plt.plot(rho_bin_center[0:irmx],powden_e[0:irmx],'r.',linewidth=linw)
plt.plot(rho_bin_center[0:irmx],powden_e[0:irmx],'k',linewidth=linw)
print 'binvol[NR-1], binvol[0] (cm^3)', binvol[NR-1],binvol[0]
if binvol[NR1-1]>binvol[0]:
    plt.ylabel('$W/cm^3$')
if binvol[NR1-1]==binvol[0]:
    plt.ylabel('$W$ $per$ $bin$')
    
plt.subplot(222) #-------------------------
plt.hold(True)
plt.grid(True)
txt="$p_{i}$"+"$=$%3.3f" %(powtot_i) +" $kW$"
plt.title(txt,y=1.03)
plt.xlim((0.,rhomax))
plt.ylim((0.,1.05*np.amax(powden_i[0:irmx])+1.e-8))
#plt.axis([0.,rhomax,0.,1.05*np.amax(powden_i[0:irmx])+1.e-8])
plt.plot(rho_bin_center[0:irmx],powden_i[0:irmx],'r.',linewidth=linw)
plt.plot(rho_bin_center[0:irmx],powden_i[0:irmx],'k',linewidth=linw)

plt.subplot(223) #-------------------------
plt.hold(True)
plt.grid(True)
plt.xlabel(r'$\rho$')
txt="$p_{cl}$"+"$=$%3.3f" %(powtot_cl) +" $kW$"
plt.title(txt,verticalalignment='center',y=1.03)
if binvol[NR1-1]>binvol[0]:
    plt.ylabel('$W/cm^3$')
if binvol[NR1-1]==binvol[0]:
    plt.ylabel('$W$ $per$ $bin$')
if n_powden_cl>0:
    plt.xlim((0.,rhomax))
    plt.ylim((0.,1.05*np.amax(powden_cl[0:irmx])+1.e-8))
    #plt.axis([0.,rhomax,0.,1.05*np.amax(powden_cl[0:irmx])+1.e-8])
    plt.plot(rho_bin_center[0:irmx],powden_cl[0:irmx],'r.',linewidth=linw)
    plt.plot(rho_bin_center[0:irmx],powden_cl[0:irmx],'k',linewidth=linw)
    
plt.subplot(224) #-------------------------
plt.hold(True)
plt.grid(True)
plt.xlabel(r'$\rho$')
txt="$p_{total}$"+"$=$%3.3f" %(power_total) +" $kW$"
plt.title(txt,verticalalignment='center',y=1.03)
plt.xlim((0.,rhomax))
plt.ylim((0.,1.05*np.amax(powden[0:irmx])+1.e-8))
#plt.axis([0.,rhomax,0.,1.05*np.amax(powden[0:irmx])+1.e-8])
plt.plot(rho_bin_center[0:irmx],powden[0:irmx],'r.',linewidth=linw)
plt.plot(rho_bin_center[0:irmx],powden[0:irmx],'k',linewidth=linw)
savefig('genray_profiles_power.png')
show() 


print 'Sum over rays:'
print 'integral(powden*binvol)=    ', sum(powden*binvol),   ' W'
print 'integral(powden_e*binvol)=  ', sum(powden_e*binvol), ' W'
print 'integral(powden_i*binvol)=  ', sum(powden_i*binvol), ' W'
print 'integral(powden_cl*binvol)= ', sum(powden_cl*binvol),' W'
print 'Note: check plots of binvol: Could be set to 1.0 if no flux surfaces'
print 'ABSORBED/LOST-in-plasma:  power_total =   ',power_total,    ' kW'
print 'INJECTED (STARTED power): power_inj_total=',power_inj_total,' kW'
print '==========================================================='
#--------------------------------------------------------------------------
#fig1=plt.figure()   # Current profiles
#plt.subplot(221) #-------------------------
#title('Current Profiles  '+'$(A/cm^2)$')
#plt.hold(True)
#plt.grid(True)
#plt.ylabel('$<j_{||}>$')
#plt.plot(rho_bin_center,s_cur_den_parallel,'r',linewidth=linw)
#plt.subplot(222) #-------------------------
#plt.hold(True)
#plt.grid(True)
#plt.title('<j.B>'+'$/B_0$')
#plt.plot(rho_bin_center,s_cur_den_onetwo,'r',linewidth=linw)
#plt.subplot(223) #-------------------------
#plt.hold(True)
#plt.grid(True)
#plt.xlabel(r'$\rho$')
#plt.title('$<j_{||}>f<1/r^2>/(<B><1/r>)$',verticalalignment='center')
#plt.plot(rho_bin_center,s_cur_den_toroidal,'r',linewidth=linw)
#plt.subplot(224) #-------------------------
#plt.hold(True)
#plt.grid(True)
#plt.xlabel(r'$\rho$')
#plt.title('$<j_{||}>B_{pol}($'+r'$\theta$'+'$=0)/<B>$',verticalalignment='center')
#plt.plot(rho_bin_center,s_cur_den_poloidal,'r',linewidth=linw)
#savefig('genray_profiles_J.png')
#show() 



#--------------------------------------------------------------------------
fig10=plt.figure()   # Power profiles over (Z,R) grid and rays over (Z,X) plane

Rgr,Zgr = np.meshgrid(Rgrid,Zgrid) # 2D grids [cm]
Nc_pwr=10 # number of contour levels in power(Z,R) plots

# The 2D plots of absorbed power in (Z,R) plane are based on data accumulated 
# over R=sqrt(X^2 +Y^2) radial grid. If ray is propagating 
# at negative X and/or negative Y, the value of R is still positive.
# In general, the path of rays in (Z,X) plane may not coincide 
# with power pattern in (Z,R) grid because of such mapping.
# For convenience, the power deposition P(Z,R) can be also plotted
# as mirror reflection over the Z axis, i.e. as P(Z,-R),
# or as the combination of the two, i.e. P(Z,R) & P(Z,-R).
# i_plot_powerRZ= 0 ==> Do not make such plots at all.
# i_plot_powerRZ=+1 ==> Plot P(Z, R) (where R is .ge.0).
# i_plot_powerRZ=-1 ==> Plot mirror image P(Z,-R) 
# i_plot_powerRZ= 2 ==> Plot the combination P(Z,R)&P(Z,-R) 

ax=plt.subplot(111) #-------------------------
ax.set_aspect(1.0)
plt.axis([Zmin,Zmax,-Rmax,Rmax])
txt="$p_{e}$"+"$=$%3.3f" %(sum_spwr_rz_e) +" $kW$"
plt.title(txt,y=1.03)
plt.ylabel('$X$  $(cm)$')
plt.xlabel('$Z$  $(cm)$')
if sum_spwr_rz_e>0:
    if i_plot_powerRZ==+1 or i_plot_powerRZ==2 :
        CSP=plt.contour(Zgr,+Rgr,spwr_rz_e[:],Nc_pwr,cmap=plt.cm.jet)
    # Mirror Reflection:
    if i_plot_powerRZ==-1 or i_plot_powerRZ==2 :
        CSP=plt.contour(Zgr,-Rgr,spwr_rz_e[:],Nc_pwr,cmap=plt.cm.jet)

levels=np.arange(PSImin,PSImax,(PSImax-PSImin)/5)
CS=plt.contour( Zr, R,PSI[:],levels,linewidth=linw,colors='k')
# Reflection:
CS=plt.contour( Zr,-R,PSI[:],levels,linewidth=linw,colors='k')

# Plot levels of omega/omega_c:
#CS_wwc=plt.contour(Zx,X,w_wc[:],levels_wwc,linewidth=linw,cmap=cm.jet)
#if contour_labels==1:
#    plt.clabel(CS_wwc,inline=1,orientation='horizontal',fmt='%1.0f',fontsize=11,colors='k',manual=True)
#if contour_labels==2:
#    plt.clabel(CS_wwc,inline=1,orientation='horizontal',fmt='%1.0f',fontsize=11,colors='k',manual=False)

# plot small circle at the launching point:
plt.plot(wz[0,0], wx[0,0],'ko')
#plt.plot(wz[0,0],-wr[0,0],'ko')
# plot walls, if any:
if n_wall>0:
    plt.plot(z_wall,r_wall,'k',linewidth=linw*2)    
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    # Rays projected to (Z,X) plane:
    plt.plot(wz[i,0:Nm],wx[i,0:Nm],color=col,linewidth=linw)
    # plot small dot at the launching point:
    plt.plot(wz[i,0],wx[i,0],'k.')
    
savefig('genray_profiles_Pe_RZ.png')
show() 

fig10=plt.figure()   # Power profiles Pi over R,Z grid
ax=plt.subplot(111) #-------------------------
ax.set_aspect(1.0)
plt.axis([Zmin,Zmax,-Rmax,Rmax])
txt="$p_{i}$"+"$=$%3.3f" %(sum_spwr_rz_i) +" $kW$"
plt.title(txt,y=1.03)
plt.ylabel('$X$  $(cm)$')
plt.xlabel('$Z$  $(cm)$')
if sum_spwr_rz_i>0:
    if i_plot_powerRZ==+1 or i_plot_powerRZ==2 :
        CSP=plt.contour(Zgr,+Rgr,spwr_rz_i[:],Nc_pwr,cmap=plt.cm.jet)
    #Reflection:
    if i_plot_powerRZ==-1 or i_plot_powerRZ==2 :
        CSP=plt.contour(Zgr,-Rgr,spwr_rz_i[:],Nc_pwr,cmap=plt.cm.jet)
#CBP=plt.colorbar(orientation='vertical', shrink=0.4,format='%.1e') # colorbar
levels=np.arange(PSImin,PSImax,(PSImax-PSImin)/5)
CS=plt.contour( Zr, R,PSI[:],levels,linewidth=linw,colors='k')
CS=plt.contour( Zr,-R,PSI[:],levels,linewidth=linw,colors='k')
# Plot levels of omega/omega_c:
CS_wwc=plt.contour(Zx,X,w_wc[:],levels_wwc,linewidth=linw,cmap=cm.jet)
if contour_labels==1:
    plt.clabel(CS_wwc,inline=1,orientation='horizontal',fmt='%1.0f',fontsize=11,colors='k',manual=True)
if contour_labels==2:
    plt.clabel(CS_wwc,inline=1,orientation='horizontal',fmt='%1.0f',fontsize=11,colors='k',manual=False)

# plot small circle at the launching point:
plt.plot(wz[0,0], wx[0,0],'ko')
#plt.plot(wz[0,0],-wr[0,0],'ko')
# plot walls, if any:
if n_wall>0:
    plt.plot(z_wall,r_wall,'k',linewidth=linw*2)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    # Rays projected to (R,Z) plane:
    plt.plot(wz[i,0:Nm],wx[i,0:Nm],color=col,linewidth=linw)
    # plot small dot at the launching point:
    plt.plot(wz[i,0],wx[i,0],'k.')
    
savefig('genray_profiles_Pi_RZ.png')
show() 

fig10=plt.figure()   # Power profiles Pcl over R,Z grid
ax=plt.subplot(111) #-------------------------
ax.set_aspect(1.0)
plt.axis([Zmin,Zmax,-Rmax,Rmax,])
txt="$p_{collisional}$"+"$=$%3.3f" %(sum_spwr_rz_cl) +" $kW$"
plt.title(txt,y=1.03)
plt.ylabel('$X$  $(cm)$')
plt.xlabel('$Z$  $(cm)$')
if sum_spwr_rz_cl>0:
    if i_plot_powerRZ==+1 or i_plot_powerRZ==2 :
        CSP=plt.contour(Zgr, Rgr,spwr_rz_cl[:],Nc_pwr,cmap=plt.cm.jet)
    #Mirror Reflection:
    if i_plot_powerRZ==-1 or i_plot_powerRZ==2 :
        CSP=plt.contour(Zgr,-Rgr,spwr_rz_cl[:],Nc_pwr,cmap=plt.cm.jet)
#CBP=plt.colorbar(orientation='vertical', shrink=0.4,format='%.1e') # colorbar
levels=np.arange(PSImin,PSImax,(PSImax-PSImin)/5)
CS=plt.contour( Zr, R,PSI[:],levels,linewidth=linw,colors='k')
CS=plt.contour( Zr,-R,PSI[:],levels,linewidth=linw,colors='k')
# Plot levels of omega/omega_c:
#CS_wwc=plt.contour(Zx,X,w_wc[:],levels_wwc,linewidth=linw,cmap=cm.jet)
#if contour_labels==1:
#    plt.clabel(CS_wwc,inline=1,orientation='horizontal',fmt='%1.0f',fontsize=11,colors='k',manual=True)
#if contour_labels==2:
#    plt.clabel(CS_wwc,inline=1,orientation='horizontal',fmt='%1.0f',fontsize=11,colors='k',manual=False)

# plot small circle at the launching point:
plt.plot(wz[0,0], wx[0,0],'ko')
#plt.plot(wz[0,0],-wr[0,0],'ko')
# plot walls, if any:
if n_wall>0:
    plt.plot(z_wall,r_wall,'k',linewidth=linw*2)
for i in range(0,Nrays,1):  # i goes from 0 to Nrays-1
    if remainder(i,6)==0: col='b'
    if remainder(i,6)==1: col='g'
    if remainder(i,6)==2: col='r'
    if remainder(i,6)==3: col='c'    
    if remainder(i,6)==4: col='m' 
    if remainder(i,6)==5: col='k'  
    Nm= np.asscalar(nrayelt[i])-1 # max number of points along a ray
    # Rays projected to (Z,X) plane:
    plt.plot(wz[i,0:Nm],wx[i,0:Nm],color=col,linewidth=linw)  
    # plot small dot at the launching point:
    plt.plot(wz[i,0],wx[i,0],'k.')
    
if sb_z[0,0]>0:   # Bz>0  plot marker pointing ->
    plt.plot(Zmax*0.9,Rmax*0.98,'r>',markersize=10)    
    text(Zmax*0.9,Rmax*1.1,r'$\rightarrow$ $B_Z$')
if sb_z[0,0]<0:   # Bz<0  plot marker pointing <-
    plt.plot(Zmax*0.9,Rmax*0.98,'r<',markersize=10)    
    text(Zmax*0.9,Rmax*1.1,r'$\leftarrow$ $B_Z$')
    
savefig('genray_profiles_Pcl_RZ.png')
show() 



dat.close() # close genray.nc
         
elapsed_time = time.time() - e0
cpu_time = time.clock() - c0
print 'elapsed and cpu time since start (sec.) =', elapsed_time, cpu_time
print 'FINISHED'