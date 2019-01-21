# 3dcore_v3.py
# Runs the 3DCORE CME flux rope model and makes synthetic in situ magnetic field, 
# plasma speed and Dst data at specified planets and spacecraft.
# This is the only program in this package.
# The publication describing this method and showing 
# first results is Moestl et al. 2018 Space Weather 
# https://arxiv.org/abs/1710.00587
# Author: C. Moestl, IWF Graz, February 2015 - July 2018

#install python 3.5 anaconda for all packages, and add sunpy and seaborn

#This work is published under the MIT LICENSE
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
#PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
#FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
#TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
#OTHER DEALINGS IN THE SOFTWARE.


from scipy import stats
import scipy.io
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LightSource
import numpy as np
import sunpy.sun
import sunpy.time
import time
import pickle
import os
import math
from mpl_toolkits.mplot3d import axes3d
from matplotlib.colors import LightSource
import seaborn as sns
from matplotlib import cm
import pdb
import copy



###########  get input parameters from file

#choose which one:

#for 1st paper targets are Earth and MESSENGER

#Earth
#simulation parameter file
inputfilename='input_files/3DCORE_init_July_2013_paper_Earth.txt'
#directory of output plots, movies, for this event
outputdirectory='output_files_Earth_L1'

#for other inclinations:
#inputfilename='input_files/3DCORE_init_July_2013_paper_Earth_inc232.txt'
#inputfilename='input_files/3DCORE_init_July_2013_paper_Earth_inc272.txt'


#MESSENGER
#inputfilename='input_files/3DCORE_init_July_2013_paper_MESSENGER.txt'
#outputdirectory='output_files_MESSENGER'














############################ procedures
 
def sphere2cart(r, phi, theta):
    x = r*np.cos(theta)*np.cos(phi)
    y = r*np.cos(theta)*np.sin(phi)
    z = r*np.sin(theta)
    return (x, y, z) 

def getpositions(filename):  
    pos=scipy.io.readsav(filename)  
    print
    print('positions file:', filename) 
    return pos

   
def getcat(filename):
  print('reading CAT '+filename)
  cat=scipy.io.readsav(filename, verbose='true')  
  print('done reading CAT')
  return cat  
  
def decode_array(bytearrin):
  #for decoding the strings from the IDL .sav file to a list of python strings, not bytes 
  #make list of python lists with arbitrary length
  bytearrout= ['' for x in range(len(bytearrin))]
  for i in range(0,len(bytearrin)-1):
   bytearrout[i]=bytearrin[i].decode()
  #has to be np array so to be used with numpy "where"
  bytearrout=np.array(bytearrout)
  return bytearrout   
  
def time_to_num_cat(time_in):  
  #for time conversion  
  #for all catalogues
  #time_in is the time in format: 2007-11-17T07:20:00
  #for times help see: 
  #http://docs.sunpy.org/en/latest/guide/time.html
  #http://matplotlib.org/examples/pylab_examples/date_demo2.html
  
  j=0
  #time_str=np.empty(np.size(time_in),dtype='S19')
  time_str= ['' for x in range(len(time_in))]
  #=np.chararray(np.size(time_in),itemsize=19)
  time_num=np.zeros(np.size(time_in))
 
  for i in time_in:
   #convert from bytes (output of scipy.readsav) to string
   time_str[j]=time_in[j][0:16].decode()+':00'
   year=int(time_str[j][0:4])
   #convert time to sunpy friendly time and to matplotlibdatetime
   #only for valid times so 9999 in year is not converted
   #pdb.set_trace()
   if year < 2100:
    	  time_num[j]=mdates.date2num(sunpy.time.parse_time(time_str[j]))
   j=j+1  
   #the date format in matplotlib is e.g. 735202.67569444
   #this is time in days since 0001-01-01 UTC, plus 1.
  return time_num

    
def eltor2cart(cphi, cpsi, cr, aspect):
    #converts ellipse torus to cartesian coordinates 
    #(Sun in Origin, X towards Earth, Y towards solar west)
    #for ellipse torus see: http://mathworld.wolfram.com/EllipticTorus.html
    #a is the r value (minor radius of a normal torus
    #rho as distance from center, and r as major axis a
    a=cr
    b=cr/aspect     #aspect ratio = a/b 

    #auxiliary function which makes torus size from small (Sun) to wide  (apex)
    C=1    
    f=C*np.sin(cpsi/2)          
    #first for B set f to 1
    #f=1
    #these are the normal ellipse torus equations with f
    x = -(rho0+a*f*np.cos(cphi))*np.cos(cpsi)  #negative sign gives correct torus configuration 
    y = (rho0+a*f*np.cos(cphi))*np.sin(cpsi) 
    z = b*f*np.sin(cphi)
    return (x, y, z)     
 
    
def rotate(xi,yi,zi,rot_angle,kx,ky,kz):
    #rotate a point about an axis according to Euler-Rodrigues    
    #rot_angle is in degrees, vector k is the rotation axis
        
    #coefficients
    a=np.cos(np.radians(rot_angle)/2)
    b=kx*np.sin(np.radians(rot_angle)/2)
    c=ky*np.sin(np.radians(rot_angle)/2)
    d=kz*np.sin(np.radians(rot_angle)/2)
    
    #rotate point    
    xo=(a**2+b**2-c**2-d**2)*xi+2*(b*c+a*d)*yi+2*(b*d-a*c)*zi
    yo=2*(b*c-a*d)*xi+(a**2+c**2-b**2-d**2)*yi+2*(c*d+a*b)*zi
    zo=2*(b*d+a*c)*xi+2*(c*d-a*b)*yi+(a**2+d**2-b**2-c**2)*zi
    
    return (xo,yo,zo)

   
def unifield2cart(xb,yb,zb,cr,b0,tau,phit,psit,handedness):  
  
  #calculate a uniform twist magnetic field at every point xb,yb,zb of the cloud
  #and convert to a cartesian cloud coordinate system
  
  #see uniform twist Farrugia et al. 1999 and Hu et al. 2014   
  #azimuthal is the phi component
  bphic=handedness*b0*tau*cr/(1+(tau**2)*(cr**2))
  #axial is the psi component
  bpsic=b0/(1+(tau**2)*(cr**2))
  #no radial component
  brc=0
  
  #get unit vectors for torus with auxiliary function by derivatives http://mathworld.wolfram.com/SphericalCoordinates.html
  #they are given for torus coordinates here http://en.m.wikipedia.org/wiki/Toroidal_and_poloidal  - here theta= our phi, their zeta= our psi
  
  #unit vectors for torus:  
  #er=[cos(phi)*cos(psi), cos(phi)*sin(psi), sin(phi)]
  #ephi=[-sin(phi)*cos(psi), -sin(phi)*sin(psi), cos(phi)]  
  #epsi=[-sin(psi), cos(psi), 0]
  
  #this equation is simply the result of dot products with the unit vectors above
  #Bvec=brc*er+Bphi*ephi+Bpsi*epsi=bx*ex+by*ey+bz*ez  this is multplied with e.g. ex to find bx
  #bx=brc* er dot ex + Bphi * ephi dot ex + Bpsi* epsi dot ex
  bxcloud=brc*np.cos(phit)*np.cos(psit)+bphic*(-np.sin(phit)*np.cos(psit))+bpsic*(-np.sin(psit))
  #by=brc* er dot ey + Bphi * ephi dot ey + Bpsi* epsi dot ey
  bycloud=brc*np.cos(phit)*np.sin(psit)+bphic*(-np.sin(phit)*np.sin(psit))+bpsic*(np.cos(psit))
  #bz=brc* er dot ez + Bphi * ephi dot ez + Bpsi* epsi dot ez
  bzcloud=brc*np.sin(phit)+bphic*np.cos(phit)

  
  return (bxcloud,bycloud,bzcloud)
  
  
 
def make_3dcore_dst(btot_in,bx_in, by_in,bz_in,v_in,time_in):

 #this makes from synthetic or observed solar wind the Dst index	
 #btot_in IMF total field, in nT
 #bx_in, by_in, bz_in are the components

 #v_in - the speed in km/s
 #vx_in - the solar wind speed x component (GSE or GSM?) in km/s
 #time_in - the time in matplotlib date format

 #these parameters are not yet derived from 3DCORE, so we set them reasonably:
 
 #vx_in is the the same as the v_in, vx needs to be positive
 vx_in=v_in
 
 #density_in is set to a generic value of 10 protons per ccm
 density_in=np.zeros(len(bz_in))+10
  

 #################################################### first 2 models
 #define variables
 Ey=np.zeros(len(bz_in))
 #dynamic pressure is constant
 pdyn1=np.zeros(len(bz_in))+1 #set pdyn to 10 
 dststar1=np.zeros(len(bz_in))
 dstcalc1=np.zeros(len(bz_in))
 dststar2=np.zeros(len(bz_in))
 dstcalc2=np.zeros(len(bz_in))
 
 #set all fields above 0 to 0 
 bz_in_negind=np.where(bz_in > 0)  
 
 #important: make a deepcopy because you manipulate the input variable
 bzneg=copy.deepcopy(bz_in)
 bzneg[bz_in_negind]=0

 #define interplanetary electric field 
 Ey=v_in*abs(bzneg)*1e-3; #now Ey is in mV/m
 
 
 #### model 1: Burton et al. 1975 
 Ec=0.5  
 a=3.6*1e-5
 b=0.2*100 #*100 wegen anderer dynamic pressure einheit in Burton
 c=20  
 d=-1.5/1000 
 for i in range(len(bz_in)-1):
  if Ey[i] > Ec:
   F=d*(Ey[i]-Ec) 
  else: F=0
  #Burton 1975 seite 4208: Dst=Dst0+bP^1/2-c   / und b und c positiv  
  #this is the ring current Dst
  deltat_sec=(time_in[i+1]-time_in[i])*86400 #timesyn is in days - convert to seconds

  dststar1[i+1]=(F-a*dststar1[i])*deltat_sec+dststar1[i];  #deltat must be in seconds
 
  #this is the Dst of ring current and magnetopause currents 
  dstcalc1[i+1]=dststar1[i+1]+b*np.sqrt(pdyn1[i+1])-c; 

 ######## model 2: OBrien and McPherron 2000 
 #constants
 Ec=0.49
 b=7.26  
 c=11  #nT
 for i in range(len(bz_in)-1):
  if Ey[i] > Ec:            #Ey in mV m
   Q=-4.4*(Ey[i]-Ec) 
  else: Q=0
  tau=2.4*np.exp(9.74/(4.69+Ey[i])) #tau in hours
  #OBrien abstract: Dst=Dst*+7.26P^1/2 -11
  #this is the ring current Dst
  #dststar2[i+1]=(Q-dststar2[i]/tau)+dststar2[i] #t is pro stunde, time intervall ist auch 1h
  deltat_hours=(time_in[i+1]-time_in[i])*24 #time_in is in days - convert to hours
  #deltat=dststar
  dststar2[i+1]=((Q-dststar2[i]/tau))*deltat_hours+dststar2[i] #t is pro stunde, time intervall ist auch 1h
  #this is the Dst of ring current and magnetopause currents 
  dstcalc2[i+1]=dststar2[i+1]+b*np.sqrt(pdyn1[i+1])-c; 
  
 
 
  dst_temerin_li_out=0

 
  #---------------- loop over
  
  

 return (dstcalc1,dstcalc2, dst_temerin_li_out)
   
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

########################################################################################
############################### main program ###########################################
########################################################################################


plt.close('all')

#check if directory of event exists
os.path.isdir(outputdirectory)
#if not make new directory
if os.path.isdir(outputdirectory) == False: os.mkdir(outputdirectory)
#also make directory for movie
if os.path.isdir(outputdirectory+'/movie') == False: os.mkdir(outputdirectory+'/movie')


#reads all lines as strings
lines = open(inputfilename).read().splitlines()

#choose if figures are saved (slow computation, 1) or not (fast computation, 0)
save_figures=int(lines[1])

#define the location of the insitu measuring spacecraft 
#availabe: Earth, MESSENGER
insitu_name=lines[4]


V0=float(lines[7]) #initial CME speed in km/s
R0=float(lines[9])*695508  #at initial distance in km
T0str=lines[11] #at initial time
T0=mdates.date2num(sunpy.time.parse_time(T0str)) #at initial time
gamma=float(lines[13])*1e-7 #DBM gamma
w=float(lines[15])  #background wind speed

#background in nT for each component
background_field_comp=float(lines[17])

#simulation start time
sim_start_str=lines[20]
sim_start=mdates.date2num(sunpy.time.parse_time(sim_start_str))

#how long the simulation is made starting with T0
simulation_duration_in_days=float(lines[22])


#direction of CME
longitude=float(lines[26]) 
latitude=float(lines[28])

#axis direction is a position angle: 0 deg solar north, 90 solar east, 180 south, 270 west to solar equatorial plane
axisangle=float(lines[30])
inclination=float(lines[30])-90

#+1 right handed, -1 left handed
handedness=int(lines[32])

if handedness == -1: Hstr='L'
if handedness == 1: Hstr='R'

#in the uniform twist model, tau is the number of turns a field line makes per AU around the axis
tau=float(lines[34])

#these should be made better with expansion laws and general B laws + flux, helicity..
#radial diameter of torus central cross section at 1 AU
D1AU=float(lines[36])
#magnetic field 
b1AU=float(lines[38])


#grid and time resolution:
dpsi=float(lines[42]) #steps in degree for psi
dphi=float(lines[44]) # steps in degree for phi
dr=float(lines[46]) #steps in radial direction for cross section in AU
dt=float(lines[48])   #time resolution in hours
#this parameter says how far the in situ location can be from the next point
#so that a measurement is initiated 
measurement_parameter=float(lines[50])	#AU


##########################################

#This is the aspect ratio of the torus-cross section
#if ellipses want to be used, the magnetic field configuration needs to 
#be included for a 2.5D uniform twist field too, currently a 
#circular cross-sections is used
aspectratio=1   

number_of_timesteps=int(24/dt*simulation_duration_in_days)

#make simulation time array
sim_time=sim_start+np.arange(number_of_timesteps)/24.*dt


####################### derived parameters from relationships:

#make new file

#open file for logging results
logfile=outputdirectory+'/results_log_lon%03d_lat%03d_inc%03d_%s_%s.txt'  %(longitude, latitude, axisangle, Hstr, insitu_name)

log=open(logfile,'w')


######################################

print('') 
print('3DCORE runs with the following parameters:')
print('') 
print( 'input:',inputfilename)
print('') 
print( 'in situ spacecraft: ', insitu_name)
print('') 
print( 'CME:')
print( ' V0, R0, T0:', V0, ' ', R0/695508,'   ', T0str)
print( ' gamma and w:', gamma/1e-7, w)
print( ' background field:', background_field_comp)
print('') 
print( 'MFR:')
print( ' longitude:', longitude)
print( ' latitude:', latitude)
print( ' Axis position angle', axisangle)
print( ' handedness:', handedness)
print( ' tau:', tau)
print( ' D1AU:', D1AU)
print( ' B1AU:', b1AU)
print('') 
print( 'simulation:')
print( ' dpsi, dph, dr, dt, m:', dpsi, ' ', dphi, ' ', dr, ' ', dt, ' ', measurement_parameter)
print( ' runs from ',mdates.num2date(sim_time[0]), ' to ', mdates.num2date(sim_time[number_of_timesteps-1]))
print( ' duration in days: ', simulation_duration_in_days)
print( ' time resolution in hours: ', dt)
print( ' number of time steps: ', number_of_timesteps)

#same to log file
log.write('')
log.write( '3DCORE runs with the following parameters:')
log.write('')
log.write( 'input:'+inputfilename)
log.write('')
log.write('in situ spacecraft: '+ insitu_name)
log.write('')
log.write('CME:')
log.write(' V0, R0, T0:'+ str(V0)+ ' '+ str(R0/695508)+'   '+ T0str)
log.write(' gamma and w:'+ str(gamma/1e-7)+ str(w))
log.write(' background field:'+ str(background_field_comp))
log.write('')
log.write('MFR:')
log.write(' longitude:'+ str(longitude))
log.write(' latitude:'+ str(latitude))
log.write(' Axis position angle:'+ str(axisangle))
log.write(' handedness:'+ str(handedness))
log.write(' tau:'+ str(tau))
log.write(' D1AU:'+ str(D1AU))
log.write(' B1AU:'+ str(b1AU))
log.write('')
log.write('simulation:')
log.write(' dpsi, dph, dr, dt, m:'+ str(dpsi)+ ' '+ str(dphi)+ ' '+ str(dr)+ ' '+ str(dt)+ ' '+ str(measurement_parameter))
log.write(' runs from '+str(mdates.num2date(sim_time[0]))+ ' to '+ str(mdates.num2date(sim_time[number_of_timesteps-1])))
log.write(' duration in days: '+ str(simulation_duration_in_days))
log.write(' time resolution in hours: '+ str(dt))
log.write(' number of time steps: '+ str(number_of_timesteps))





##################################### initialize everything

AU=149597870.700 #km

#get spacecraft and planet positions
pos=getpositions('DATACAT/positions_2007_2018_HEEQ_6hours.sav')
#convert time to matplotlib format
#available as pos.mercury etc.
pos_time_num=time_to_num_cat(pos.time)


#closes window
plt.close(1)


#ARRAYS for SYNTHETIC IN SITU DATA 
#synthetic magnetic field components at in situ location
#these have the maximum size of the number timesteps in the simulation
bxi=(np.random.random(number_of_timesteps)*background_field_comp-background_field_comp/2)*2
byi=(np.random.random(number_of_timesteps)*background_field_comp-background_field_comp/2)*2
bzi=(np.random.random(number_of_timesteps)*background_field_comp-background_field_comp/2)*2
bti=np.sqrt(bxi**2+byi**2+bzi**2)

#these are for the measurements in solar equatorial coordinates
bxiseq=np.copy(bxi)
byiseq=np.copy(byi)
bziseq=np.copy(bzi)
btiseq=np.copy(bti)

#these are for the measurements in RTN coordinates
bxirtn=np.copy(bxi)
byirtn=np.copy(byi)
bzirtn=np.copy(bzi)
btirtn=np.copy(bti)


#these are for the measurements in GSE and below GSM

bxigse=np.copy(bxi)
byigse=np.copy(byi)
bzigse=np.copy(bzi)
btigse=np.copy(bti)

bxigsm=np.copy(bxi)
byigsm=np.copy(byi)
bzigsm=np.copy(bzi)
btigsm=np.copy(bti)



#speed is set to background wind
speedi=np.random.normal(loc=w,scale=10.0,size=number_of_timesteps)
#this array will then contain the times only of the MFR measurements
time_mfr=np.zeros(number_of_timesteps)

#this variable counts the frames and always starts with 0
frame_counter=0




################################ START LOOP FOR EACH FRAME


for p in range(0,number_of_timesteps,1):

 
 #clear windows
 plt.clf()
 
 print 
 #print( 'frame: ', p)
 complete=float(p)/number_of_timesteps*100
 print( 'percentage completed: %2.1f' %complete)
 print( 'time: ', mdates.num2date(sim_time[p]))
 
 sim_time_from_t0_sec=(sim_time[p]-T0)*86400
 print( 'time since T0 in seconds:', int(sim_time_from_t0_sec))
 
 
 ############## get parameters of torus for timestep:

 #sign for deceleration regime is +1
 sign=1 
 #for acceleration (background wind faster than CME) its -1
 if V0 < w: sign=-1
  
 #DBM equation for distance Vrsnak et al. 2013
 rapex=(sign/gamma*math.log1p(1+sign*gamma*(V0-w)*sim_time_from_t0_sec)+w*sim_time_from_t0_sec+R0)/AU
 print( 'Rapex(t) in AU: %.3f' %rapex)
 
 Vapex=(V0-w)/(1+sign*gamma*(V0-w)*sim_time_from_t0_sec)+w
 print( 'Vapex(t) in km/s: ', int(Vapex))
  
 rho1=(D1AU*((rapex)**1.14))/2  # Leitner 2007 1.14
 print( 'rho1(t) in AU: %.3f' %rho1)
 
 rho0=(rapex-rho1)/2
 print( 'rho0(t) in AU: %.3f' %rho0)
 
 b0=b1AU*(2*rho0)**(-1.64)## Leitner 2007
 print( 'B1AU: %.3f' %b1AU)
 print( 'B0(t): %.3f' %b0)


 ########################## define coordinates of the MC 
 psi=np.radians(np.r_[0:360.:dpsi]); #convert array of angles to radians
 phi=np.radians(np.r_[0:360.:dphi]); #convert array of angles to radians
 rcoord=np.r_[0.01:rho1:dr]     # in AU
 
 print(' ')
 print('----------- simulation step number ', p)

 #make cartesian arrays
 sizeof_cartesian_array=np.size(phi)*np.size(psi)#*size(eta)
 #print 'number of grid points', sizeof_cartesian_array
 x=np.zeros(sizeof_cartesian_array)
 y=np.zeros(sizeof_cartesian_array)
 z=np.zeros(sizeof_cartesian_array)


 

 sizeof_volume_array=np.size(phi)*np.size(psi)*np.size(rcoord)

 #print 'number of grid points', sizeof_cartesian_array
 xv=np.zeros(sizeof_volume_array)
 yv=np.zeros(sizeof_volume_array)
 zv=np.zeros(sizeof_volume_array)
 bx=np.zeros(sizeof_volume_array)
 by=np.zeros(sizeof_volume_array)
 bz=np.zeros(sizeof_volume_array)
 speed=np.zeros(sizeof_volume_array)


 ########################### #calculate torus surface
 #go through all points of the torus surface, this is a surface with rho1 as minor radius
 i=0
 #psiall=np.zeros(np.size(phi)*np.size(psi))
 #phiall=np.zeros(np.size(phi)*np.size(psi))

 for phiind in range(0,np.size(phi)):
 	for psiind in range(0,np.size(psi)): 	
							x[i],y[i],z[i]=eltor2cart(phi[phiind],psi[psiind],rho1, aspectratio)
							#phiall[i]=phi[phiind]
							#psiall[i]=psi[psiind]
							i=i+1
							
	######## calculate all points in torus volume
 ##### includes now all r values = array rcoord
 i=0
 for psiind in range(0,np.size(psi)):
  for rind in range(0,np.size(rcoord)):
   for phiind in range(0,np.size(phi)):
			     xv[i],yv[i],zv[i]=eltor2cart(phi[phiind],psi[psiind],rcoord[rind],aspectratio)
			     # make cloud vector field -> each point of xv,yv,zv gets the lokal uniform twist components 
			     # conversion of cloud vector field to cartesian coordinates for a non-rotated MFR
			     bx[i],by[i],bz[i]=unifield2cart(xv[i],yv[i],zv[i],rcoord[rind],b0,tau,phi[phiind],psi[psiind],handedness)
			     i=i+1
						 	
 #shift center of torus so that the torus ends on the Sun
 x=x+rho0
 xv=xv+rho0
		 	
	
 ############### plot everything 
 sns.set_style("white") 
 fig2=plt.figure(1,figsize=(20,20),dpi=300)
 ax = fig2.add_subplot(1,1,1, projection='3d')
 fig2.set_size_inches(10, 10, forward=True)
 

 #trick for setting the aspect ratio correct in 3D plots 
 ax.set_aspect('equal')
 MAX = 0.9
 for direction in (-1, 1):
  for point in np.diag(direction * MAX * np.array([1,1,1])):
   ax.plot([point[0]], [point[1]], [point[2]], 'w')

 # plot solar equatorial plane
 xe = np.linspace(-1.5, 1.5, 4)
 ye = np.linspace(-1.5, 1.5, 4)
 zecl=np.zeros([np.size(xe),np.size(ye)])
 xecl,yecl = np.meshgrid(xe, ye)
 ax.plot_wireframe(xecl,yecl,zecl,color='g', linestyle='dashed',lw=0.5, zorder=0.5)

 #view
 #ax.view_init(25,-50+p/10.)
 ax.view_init(42,33+p/10)
 
 
 #make this larger for viewpoint thats further away
 ax.dist = 6
 
 
 
 

 ############## plot spacecraft positions in HEEQ (time dependent)


 #check which index is closest in positions to current time
 #pos_time_num vs. sim_time 
 sim_timend=np.argmin(abs(sim_time[p]-pos_time_num))
 
 sunsize=15
 planetsize=8
 spacecraftsize=4


 #Sun position
 Sun=[0,0,0]
 #ax.plot(Sun[0],Sun[1],Sun[2])#,marker='o',markersize=10	)
 ax.plot([Sun[0],Sun[0]],[Sun[1],Sun[1]],[Sun[2],Sun[2]],color='y',markersize=sunsize,marker='o')

 
 
 #Earth_L1 position
 earth_L1=sphere2cart(pos.earth_L1[0,sim_timend], pos.earth_L1[1,sim_timend], pos.earth_L1[2,sim_timend])
 earth_L1x=float(earth_L1[0])
 earth_L1y=float(earth_L1[1])
 earth_L1z=float(earth_L1[2])
 ax.plot([earth_L1x, earth_L1x],[earth_L1y,earth_L1y],[earth_L1z,earth_L1z],color='mediumseagreen',markersize=planetsize,marker='o',zorder=0.5)

 #Sun-Earth line
 Earth=[earth_L1x,earth_L1y,earth_L1z]
 ax.plot([Sun[0],earth_L1x],[Sun[1],earth_L1y],[Sun[2],earth_L1z],color='k',lw=0.3, zorder=0.5)

 
  
 #Venus position
 Venus=sphere2cart(pos.venus[0,sim_timend], pos.venus[1,sim_timend], pos.venus[2,sim_timend])
 venusx=float(Venus[0])
 venusy=float(Venus[1])
 venusz=float(Venus[2])
 Venus=[venusx,venusy,venusz]
 ax.plot([venusx, venusx],[venusy,venusy],[venusz,venusz],c='orange',markersize= planetsize,marker='o')

 #MESSENGER position
 mes=sphere2cart(pos.messenger[0,sim_timend], pos.messenger[1,sim_timend], pos.messenger[2,sim_timend])
 mesx=float(mes[0])
 mesy=float(mes[1])
 mesz=float(mes[2])
 MESSENGER=[mesx,mesy,mesz]
 ax.plot([mesx, mesx],[mesy,mesy],[mesz,mesz],color='dimgrey',markersize=spacecraftsize,marker='o',zorder=1)

 #Mercury position
 mercury=sphere2cart(pos.mercury[0,sim_timend], pos.mercury[1,sim_timend], pos.mercury[2,sim_timend])
 mercuryx=float(mercury[0])
 mercuryy=float(mercury[1])
 mercuryz=float(mercury[2])
 ax.plot([mercuryx, mercuryx],[mercuryy,mercuryy],[mercuryz,mercuryz],color='dimgrey',markersize= planetsize,marker='o',zorder=1)

 #STEREO-A position
 sta=sphere2cart(pos.sta[0,sim_timend], pos.sta[1,sim_timend], pos.sta[2,sim_timend])
 stax=float(sta[0])
 stay=float(sta[1])
 staz=float(sta[2])
 ax.plot([stax, stax],[stay,stay],[staz,staz],color='red',markersize=spacecraftsize,marker='o')

 #STEREO-B position
 stb=sphere2cart(pos.stb[0,sim_timend], pos.stb[1,sim_timend], pos.stb[2,sim_timend])
 stbx=float(stb[0])
 stby=float(stb[1])
 stbz=float(stb[2])
 ax.plot([stbx, stbx],[stby,stby],[stbz,stbz],color='royalblue',markersize=spacecraftsize,marker='o')
 

 xcirc=np.zeros(np.size(np.arange(0,2*np.pi+0.5,0.1)))
 ycirc=np.zeros(np.size(np.arange(0,2*np.pi+0.5,0.1)))
 xmcirc=np.zeros(np.size(np.arange(0,2*np.pi+0.5,0.1)))
 ymcirc=np.zeros(np.size(np.arange(0,2*np.pi+0.5,0.1)))
 zcirc=np.zeros(np.size(np.arange(0,2*np.pi+0.5,0.1)))
 
 qi=0
 #1 AU circle
 for q in np.arange(0, 2* np.pi+0.5, 0.1):
    xcirc[qi]=np.linalg.norm(earth_L1)*np.sin(q)
    ycirc[qi]=np.linalg.norm(earth_L1)*np.cos(q)
    xmcirc[qi]=np.linalg.norm(mes)*np.sin(q)
    ymcirc[qi]=np.linalg.norm(mes)*np.cos(q)
    zcirc[qi]=0          
    qi=qi+1
    
 ax.plot(xcirc,ycirc,zcirc,color='mediumseagreen',lw=0.3, zorder=1)
 ax.plot(xmcirc,ymcirc,zcirc,color='dimgrey',lw=0.3, zorder=1)
 

 #these are further available
 #index 0 is rdist, index 1 is longitude, index 2 latitude

 #ax.scatter(pos.earth[1,sim_timend], pos.earth[0,sim_timend], s=50, c='mediumseagreen', alpha=1,lw=0)
 #ax.scatter(pos.mars[1,sim_timend], pos.mars[0,sim_timend], s=50, c='orangered', alpha=1,lw=0)
 #ax.scatter(pos.ulysses[1,sim_timend], pos.ulysses[0,sim_timend], s=25, c='darkolivegreen', alpha=1,lw=0,marker='s')
 #ax.scatter(pos.msl[1,sim_timend], pos.msl[0,sim_timend], s=25, c='magenta', alpha=1,lw=0,marker='s')
 #ax.scatter(pos.maven[1,sim_timend], pos.maven[0,sim_timend], s=25, c='steelblue', alpha=1,lw=0, marker='s')
 #ax.scatter(pos.rosetta[1,sim_timend], pos.rosetta[0,sim_timend], s=25, c='black', alpha=1,lw=0, marker='s')

 

 ######################################### Rotation of the cloud

 #by convention at 0 deg longitude and latitude the MFR points towards the Earth, 
 #with 0 inclination the axis points to solar north, rotates counterclockwise as a position angle

 #1 Direction longitude (rotation around z axis)
 #negative sign so that positive angles are to right or solar west

 #surface
 x,y,z=rotate(x,y,z,-longitude,0,0,1) 
 #volume
 xv,yv,zv=rotate(xv,yv,zv,-longitude,0,0,1) 

 #the axis of rotation changes now depending on the longitude!

 #2 Direction latitude (the y axis is shifted by the angle longitude)
 #surface
 x,y,z=rotate(x,y,z,latitude,-np.sin(np.radians(longitude)),np.cos(np.radians(longitude)),0)
 #volume
 xv,yv,zv=rotate(xv,yv,zv,latitude,-np.sin(np.radians(longitude)),np.cos(np.radians(longitude)),0)

 #latitude from  (south) to 90 (north), theta from 0 (north) to 180 (south)
 theta=abs(latitude-90)
 

 #3 Orientation  (the inclination must be rotated around the "x"-axis 
 # which is now shifted by both angles latitude and longitude according 
 # to the usual spherical coordinates transformation
 #surface
 x,y,z=rotate(x,y,z,-inclination,np.cos(np.radians(longitude))*np.sin(np.radians(theta)),np.sin(np.radians(longitude))*np.sin(np.radians(theta)),np.cos(np.radians(theta)))
 #volume
 xv,yv,zv=rotate(xv,yv,zv,-inclination,np.cos(np.radians(longitude))*np.sin(np.radians(theta)),np.sin(np.radians(longitude))*np.sin(np.radians(theta)),np.cos(np.radians(theta)))



 ###################################### Plot cloud


 #plot torus as wire works
 #light = LightSource(90, 20)
 #illuminated_surface = light.shade(z, cmap=cm.coolwarm)

 #**this does not work in the new matplotlib version 2.2.2
 #ax.plot_wireframe(x,y,z,color='coral',zorder=3, lw=3)#
 
 #plot grid points inside as scatter works
 ax.scatter(xv,yv,zv,color='coral')
  
 #does not work -need to make psi and phi for all x, y, z
 #tritor=matplotlib.tri.Triangulation(phiall,psiall)
 #ax.plot_trisurf(xv, yv, zv, triangles=tritor.triangles, cmap=plt.cm.Spectral)
 
 #ax.plot_surface(xm, ym, zm, rstride=1, cstride=1, linewidth=0, antialiased=False)#, facecolors=illuminated_surface)
 #just tries below
 #X = np.outer(np.sin(u), np.sin(v))
 
 #create an rgb array for single-color surfaces.
 #white = np.ones((Z.shape[0], Z.shape[1], 3))
 #red = white * np.array([1,0,0])
 #green = white * np.array([0,1,0])
 #blue = white * np.array([0,0,1])

                
 #light = LightSource(90, 20)
 #illuminated_surface = light.shade(Z, cmap=cm.coolwarm)

 #rgb = np.ones((Z.shape[0], Z.shape[1], 3))

 #illuminated_surface = light.shade_rgb(rgb * red, Z)
 #ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, antialiased=False, facecolors=illuminated_surface)

 #tritor=matplotlib.tri.Triangulation(x,z)

 #ax.plot_trisurf(x, y, z, triangles=tritor.triangles, cmap=plt.cm.Spectral)


 
 plt.figtext(0.2,0.95,'3DCORE  lat: '+str(int(latitude))+'  long: '+str(int(longitude))+'  inc: '+str(int(axisangle)) , fontsize=15, ha='center')
   
 plt.xlabel('x [AU]')
 plt.ylabel('y [AU]')
 
 #plot text for date extra so it does not move 
 frame_time_str=str(mdates.num2date(sim_time[p]))
 plt.figtext(0.5,0.95,frame_time_str[0:16], fontsize=15, ha='center')
 #plot main parameters of Apex and B0
 rapexstr='%.3f' %rapex
 plt.figtext(0.8,0.95,'B: '+ str(int(b0))+ ' nT      V: '+ str(int(Vapex))+' km/s     R: '+ str(rapexstr)+' AU' , fontsize=15, ha='center')  

 #year
 #plt.figtext(0.40,0.77,frame_time_str[0:4], fontsize=15, ha='center')
 #month
 #plt.figtext(0.48,0.77,frame_time_str[5:7], fontsize=15, ha='center')
 #day
 #plt.figtext(0.50,0.77,frame_time_str[8:10], fontsize=15, ha='center')
 #hours and minutes
 #plt.figtext(0.53,0.77,frame_time_str[11:13]+':'+frame_time_str[14:16], fontsize=15, ha='center')

  
 #ax.legend()
 plt.show()  
 #frame numbers for ffmpeg need to start at 00
 filename=outputdirectory+'/movie/B_poc_%04d.png' %frame_counter
 if save_figures == 1:
   plt.savefig(filename, dpi=150)
 frame_counter=frame_counter+1
 



 ########################################### Extract magnetic field
 
 #xv,yv,zv are positions of B vectors

 if insitu_name == 'Earth': insitu=Earth
 if insitu_name == 'Venus': insitu=Venus
 if insitu_name == 'MESSENGER': insitu=MESSENGER
 
 #get components of location
 xi=insitu[0]
 yi=insitu[1]
 zi=insitu[2]
 
 #array where distances from in situ location to each cloud point are saved
 dist2insitu=np.zeros(np.size(xv))
 
 for k in range(0,np.size(xv)):
  dist2insitu[k]=np.sqrt((xi-xv[k])**2+(yi-yv[k])**2+(zi-zv[k])**2)
 
 #get index in array of minimum distance
 mindist=np.argmin(dist2insitu) 
 
 #check if in situ is outside of MFR
 if dist2insitu[mindist] > measurement_parameter:
  print(insitu_name,'is outside of MFR')
 else: 
  #in case in situ spacecraft is inside, save magnetic field components of nearest point
  #bxi,byi,bzi are the observed magnetic field components at each timestep at the in situ location (already defined)
  #this is still a cloud cartesian system
  bxi[p]=bx[mindist]
  byi[p]=by[mindist]
  bzi[p]=bz[mindist]
  bti[p]=np.sqrt(bx[mindist]**2+by[mindist]**2+bz[mindist]**2)
  
  #now rotate to cartesian Sun-origin system 
  #1 Direction longitude (rotation around z axis)
  bxi[p],byi[p],bzi[p]=rotate(bxi[p],byi[p],bzi[p],-longitude,0,0,1) 
  #2 Direction latitude (the y axis is shifted by the angle longitude)
  bxi[p],byi[p],bzi[p]=rotate(bxi[p],byi[p],bzi[p],latitude,-np.sin(np.radians(longitude)),np.cos(np.radians(longitude)),0)
  #latitude from -90 (south) to 90 (north), theta from 0 (north) to 180 (south)
  theta=abs(latitude-90)
  #3 Inclination  (the inclination must be rotated around the "x"-axis)
  bxi[p],byi[p],bzi[p]=rotate(bxi[p],byi[p],bzi[p],-inclination,np.cos(np.radians(longitude))*np.sin(np.radians(theta)),np.sin(np.radians(longitude))*np.sin(np.radians(theta)),np.cos(np.radians(theta)))

  
  #magnetic field vector at this timestep
  Bi=[bxi[p],byi[p],bzi[p]]
  print('B origin x,y,z: %4.1f %4.1f %4.1f nT' %(Bi[0],Bi[1],Bi[2]))
 
  #calculate speed by assumption of self-similar expansion
  #first calculate heliocentric distance of measurement point
  rmindist=np.sqrt((xv[mindist])**2+(yv[mindist])**2+(zv[mindist])**2)
  #then use self-similar assumption as in Moestl and Davies 2013
  speedi[p]=rmindist*Vapex/rapex
  
  #****** ad hoc fix that the speed goes not too low in the flux rope
  #if speedi[p] < w: 
  #      speedi[p]=w
  #save time of mfr observation
  time_mfr[p]=sim_time[p]
  
  ########## conversion to coordinate system for comparison to observations
  
  #(1)
  #convert to SEQ system for in situ spacecraft position:   
  #Xseq from Sun to spacecraft, but in SEQ plane, this is a unit vector    
  Xseq=[insitu[0], insitu[1], 0]/np.linalg.norm([insitu[0], insitu[1],0])
  #Yseq is cross product with solar rotation axis at 0, 0, 1, which defines Zseq
  Zseq=[0,0,1]  
  Yseq=np.cross(Zseq, Xseq)/np.linalg.norm(np.cross(Zseq, Xseq))
  #Xseq, Yseq, and Zseq form a right handed system
  
  #project into new system
  bxiseq[p]=np.dot(Bi,Xseq)
  byiseq[p]=np.dot(Bi,Yseq)
  bziseq[p]=np.dot(Bi,Zseq)
  btiseq[p]=np.sqrt(bxiseq[p]**2+byiseq[p]**2+bziseq[p]**2)
  Biseq=[bxiseq[p],byiseq[p],bziseq[p]]
  
  print( 'B SEQ x,y,z: %4.1f %4.1f %4.1f nT' %(Biseq[0],Biseq[1],Biseq[2]))

 
  #(2)  
  #convert magnetic field components to RTN
  
  #R from Sun to spacecraft
  Xrtn=[insitu[0], insitu[1],insitu[2]]/np.linalg.norm([insitu[0], insitu[1],insitu[2]])
  #solar rotation axis at 0, 0, 1
  solrot=[0,0,1]
  Yrtn=np.cross(solrot,Xrtn)/np.linalg.norm(np.cross(solrot,Xrtn))
  Zrtn=np.cross(Xrtn, Yrtn)/np.linalg.norm(np.cross(Xrtn, Yrtn))


  #project into new system
  bxirtn[p]=np.dot(Bi,Xrtn)
  byirtn[p]=np.dot(Bi,Yrtn)
  bzirtn[p]=np.dot(Bi,Zrtn)
  btirtn[p]=np.sqrt(bxirtn[p]**2+byirtn[p]**2+bzirtn[p]**2)
  
  Birtn=[bxirtn[p],byirtn[p],bzirtn[p]]
  
  print( 'B RTN r,t,n: %4.1f %4.1f %4.1f nT' %(Birtn[0],Birtn[1],Birtn[2]))




  #for Earth, additionally make GSM and GSE
  if insitu_name == 'Earth': 
   #(3)  
   #convert magnetic field components in the origin system HEEQ (Bi) to HAE, then HAE to HEE, and then HEE to GSE  
   #Hapgood 1992, use a conversion matrix S2^-1 x S1 obtained from sections 5.1 and 5.2
  

   #get MJD in Hapgood 1992
   sim_time_sunpy=sunpy.time.break_time(mdates.num2date(sim_time[p]))
   jd=sunpy.time.julian_day(sim_time_sunpy)
   mjd=float(int(jd-2400000.5)) #use modified julian date    
   
   #then lambda_sun
   T00=(mjd-51544.5)/36525.0
   dobj=mdates.num2date(sim_time[p])
   UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours   
   LAMBDA=280.460+36000.772*T00+0.04107*UT
   M=357.528+35999.050*T00+0.04107*UT
   #lt2 is lambdasun in Hapgood, equation 5, here in rad
   lt2=(LAMBDA+(1.915-0.0048*T00)*np.sin(M*np.pi/180)+0.020*np.sin(2*M*np.pi/180))*np.pi/180
   #note that some of these equations are repeated later for the GSE to GSM conversion
   
   ################### lambdasun und z, S1 matrix for HAE to HEE
   S1=np.matrix([[np.cos(lt2+np.pi), np.sin(lt2+np.pi),  0], [-np.sin(lt2+np.pi) , np.cos(lt2+np.pi) , 0], [0,  0,  1]])
   
   #create S2 matrix with angles with reversed sign for transformation HEEQ to HAE
   omega_node=(73.6667+0.013958*((mjd+3242)/365.25))*np.pi/180 #in rad
   S2_omega=np.matrix([[np.cos(-omega_node), np.sin(-omega_node),  0], [-np.sin(-omega_node) , np.cos(-omega_node) , 0], [0,  0,  1]])
   inclination_ecl=7.25*np.pi/180
   S2_incl=np.matrix([[1,0,0],[0,np.cos(-inclination_ecl), np.sin(-inclination_ecl)], [0, -np.sin(-inclination_ecl), np.cos(-inclination_ecl)]])
   
   #calculate theta
   theta_node=np.arctan(np.cos(inclination_ecl)*np.tan(lt2-omega_node)) 
   
   #quadrant of theta must be opposite lt2 - omega_node Hapgood 1992 end of section 5   
   #get lambda-omega angle in degree mod 360   
   lambda_omega_deg=np.mod(lt2-omega_node,2*np.pi)*180/np.pi
   #get theta_node in deg
   theta_node_deg=theta_node*180/np.pi

   #if in same quadrant, then theta_node = theta_node +pi   
   if abs(lambda_omega_deg-theta_node_deg) < 180: theta_node=theta_node+np.pi

   S2_theta=np.matrix([[np.cos(-theta_node), np.sin(-theta_node),  0], [-np.sin(-theta_node) , np.cos(-theta_node) , 0], [0,  0,  1]])
   
   #make S2 matrix
   S2=np.dot(np.dot(S2_omega,S2_incl),S2_theta)

   #this is the matrix S2^-1 x S1
   HEEQ_to_HEE_matrix=np.dot(S1, S2)
  
   #convert HEEQ components to HEE
   
   HEEQ=np.matrix([[Bi[0]],[Bi[1]],[Bi[2]]]) 
   HEE=np.dot(HEEQ_to_HEE_matrix,HEEQ)
     
   #change of sign HEE X / Y to GSE is needed
   bxigse[p]=-HEE.item(0)
   byigse[p]=-HEE.item(1)
   bzigse[p]=HEE.item(2)
   btigse[p]=np.sqrt(bxigse[p]**2+byigse[p]**2+bzigse[p]**2)

   Bigse=[bxigse[p],byigse[p],bzigse[p]]
   print( 'B GSE x,y,z: %4.1f %4.1f %4.1f nT' %(Bigse[0],Bigse[1],Bigse[2]))




   #(4) #for Earth also calculate GSM coordinates, for Dst prediction
   #need to get angle psigsm after Hapgood 1992/1997, section 4.3
   #define position of geomagnetic pole in GEO coordinates
   #convert to sunpy time
   
   sim_time_sunpy=sunpy.time.break_time(mdates.num2date(sim_time[p]))
   jd=sunpy.time.julian_day(sim_time_sunpy)
   mjd=float(int(jd-2400000.5)) #use modified julian date    
   
   pgeo=78.8+4.283*((mjd-46066)/365.25)*0.01 #in degrees
   lgeo=289.1-1.413*((mjd-46066)/365.25)*0.01 #in degrees
   #GEO vector
   Qg=[np.cos(pgeo*np.pi/180)*np.cos(lgeo*np.pi/180), np.cos(pgeo*np.pi/180)*np.sin(lgeo*np.pi/180), np.sin(pgeo*np.pi/180)]
   #now move to equation at the end of the section, which goes back to equations 2 and 4:
   #CREATE T1
   T00=(mjd-51544.5)/36525.0
   dobj=mdates.num2date(sim_time[p])
   UT=dobj.hour + dobj.minute / 60. + dobj.second / 3600. #time in UT in hours   
   zeta=(100.461+36000.770*T00+15.04107*UT)*np.pi/180
   ################### theta und z
   T1=np.matrix([[np.cos(zeta), np.sin(zeta),  0], [-np.sin(zeta) , np.cos(zeta) , 0], [0,  0,  1]]) #angle for transpose
   #CREATE T2
   LAMBDA=280.460+36000.772*T00+0.04107*UT
   M=357.528+35999.050*T00+0.04107*UT
   lt2=(LAMBDA+(1.915-0.0048*T00)*np.sin(M*np.pi/180)+0.020*np.sin(2*M*np.pi/180))*np.pi/180
   ##################### lamdbda und Z
   t2z=np.matrix([[np.cos(lt2), np.sin(lt2),  0], [-np.sin(lt2) , np.cos(lt2) , 0], [0,  0,  1]])
   et2=(23.439-0.013*T00)*np.pi/180
   ###################### epsilon und x
   t2x=np.matrix([[1,0,0],[0,np.cos(et2), np.sin(et2)], [0, -np.sin(et2), np.cos(et2)]])
   T2=np.dot(t2z,t2x)  #equation 4 in Hapgood 1992
   #matrix multiplications   
   T2T1t=np.dot(T2,np.matrix.transpose(T1))
   ################
   Qe=np.dot(T2T1t,Qg) #Q=T2*T1^-1*Qq
   psigsm=np.arctan(Qe.item(1)/Qe.item(2)) #arctan(ye/ze) in between -pi/2 to +pi/2
   
   T3=np.matrix([[1,0,0],[0,np.cos(-psigsm), np.sin(-psigsm)], [0, -np.sin(-psigsm), np.cos(-psigsm)]])
   GSE=np.matrix([[bxigse[p]],[byigse[p]],[bzigse[p]]]) 
   GSM=np.dot(T3,GSE)   #equation 6 in Hapgood
   bxigsm[p]=GSM.item(0)
   byigsm[p]=GSM.item(1)
   bzigsm[p]=GSM.item(2)
  	##########################################################    
   
   #other version: rotate B in GSE around XGSE by angle psigsm - needs to be in degrees
   #bxigsm[p],byigsm[p],bzigsm[p]=rotate(bxigse[p],byigse[p],bzigse[p],psigsm,Xgse[0],Xgse[1],Xgse[2])

   Bigsm=[bxigsm[p],byigsm[p],bzigsm[p]]
   btigsm[p]=np.sqrt(bxigsm[p]**2+byigsm[p]**2+bzigsm[p]**2)

   print( 'B GSM x,y,z: %4.1f %4.1f %4.1f nT' %(Bigsm[0],Bigsm[1],Bigsm[2]))


  
 
########################################### end of for loop for MC evolution











#save figure of last frame for reference
filename=outputdirectory+'/3D_MFR_paper1_lon%03d_lat%03d_inc%03d_%s_%s.eps'  %(longitude, latitude, axisangle, Hstr, insitu_name)
plt.savefig(filename)
filename=outputdirectory+'/3D_MFR_paper1_lon%03d_lat%03d_inc%03d_%s_%s.png'  %(longitude, latitude, axisangle, Hstr, insitu_name)
plt.savefig(filename, dpi=300)


print()

print( '------------------------------------')


print( '3DCORE Run done.')

print()


#if MFR was observed, write out parameters as measured from observations
if np.size(time_mfr) > 0:
   #this contains only the indices when the MFR was observed
   time_mfr_index=np.where(time_mfr > 0)
   #trim zeros of observation time so that time_mfr only contains the observed mfr times
   time_mfr=np.trim_zeros(time_mfr, trim='fb')
   print( 'The MFR impacts ',insitu_name,' on ', mdates.num2date(time_mfr[0]))
   print( 'The MFR exits ',insitu_name,' on ', mdates.num2date(time_mfr[np.size(time_mfr)-1]))
   mfr_duration=(time_mfr[np.size(time_mfr)-1]-time_mfr[0])*24
   print( 'The MFR duration is', mfr_duration, ' hours')
   print()
   print( 'Max of total B during MFR GSM: %4.1f' %np.max(btirtn[time_mfr_index]))
   print( 'Mean of total B during MFR GSM: %4.1f' %np.mean(btirtn[time_mfr_index]))
   print( 'Std of total B during MFR GSM: %4.1f' %np.std(btirtn[time_mfr_index]))
   print( 'Min of Bz during MFR GSM: %4.1f' %np.min(bzirtn[time_mfr_index]))

   #same to logfile   
   log.write(' ')
   log.write('The MFR impacts '+insitu_name+' on '+ str(mdates.num2date(time_mfr[0])))
   log.write('The MFR exits '+insitu_name+' on '+ str(mdates.num2date(time_mfr[np.size(time_mfr)-1])))
   log.write('The MFR duration is'+ str(mfr_duration)+ ' hours')
   log.write(' ')
   log.write('Max of total B during MFR GSM:: %4.1f'+str(np.max(btigsm[time_mfr_index])))
   log.write('Mean of total B during MFR GSM:: %4.1f'+str(np.mean(btigsm[time_mfr_index])))
   log.write('Std of total B during MFR GSM:: %4.1f'+str(np.std(btigsm[time_mfr_index])))
   log.write('Min of Bz during MFR GSM:: %4.1f' +str(np.min(bzigsm[time_mfr_index])))
   
   
  


##################### make plot of synthetic magnetic field components and speed



#load observations for comparison for the July 13 2013 event, here Wind and MESSENGER are used
[datames,datawin,datamesn,datames_time,datawin_time,datamesn_time]=pickle.load( open( "DATACAT/mes_win_july2013_3dcore_paper1.p", "rb" ) )

#load OMNI observed B in GSE and GSM and Dst, here Dst is used; OMNI is NOT used for the data plot, which is Wind
[bxe,bye,bze,bxm,bym,bzm,dst,speed, times1]= pickle.load( open( "DATACAT/omni2_gse_gsm_dst_small.p", "rb" ) )


#for comparison to real data: extract data for simulation start time +10 days
ndays=10
s=mdates.date2num(sunpy.time.parse_time(mdates.num2date(sim_time[0])))
#take this starting date + 10 days for a smaller array for calculating Dst
#"i" stands for interval
ind=np.where(s==times1)[0][0]

bxmi=bxm[ind:ind+24*ndays]
bymi=bym[ind:ind+24*ndays]
bzmi=bzm[ind:ind+24*ndays]
vi=speed[ind:ind+24*ndays]
timesi=times1[ind:ind+24*ndays]
dsti=dst[ind:ind+24*ndays]
btmi=(bxmi**2+bymi**2+bzmi**2)**0.5  



print( 'plot observed and synthetic magnetic field components B and speed V')

#plot controls from seaborn
sns.set_context("talk")     
sns.set_style("darkgrid")  
sns.set_style("ticks")  



#make plot depending on target

#plot from start of MFR observations -1 day to + 2 days after end
plotstartdate=mdates.num2date(time_mfr[0]-1)
plotenddate=mdates.num2date(time_mfr[np.size(time_mfr)-1]+2)

myformat = mdates.DateFormatter('%d %h %Hh')


###############load ICMECAT for ICME times and parameters
filename_icmecat='ICMECAT/HELCATS_ICMECAT_v10_SCEQ.sav'
i=getcat(filename_icmecat)
sc=i.icmecat['sc_insitu'] #string
sc=decode_array(sc)

#get indices of events in different spacecraft
vexind=np.where(sc == 'VEX')
staind=np.where(sc == 'STEREO-A')
stbind=np.where(sc == 'STEREO-B')
winind=np.where(sc == 'Wind')
mesind=np.where(sc == 'MESSENGER')
ulyind=np.where(sc == 'ULYSSES')

#get boundary times
icme_start_time_str=i.icmecat['icme_start_time']
mo_start_time_str=i.icmecat['mo_start_time']
mo_end_time_str=i.icmecat['mo_end_time']

#save them as matplotlib date number
icme_start_time_num=time_to_num_cat(icme_start_time_str)
mo_start_time_num=time_to_num_cat(mo_start_time_str)
mo_end_time_num=time_to_num_cat(mo_end_time_str)


###################

if insitu_name == 'Earth':


  plt.figure(2, figsize=(10,10), dpi=100, facecolor='w', edgecolor='w')

  ######### 1 magnetic field plot
  ax1 = plt.subplot2grid((3,1), (0, 0))

  #set mathtext font so it looks the same as normal text 
  #matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
  plt.title('MFR magnetic field components at '+insitu_name)

  #obs
  ax1.plot_date(datawin_time, datawin.bx,'r',lw=1)
  ax1.plot_date(datawin_time, datawin.by,'g', lw=1)
  ax1.plot_date(datawin_time, datawin.bz,'b', lw=1)
  ax1.plot_date(datawin_time, datawin.btot,'k', label='observation',lw=1)
  
  #sim
  ax1.plot_date(sim_time, bxiseq,'r--')
  ax1.plot_date(sim_time, byiseq,'g--')
  ax1.plot_date(sim_time, bziseq,'b--')
  ax1.plot_date(sim_time, btiseq,'k--', label='simulation')
 
  plt.ylabel('magnetic field B (HEEQ) [nT]')
 
  #plt.ylabel('magnetic field B [nT]')
  plt.ylim(-max(btirtn)-20, max(btirtn)+20)
  plt.legend(loc=1,fontsize=10)
  plt.xlim((plotstartdate, plotenddate))
  ax1.xaxis.set_major_formatter(myformat)
  #ax1.grid(True)
  ax1.xaxis.set_major_locator(mdates.DayLocator())
  plt.tick_params( axis='x', labelbottom='off')
  
  #overplot MO boundaries at Wind
  ax1.plot_date([mo_start_time_num[winind],mo_start_time_num[winind]], [-max(btirtn)-20, max(btirtn)+20],'-k', lw=1, alpha=0.8)
  ax1.plot_date([mo_end_time_num[winind],mo_end_time_num[winind]], [-max(btirtn)-20, max(btirtn)+20],'-k', lw=1, alpha=0.8)

  
  

  ##### 2 speed plot
  ax2 = plt.subplot2grid((3,1), (1, 0))
  plt.title(r'Bulk plasma speed at '+insitu_name)
  
  #delete spikes for speed plot July 2013 event

  #spikes on 2013 12 July 9:47, 9:48 ; 13 July 1:21, 1:22, 1:23, 15 July 9:47, 9:48
  
  spike1=np.where(mdates.date2num(sunpy.time.parse_time('2013-07-12T09:47:00')) == datawin_time)
  spike2=np.where(mdates.date2num(sunpy.time.parse_time('2013-07-12T09:48:00')) == datawin_time)
  spike3=np.where(mdates.date2num(sunpy.time.parse_time('2013-07-13T01:21:00')) == datawin_time)
  spike4=np.where(mdates.date2num(sunpy.time.parse_time('2013-07-13T01:22:00')) == datawin_time)
  spike5=np.where(mdates.date2num(sunpy.time.parse_time('2013-07-13T01:23:00')) == datawin_time)
  spike6=np.where(mdates.date2num(sunpy.time.parse_time('2013-07-15T09:47:00')) == datawin_time)
  spike7=np.where(mdates.date2num(sunpy.time.parse_time('2013-07-15T09:48:00')) == datawin_time)
  spike8=np.where(mdates.date2num(sunpy.time.parse_time('2013-07-15T13:36:00')) == datawin_time)
  spike9=np.where(mdates.date2num(sunpy.time.parse_time('2013-07-15T13:37:00')) == datawin_time)
  spike10=np.where(mdates.date2num(sunpy.time.parse_time('2013-07-15T17:37:00')) == datawin_time)
  spike11=np.where(mdates.date2num(sunpy.time.parse_time('2013-07-15T17:38:00')) == datawin_time)


  datawin.vtot[spike1]=np.nan
  datawin.vtot[spike2]=np.nan
  datawin.vtot[spike3]=np.nan
  datawin.vtot[spike4]=np.nan
  datawin.vtot[spike5]=np.nan
  datawin.vtot[spike6]=np.nan
  datawin.vtot[spike7]=np.nan
  datawin.vtot[spike8]=np.nan
  datawin.vtot[spike9]=np.nan
  datawin.vtot[spike10]=np.nan
  datawin.vtot[spike11]=np.nan


  
  #observations
  ax2.plot_date(datawin_time,datawin.vtot, 'k', lw=1, alpha=0.9, label='observation')
  
  #simulation 
  ax2.plot_date(sim_time,speedi,'b--', lw=2, label='simulation')
  
  #overplot MO boundaries 
  ax2.plot_date([mo_start_time_num[winind],mo_start_time_num[winind]], [200, max(speedi)+200],'-k', lw=1, alpha=0.8)
  ax2.plot_date([mo_end_time_num[winind],mo_end_time_num[winind]], [200, max(speedi)+200],'-k', lw=1, alpha=0.8)


 
   
  plt.ylim((200, max(speedi)+200))
  plt.ylabel('plasma speed V [$\mathregular{km} \mathregular{\; s^{-1}}}$]')
  #plt.xlabel(r'Year '+frame_time_str[0:4])
  plt.xlim((plotstartdate, plotenddate))
  plt.legend(loc=1,fontsize=10)

  ax2.xaxis.set_major_locator(mdates.DayLocator())
  ax2.xaxis.set_major_formatter(myformat)

  plt.tick_params( axis='x', labelbottom='off')






  ######### 3  synthethic Dst plot only for Earth
  ax3 = plt.subplot2grid((3,1), (2, 0))

  #Observed Dst
  ax3.plot_date(times1, dst,'ko',label='observed Dst', markersize=5)


  #make OMNI data total field
  btm=np.sqrt(bxm**2+bym**2+bzm**2)
  #cut time from OMNI data  
  interval_omni_ind=np.where(np.logical_and(times1 > mdates.date2num(sunpy.time.parse_time('2013-07-05T00:00:00')),times1 < mdates.date2num(sunpy.time.parse_time('2013-07-20T00:00:00'))))
  
  [obs_dst_burton, obs_dst_obrien,obs_dst_temerin_li]=make_3dcore_dst(btm[interval_omni_ind],bxm[interval_omni_ind], bym[interval_omni_ind], bzm[interval_omni_ind],speed[interval_omni_ind], times1)
  #from obs OMNI
  ax3.plot_date(times1[interval_omni_ind], obs_dst_burton,'b-',label='obs. Burton et al. ')
  ax3.plot_date(times1[interval_omni_ind], obs_dst_obrien,'r-',label='obs. OBrien & McPherron')

  #with synthetic data
  [sim_dst_burton, sim_dst_obrien,sim_dst_temerin_li]=make_3dcore_dst(btigsm,bxigsm, byigsm, bzigsm,speedi,sim_time)
  #from sim
  ax3.plot_date(sim_time, sim_dst_burton,'b--',label='sim. Burton et al. ')
  ax3.plot_date(sim_time, sim_dst_obrien,'r--',label='sim. OBrien & McPherron')


  plt.title('Dst index from L1 observation / 3DCORE simulation')



  #plt.ylim((min(dst)-30, max(dst)+20))
  plt.ylim((-130, 50))

  plt.ylabel('Dst index [nT]')
  plt.legend(loc=3,fontsize=10)
  plt.xlim((plotstartdate, plotenddate))

  ax3.xaxis.set_major_locator(mdates.DayLocator())
  ax3.xaxis.set_major_formatter(myformat)
  #only for last panel
  plt.xlabel(r'Year '+frame_time_str[0:4])





#for MESSENGER no Dst, only field and speed
if insitu_name == 'MESSENGER':
 

  plt.figure(2, figsize=(10,8), dpi=100, facecolor='w', edgecolor='w')

  ######### 1 magnetic field plot
  ax1 = plt.subplot2grid((2,1), (0, 0))

  #set mathtext font so it looks the same as normal text 
  #matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
  plt.title('MFR magnetic field components at '+insitu_name)

  
  #observations
  ax1.plot_date(datames_time, datames.bx,'r', lw=1)
  ax1.plot_date(datames_time, datames.by,'g', lw=1)
  ax1.plot_date(datames_time, datames.bz,'b', lw=1)
  ax1.plot_date(datames_time, datames.btot,'k', label='observation', lw=1)

  #simulation
  ax1.plot_date(sim_time, bxiseq,'r--')
  ax1.plot_date(sim_time, byiseq,'g--')
  ax1.plot_date(sim_time, bziseq,'b--')
  ax1.plot_date(sim_time, btiseq,'k--',label='simulation')

  #overplot MO boundaries at MESSENGER
  #ax1.plot_date([icme_start_time_num[mesind],icme_start_time_num[mesind]], [-max(btirtn)-20, max(btirtn)+20],'-k', lw=0.5, alpha=0.8)
  ax1.plot_date([mo_start_time_num[mesind],mo_start_time_num[mesind]], [-max(btirtn)-20, max(btirtn)+20],'-k', lw=1, alpha=0.8)
  ax1.plot_date([mo_end_time_num[mesind],mo_end_time_num[mesind]], [-max(btirtn)-20, max(btirtn)+20],'-k', lw=1, alpha=0.8)

 
  plt.ylabel('magnetic field B (SCEQ) [nT]')

  #plt.ylabel('magnetic field B [nT]')
  plt.ylim(-max(btirtn)-20, max(btirtn)+20)
  plt.legend(loc=1,fontsize=10)
  plt.xlim((plotstartdate, plotenddate))
  ax1.xaxis.set_major_formatter(myformat)
  #ax1.grid(True)
  ax1.xaxis.set_major_locator(mdates.DayLocator())
  plt.tick_params( axis='x', labelbottom='off')
  

  ##### 2 speed plot
  ax2 = plt.subplot2grid((2,1), (1, 0))
  plt.title(r'Bulk plasma speed at '+insitu_name)
  
  #simulation 
  ax2.plot_date(sim_time,speedi,'b--', lw=2, label='simulation')
 
   
  plt.ylim((200, max(speedi)+200))
  plt.ylabel('plasma speed V [$\mathregular{km} \mathregular{\; s^{-1}}}$]')
  #plt.xlabel(r'Year '+frame_time_str[0:4])
  plt.xlim((plotstartdate, plotenddate))
  plt.legend(loc=1,fontsize=10)

  ax2.xaxis.set_major_locator(mdates.DayLocator())
  ax2.xaxis.set_major_formatter(myformat)

  plt.xlabel(r'Year '+frame_time_str[0:4])



################################


#save figure
filename=outputdirectory+'/syn_fields_paper1_lon%03d_lat%03d_inc%03d_%s_%s_dst.eps'  %(longitude, latitude, axisangle, Hstr, insitu_name)
plt.savefig(filename)
filename=outputdirectory+'/syn_fields_paper1_lon%03d_lat%03d_inc%03d_%s_%s_dst.png'  %(longitude, latitude, axisangle, Hstr, insitu_name)
plt.savefig(filename, dpi=300)


#make movie with ffmpeg if all figures are saved
if save_figures > 0:
  #first convert all convert to jpg
  os.system(os.getcwd()+'/ffmpeg -i '+outputdirectory+'/movie/B_poc_%04d.png '+outputdirectory+'/movie/B_poc_%04d.jpg')
  os.system(os.getcwd()+'/ffmpeg -r 20 -i '+outputdirectory+'/movie/B_poc_%04d.jpg -b 5000k -r 20 '+outputdirectory+'/3DCORE_movie_paper1.mp4 -y')


############ end of program



#This work is published under the MIT LICENSE
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
#PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
#FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
#TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
#OTHER DEALINGS IN THE SOFTWARE.


