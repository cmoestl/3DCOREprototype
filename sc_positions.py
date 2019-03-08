# Spacecraft trajectories for 3DCORE incl. Bepi Colombo, PSP, Solar Orbiter

#Author: C. Moestl, IWF Graz, Austria
#twitter @chrisoutofspace, https://github.com/cmoestl
#started December 2018

#needs python 3.6 with sunpy, heliopy, spiceypy, cdflib, seaborn

## MIT LICENSE
## Copyright 2018, Christian Moestl 
## Permission is hereby granted, free of charge, to any person obtaining a copy of this 
## software and associated documentation files (the "Software"), to deal in the Software
## without restriction, including without limitation the rights to use, copy, modify, 
## merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
## permit persons to whom the Software is furnished to do so, subject to the following 
## conditions:
## The above copyright notice and this permission notice shall be included in all copies 
## or substantial portions of the Software.
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
## PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
## HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
## CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
## OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



##########################################################################################
######################################### CODE START #####################################
##########################################################################################


import scipy.io
import os
import datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pdb
import urllib
import json
import pickle
import sunpy.time
import seaborn as sns
import sys
import heliopy.data.spice as spicedata
import heliopy.spice as spice
import astropy

from datetime import datetime, timedelta
#ignore warnings
#import warnings
#warnings.filterwarnings('ignore')
plt.close('all')



def sphere2cart(r, phi, theta):
    x = r*np.cos(theta)*np.cos(phi)
    y = r*np.cos(theta)*np.sin(phi)
    z = r*np.sin(theta)
    return (x, y, z) 
    
def cart2sphere(x,y,z):
    r = np.sqrt(x**2+ y**2 + z**2)            # r
    theta = np.arctan2(z,np.sqrt(x**2+ y**2))     # theta
    phi = np.arctan2(y,x)                        # phi
    return (r, theta, phi)

def getcat(filename):  
  print( 'reading positions in '+filename)
  pos=scipy.io.readsav(filename, verbose='true')  
  print( 'done reading positions')
  return pos 

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)
 
def time_to_num_cat(time_in):  

  #for time conversion from catalogue .sav to numerical time
  #this for 1-minute data or lower time resolution

  #for all catalogues
  #time_in is the time in format: 2007-11-17T07:20:00 or 2007-11-17T07:20Z
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
   time_str[j]
   #convert time to sunpy friendly time and to matplotlibdatetime
   #only for valid times so 9999 in year is not converted
   #pdb.set_trace()
   if year < 2100:
    	  time_num[j]=mdates.date2num(sunpy.time.parse_time(time_str[j]))
   j=j+1  
   #the date format in matplotlib is e.g. 735202.67569444
   #this is time in days since 0001-01-01 UTC, plus 1.
   
   #return time_num which is already an array and convert the list of strings to an array
  return time_num, np.array(time_str)




##################################################### MAIN ###############################


#get positions old style:
pos=getcat('DATACAT/positions_2007_2023_HEEQ_6hours.sav')
pos_time_num=time_to_num_cat(pos.time)[0]

#https://heliopy.readthedocs.io/en/stable/api/heliopy.spice.Trajectory.html#heliopy.spice.Trajectory.generate_positions

#load kernels if not here
spicedata.get_kernel('psp')
spicedata.get_kernel('planet_trajectories')

#Solar orbiter 'solo_2020'
#for PSP NAIF CODE is -96 (search for solar probe plus)
#-144        'SOLAR ORBITER'
#Earth 399






############ PSP

spice.furnish(spicedata.get_kernel('psp'))
psp=spice.Trajectory('-96')

starttime =datetime(2018, 10, 21)
endtime = datetime(2025, 8, 1)
psp_time = []
while starttime < endtime:
    psp_time.append(starttime)
    starttime += timedelta(days=1)
    
psp_time_num=mdates.date2num(psp_time)     
#frames https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html  Appendix. ``Built in'' Inertial Reference Frames
psp.generate_positions(psp_time,'Sun','ECLIPJ2000')
#pspe.generate_positions(time2,'Sun','ECLIPJ2000')
psp.change_units(astropy.units.AU)  
print('PSP')
print(psp.r)
print(psp.x)
print(psp.y)
print(psp.z)  


[psp_r, psp_theta, psp_phi]=cart2sphere(psp.x,psp.y,psp.z)


plt.figure(1, figsize=(12,9))
plt.plot_date(psp_time,psp_r,'-', label='R')

plt.plot_date(psp_time,psp_theta,'-',label='theta')

plt.plot_date(psp_time,psp_phi,'-',label='phi')


#plt.plot_date(psp_time,psp.x,'-',label='X')
#plt.plot_date(psp_time,psp.y,'-',label='Y')
#plt.plot_date(psp_time,psp.z,'-',label='Z')
plt.legend()

plt.title('PSP position ECLIPJ2000')
plt.ylabel('AU')


################## Earth position in ECLIJ2000

print()
print()
earth=spice.Trajectory('399')  
starttime =datetime(2018, 10, 21)
endtime = datetime(2025, 8, 1)
earth_time = []
while starttime < endtime:
    earth_time.append(starttime)
    starttime += timedelta(days=1)
    
earth_time_num=mdates.date2num(earth_time)     
#frames https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html  Appendix. ``Built in'' Inertial Reference Frames
earth.generate_positions(earth_time,'Sun','ECLIPJ2000')
earth.change_units(astropy.units.AU)  
print('Earth')
print(earth.r)
print(earth.x)
print(earth.y)
print(earth.z)  #zero anyway

print()
print()



################ CHANGE PSP coordinates FROM ECLIPJ2000 to HEE with Earth position
pspx_hee=np.zeros(len(earth_time))
pspy_hee=np.zeros(len(earth_time))
pspz_hee=np.zeros(len(earth_time))

for i in np.arange(np.size(earth_time)):
   ex=[np.array(earth.x[i]),np.array(earth.y[i]),0]    
   ez=np.array([0,0,1]) 
   ey=np.cross(ez,ex)  
   
   pspvec=[np.array(psp.x[i]),np.array(psp.y[i]),np.array(psp.z[i])]  
   pspx_hee[i]=np.dot(pspvec,ex) 
   pspy_hee[i]=np.dot(pspvec,ey) 
   pspz_hee[i]=np.dot(pspvec,ez) 
   
print('PSP HEE') 
print(pspx_hee)
print(pspy_hee)
print(pspz_hee)
[psp_r_hee, psp_theta_hee, psp_phi_hee]=cart2sphere(pspx_hee,pspy_hee,pspz_hee)






###################BepiColombo


starttime =datetime(2018, 10, 21)
endtime = datetime(2025, 8, 1)
bepi_time = []
while starttime < endtime:
    bepi_time.append(starttime)
    starttime += timedelta(days=1)

#    -68         'BEPICOLOMBO MMO'
print('Bepi')

spice.furnish('/Users/chris/heliopy/data/spice/bc_mpo_fcp_00040_20181020_20251102_v01.bsp')

bepi=spice.Trajectory('BEPICOLOMBO MPO')

bepi.generate_positions(bepi_time,'Sun','J2000')
#bepie.generate_positions(time2,'Sun','ECLIPJ2000')

bepi.change_units(astropy.units.AU)  

print(bepi.r)
print(bepi.x)
print(bepi.y)
print(bepi.z)  

plt.figure(2)
plt.plot_date(bepi_time,bepi.r,'-')
plt.plot_date(bepi_time,bepi.x,'-')
plt.plot_date(bepi_time,bepi.y,'-')
plt.plot_date(bepi_time,bepi.z,'-')



[bepi_r, bepi_theta, bepi_phi]=cart2sphere(bepi.x,bepi.y,bepi.z)


################ CHANGE bepi coordinates FROM ECLIPJ2000 to HEE with Earth position
bepix_hee=np.zeros(len(earth_time))
bepiy_hee=np.zeros(len(earth_time))
bepiz_hee=np.zeros(len(earth_time))

for i in np.arange(np.size(earth_time)):
   ex=[np.array(earth.x[i]),np.array(earth.y[i]),0]    
   ez=np.array([0,0,1]) 
   ey=np.cross(ez,ex)  
   
   bepivec=[np.array(bepi.x[i]),np.array(bepi.y[i]),np.array(bepi.z[i])]  
   bepix_hee[i]=np.dot(bepivec,ex) 
   bepiy_hee[i]=np.dot(bepivec,ey) 
   bepiz_hee[i]=np.dot(bepivec,ez) 
   
print('bepi HEE') 
print(bepix_hee)
print(bepiy_hee)
print(bepiz_hee)
[bepi_r_hee, bepi_theta_hee, bepi_phi_hee]=cart2sphere(bepix_hee,bepiy_hee,bepiz_hee)






###################Solar Orbiter

###############################################################################
# Load the solar orbiter spice kernel. HelioPy will automatically fetch the
# latest kernel
orbiter_kernel = spicedata.get_kernel('solo_2020')
spice.furnish(orbiter_kernel)
orbiter = spice.Trajectory('Solar Orbiter')

###############################################################################
# Generate a time for every day between starttime and endtime
starttime = datetime(2020, 3, 1)
endtime = datetime(2026, 1, 1)
solo_time = []
while starttime < endtime:
    solo_time.append(starttime)
    starttime += timedelta(days=1)

###############################################################################
# Generate positions
orbiter.generate_positions(solo_time, 'Sun', 'ECLIPJ2000')
orbiter.change_units(astropy.units.AU)

plt.figure(3, figsize=(12,8))
plt.plot_date(psp_time,psp.r,'-',label='PSP')
plt.plot_date(bepi_time,bepi.r,'-',label='Bepi Colombo')
plt.plot_date(solo_time,orbiter.r,'-',label='Solar Orbiter')
plt.legend()
plt.title('Heliocentric distance of heliospheric observatories')
plt.ylabel('AU')
plt.savefig('bepi_psp_solo.png')



######################## Animation


frame_time_num=mdates.date2num(sunpy.time.parse_time('2020-Jan-25 18:00:00'))

frame_time_num=mdates.date2num(sunpy.time.parse_time('2021-Apr-29 00:00:00'))

#frame_time_num=mdates.date2num(sunpy.time.parse_time('2020-Jun-03 00:00:00'))

#frame_time_num=mdates.date2num(sunpy.time.parse_time('2024-Dec-25 18:00:00'))


k=0

sns.set_context('talk')
sns.set_style('darkgrid')

plt.figure(4, figsize=(12,9), dpi=100, facecolor='w', edgecolor='w')

ax = plt.subplot(111,projection='polar')

dct=frame_time_num+k-pos_time_num
#get index of closest to 0, use this for position
timeind=np.argmin(abs(dct))


dct=frame_time_num+k-psp_time_num
#get index of closest to 0, use this for position
psp_timeind=np.argmin(abs(dct))

frame_time_str=str(mdates.num2date(frame_time_num+k))
print( 'current frame_time_num', frame_time_str)

#index 1 is longitude, 0 is rdist
symsize=100
ax.scatter(pos.venus[1,timeind], pos.venus[0,timeind], s=symsize, c='orange', alpha=1, lw=0, zorder=3)
ax.scatter(pos.mercury[1,timeind], pos.mercury[0,timeind], s=symsize, c='dimgrey', alpha=1,lw=0, zorder=3)
ax.scatter(pos.sta[1,timeind], pos.sta[0,timeind], s=symsize, c='red', alpha=1,marker='s',lw=0,zorder=3)
ax.scatter(pos.earth[1,timeind], pos.earth[0,timeind], s=symsize, c='mediumseagreen', alpha=1,lw=0,zorder=3)
ax.scatter(pos.mars[1,timeind], pos.mars[0,timeind], s=symsize, c='orangered', alpha=1,lw=0,zorder=3)
#ax.scatter(psp_phi[psp_timeind], psp_r[psp_timeind], s=symsize, c='black', marker='s', alpha=1,lw=0,zorder=3)
ax.scatter(psp_phi_hee[psp_timeind], psp_r_hee[psp_timeind], s=symsize, c='black', marker='s', alpha=1,lw=0,zorder=3)
ax.scatter(bepi_phi_hee[psp_timeind], bepi_r_hee[psp_timeind], s=symsize, c='blue', marker='s', alpha=1,lw=0,zorder=3)
 

 
plt.suptitle('Spacecraft trajectories')	
 
#Sun
ax.scatter(0,0,s=100,c='yellow',alpha=0.8, edgecolors='yellow')
plt.figtext(0.51,0.5,'Sun', fontsize=10, ha='center')

#Earth
plt.figtext(0.51,0.28,'Earth', fontsize=10, ha='center')

#units
plt.figtext(0.525,0.0735,'HEEQ longitude', fontsize=10, ha='left')
#	plt.figtext(0.64,0.213,'AU', fontsize=10, ha='center')

#----------------- legend

plt.figtext(0.1-0.02,0.02,'Mercury', color='dimgrey', ha='center')
plt.figtext(0.2-0.02	,0.02,'Venus', color='orange', ha='center')
plt.figtext(0.3-0.02,0.02,'Earth', color='mediumseagreen', ha='center')
plt.figtext(0.4-0.02,0.02,'Mars', color='orangered', ha='center')
plt.figtext(0.5-0.02,0.02,'STEREO-A', color='red', ha='center')
plt.figtext(0.8-0.02,0.02,'Parker Solar Probe', color='black', ha='center')

 
#set axes
plt.thetagrids(range(0,360,45),(u'0\u00b0',u'45\u00b0',u'90\u00b0',u'135\u00b0',u'+/- 180\u00b0',u'- 135\u00b0',u'- 90\u00b0',u'- 45\u00b0'), fmt='%d')
ax.set_theta_zero_location('S')
plt.rgrids((0.25,0.5,0.75, 1.0,1.25, 1.5, 1.75, 2.0),('0.25','0.5','0.75','1.0','1.25','1.5','1.75','2.0 AU'))
ax.set_ylim(0, 1.5)

#plot text for date extra so it does not move 
#year
plt.figtext(0.47,0.85,frame_time_str[0:4], fontsize=13, ha='center')
#month
plt.figtext(0.51,0.85,frame_time_str[5:7], fontsize=13, ha='center')
#day
plt.figtext(0.54,0.85,frame_time_str[8:10], fontsize=13, ha='center')
#hours
plt.figtext(0.57,0.85,frame_time_str[11:13], fontsize=13, ha='center')







