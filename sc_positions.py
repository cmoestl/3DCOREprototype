# Spacecraft trajectories for 3DCORE incl. Bepi Colombo, PSP, Solar Orbiter

#Author: C. Moestl, IWF Graz, Austria
#twitter @chrisoutofspace, https://github.com/cmoestl
#started December 2018

#needs python 3.6 with sunpy, heliopy, astropy, spiceypy, cdflib, seaborn, urllib, json, pickle, sys

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








#TO DO: 3d plotting for using latitude correctly
#positions sta, mercury, venus with python spice


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
import time

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
  print( 'done reading IDL SPICE positions')
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



def load_url_current_directory(filename,url):
#loads a file from any url to the current directory
#I use owncloud for the direct url links, 
#also works for dropbox when changing the last 0 to 1 in the url-> gives a direct link to files

 if not os.path.exists(filename):
  print('download file ', filename, ' from')
  print(url)
  try: 
    urllib.request.urlretrieve(url, filename)
    print('done')
  except urllib.error.URLError as e:
    print(' ', data_url,' ',e.reason)



##################################################### MAIN ###############################


start=time.time()



res_in_days=0.5

#1 hour res
#res_in_days=1/24.

#1min
res_in_days=1/3600.

if res_in_days < 0.2: high_res_mode=True
else:high_res_mode=False


outputdirectory='positions_animation'

if high_res_mode:
   outputdirectory='positions_animation_flyby_high_res'


#get positions old style for comparison from IDL SPICE:
pos=getcat('DATACAT/positions_2007_2023_HEEQ_6hours.sav')
pos_time_num=time_to_num_cat(pos.time)[0]

#https://heliopy.readthedocs.io/en/stable/api/heliopy.spice.Trajectory.html#heliopy.spice.Trajectory.generate_positions

#load kernels if not here
psp_kernel=spicedata.get_kernel('psp')
solo_kernel=spicedata.get_kernel('solo_2020')
planet_kernel=spicedata.get_kernel('planet_trajectories')
#sta_kernel=spicedata.get_kernel('stereo_a')
#bepi colombo is loaded directly from file

#for PSP NAIF CODE is -96 (search for solar probe plus)
#-144        'SOLAR ORBITER'
#Earth 399






##########################################  PSP

spice.furnish(psp_kernel)
psp=spice.Trajectory('-96')

starttime =datetime(2018, 8,13)
endtime = datetime(2025, 8, 31)
psp_time = []
while starttime < endtime:
    psp_time.append(starttime)
    starttime += timedelta(days=res_in_days)
    
psp_time_num=mdates.date2num(psp_time)     
#frames https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html  Appendix. ``Built in'' Inertial Reference Frames
psp.generate_positions(psp_time,'Sun','ECLIPJ2000')
psp.change_units(astropy.units.AU)  
print('PSP')
print(psp.r)
print(psp.x)
print(psp.y)
print(psp.z)  


[psp_r, psp_theta, psp_phi]=cart2sphere(psp.x,psp.y,psp.z)


#plt.figure(1, figsize=(12,9))
#plt.plot_date(psp_time,psp_r,'-', label='R')
#plt.plot_date(psp_time,psp_theta,'-',label='theta')
#plt.plot_date(psp_time,psp_phi,'-',label='phi')
#plt.title('PSP position ECLIPJ2000')
#plt.ylabel('AU')


#plt.plot_date(psp_time,psp.x,'-',label='X')
#plt.plot_date(psp_time,psp.y,'-',label='Y')
#plt.plot_date(psp_time,psp.z,'-',label='Z')
#plt.legend()



################## Earth position in ECLIJ2000

print()
print()
earth=spice.Trajectory('399')  
starttime =datetime(2018, 8, 13)
endtime = datetime(2025, 8, 31)
earth_time = []
while starttime < endtime:
    earth_time.append(starttime)
    starttime += timedelta(days=res_in_days)
    
earth_time_num=mdates.date2num(earth_time)     
#frames https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html  Appendix. ``Built in'' Inertial Reference Frames
earth.generate_positions(earth_time,'Sun','ECLIPJ2000')
earth.change_units(astropy.units.AU)  
print('Earth')
print(earth.r)
print(earth.x)
print(earth.y)
print(earth.z)  #is zero anyway

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
[psp_r_hee, psp_lat_hee, psp_lon_hee]=cart2sphere(pspx_hee,pspy_hee,pspz_hee)


plt.figure(1, figsize=(12,9))
plt.plot_date(psp_time,psp_r_hee,'-', label='R')
plt.plot_date(psp_time,psp_lat_hee,'-',label='lat')
plt.plot_date(psp_time,psp_lon_hee,'-',label='lon')
plt.title('PSP position HEE')
plt.ylabel('AU / RAD')
plt.legend()



############ generate PSP trajectory high res end of January 2020 flyby
















print()
print()
print()

############################################## BepiColombo


starttime =datetime(2018, 10, 21)
endtime = datetime(2025, 11, 2)
bepi_time = []
while starttime < endtime:
    bepi_time.append(starttime)
    starttime += timedelta(days=res_in_days)

bepi_time_num=mdates.date2num(bepi_time) 


#    -68         'BEPICOLOMBO MMO'
print('Bepi')

#download from url to local directory: https://repos.cosmos.esa.int/socci/projects/SPICE_KERNELS/repos/bepicolombo/browse/kernels/spk
spice.furnish('bc_mpo_fcp_00054_20181020_20251102_v01.bsp')

bepi=spice.Trajectory('BEPICOLOMBO MPO')
bepi.generate_positions(bepi_time,'Sun','ECLIPJ2000')
bepi.change_units(astropy.units.AU)  



#[bepi_r, bepi_theta, bepi_phi]=cart2sphere(bepi.x,bepi.y,bepi.z)


################ CHANGE bepi coordinates FROM ECLIPJ2000 to HEE with Earth position




print()
print()
earth=spice.Trajectory('399')  
starttime =datetime(2018, 10, 21)
endtime = datetime(2025, 11, 2)
earth_time = []
while starttime < endtime:
    earth_time.append(starttime)
    starttime += timedelta(days=res_in_days)
    
earth_time_num=mdates.date2num(earth_time)     
#frames https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html  Appendix. ``Built in'' Inertial Reference Frames
earth.generate_positions(earth_time,'Sun','ECLIPJ2000')
earth.change_units(astropy.units.AU)  
print('Earth')
print(earth.r)
print(earth.x)
print(earth.y)
print(earth.z)  #is zero anyway

print()
print()






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
[bepi_r_hee, bepi_lat_hee, bepi_lon_hee]=cart2sphere(bepix_hee,bepiy_hee,bepiz_hee)



plt.figure(2, figsize=(12,9))
plt.plot_date(bepi_time,bepi_r_hee,'-', label='R')
plt.plot_date(bepi_time,bepi_lat_hee,'-',label='lat')
plt.plot_date(bepi_time,bepi_lon_hee,'-',label='lon')
plt.title('bepi position HEE')
plt.ylabel('AU / RAD')
plt.legend()






#################################################### Solar Orbiter

# Load the solar orbiter spice kernel. HelioPy will automatically fetch the
# latest kernel
spice.furnish(solo_kernel)
orbiter = spice.Trajectory('Solar Orbiter')

# Generate a time for every day between starttime and endtime
starttime = datetime(2020, 3, 1)
endtime = datetime(2026, 1, 1)
solo_time = []
while starttime < endtime:
    solo_time.append(starttime)
    starttime += timedelta(days=res_in_days)

solo_time_num=mdates.date2num(solo_time) 


solo=spice.Trajectory('Solar Orbiter')
# Generate positions
solo.generate_positions(solo_time, 'Sun', 'ECLIPJ2000')
solo.change_units(astropy.units.AU)

################ CHANGE solo coordinates FROM ECLIPJ2000 to HEE with Earth position
print()
print()
earth=spice.Trajectory('399')  
starttime =datetime(2020, 3, 1)
endtime = datetime(2026, 1, 1)
earth_time = []
while starttime < endtime:
    earth_time.append(starttime)
    starttime += timedelta(days=res_in_days)
    
earth_time_num=mdates.date2num(earth_time)     
#frames https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html  Appendix. ``Built in'' Inertial Reference Frames
earth.generate_positions(earth_time,'Sun','ECLIPJ2000')
earth.change_units(astropy.units.AU)  

solox_hee=np.zeros(len(earth_time))
soloy_hee=np.zeros(len(earth_time))
soloz_hee=np.zeros(len(earth_time))

for i in np.arange(np.size(earth_time)):
   ex=[np.array(earth.x[i]),np.array(earth.y[i]),0]    
   ez=np.array([0,0,1]) 
   ey=np.cross(ez,ex)  
   
   solovec=[np.array(solo.x[i]),np.array(solo.y[i]),np.array(solo.z[i])]  
   solox_hee[i]=np.dot(solovec,ex) 
   soloy_hee[i]=np.dot(solovec,ey) 
   soloz_hee[i]=np.dot(solovec,ez) 
   
print('solo HEE') 
print(solox_hee)
print(soloy_hee)
print(soloz_hee)
[solo_r_hee, solo_lat_hee, solo_lon_hee]=cart2sphere(solox_hee,soloy_hee,soloz_hee)



plt.figure(3, figsize=(12,9))
plt.plot_date(solo_time,solo_r_hee,'-', label='R')
plt.plot_date(solo_time,solo_lat_hee,'-',label='lat')
plt.plot_date(solo_time,solo_lon_hee,'-',label='lon')
plt.title('solo position HEE')
plt.ylabel('AU / RAD')
plt.legend()









plt.figure(4, figsize=(16,10))
plt.plot_date(psp_time,psp.r,'-',label='PSP')
plt.plot_date(bepi_time,bepi.r,'-',label='Bepi Colombo')
plt.plot_date(solo_time,solo.r,'-',label='Solar Orbiter')
plt.legend()
plt.title('Heliocentric distance of heliospheric observatories')
plt.ylabel('AU')
plt.savefig('bepi_psp_solo_R.png')




plt.figure(5, figsize=(16,10))
plt.plot_date(psp_time,psp_lon_hee*180/np.pi,'-',label='PSP')
plt.plot_date(bepi_time,bepi_lon_hee*180/np.pi,'-',label='Bepi Colombo')
plt.plot_date(solo_time,solo_lon_hee*180/np.pi,'-',label='Solar Orbiter')
plt.legend()
plt.title('HEE Longitude')
plt.ylabel('DEG')
plt.savefig('bepi_psp_solo_longitude.png')



#save hee
if high_res_mode:
 pickle.dump([psp_time,psp_time_num,psp_r_hee,psp_lon_hee,psp_lat_hee,bepi_time,bepi_time_num,bepi_r_hee,bepi_lon_hee,bepi_lat_hee,solo_time,solo_time_num,solo_r_hee,solo_lon_hee,solo_lat_hee], open( "psp_solo_bepi_high_res_HEE_1min.p", "wb" ) )
 sys.exit()

#load
#psp_time,psp_time_num,psp_r_hee,psp_lon_hee,psp_lat_hee,bepi_time,bepi_time_num,bepi_r_hee,bepi_lon_hee,bepi_lat_hee,solo_time,solo_time_num,solo_r_hee,solo_lon_hee,solo_lat_hee]==pickle.load( open( "psp_solo_bepi_high_res_HEE.p", "rb" ) )

end=time.time()
print( 'generate position in HEE took time in seconds:', (end-start) )




######################## Animation


print()
print()
print()


print('make animation')


#**************** in 3D plotten und von oben anschauen wg latitude!!

#from psp start
frame_time_num=mdates.date2num(sunpy.time.parse_time('2018-Aug-13 00:00:00'))
kend=365*2*7

#frame_time_num=mdates.date2num(sunpy.time.parse_time('2021-Apr-29 00:00:00'))
#frame_time_num=mdates.date2num(sunpy.time.parse_time('2020-Jun-03 00:00:00'))
#frame_time_num=mdates.date2num(sunpy.time.parse_time('2024-Dec-25 18:00:00'))




#high res flyby
if high_res_mode:
 frame_time_num=mdates.date2num(sunpy.time.parse_time('2020-Jan-20 00:00:00'))
 kend=500


k=0


if os.path.isdir(outputdirectory) == False: os.mkdir(outputdirectory)

sns.set_context('talk')
sns.set_style('darkgrid')

plt.figure(6, figsize=(14,10), dpi=100, facecolor='w', edgecolor='w')



#lowres
#for k in np.arange(0,2000):

#highres
for k in np.arange(0,kend):


 ax = plt.subplot(111,projection='polar') #eigentlich 3D bzw projection to 2D...
 plt.suptitle('Spacecraft trajectories')	


 dct=frame_time_num+k*res_in_days-pos_time_num
 #get index of closest to 0, use this for position
 timeind=np.argmin(abs(dct))


 dct=frame_time_num+k*res_in_days-psp_time_num
 psp_timeind=np.argmin(abs(dct))

 dct=frame_time_num+k*res_in_days-bepi_time_num
 bepi_timeind=np.argmin(abs(dct))


 dct=frame_time_num+k*res_in_days-solo_time_num
 solo_timeind=np.argmin(abs(dct))



 frame_time_str=str(mdates.num2date(frame_time_num+k*res_in_days))
 print( 'current frame_time_num', frame_time_str)

 #index 1 is longitude, 0 is rdist
 symsize_planet=70
 symsize_spacecraft=40

 ax.scatter(pos.venus[1,timeind], pos.venus[0,timeind], s=symsize_planet, c='orange', alpha=1, lw=0, zorder=3)
 ax.scatter(pos.mercury[1,timeind], pos.mercury[0,timeind], s=symsize_planet, c='dimgrey', alpha=1,lw=0, zorder=3)
 ax.scatter(pos.earth[1,timeind], pos.earth[0,timeind], s=symsize_planet, c='mediumseagreen', alpha=1,lw=0,zorder=3)
 ax.scatter(pos.mars[1,timeind], pos.mars[0,timeind], s=symsize_planet, c='orangered', alpha=1,lw=0,zorder=3)
 
 
 
 ax.scatter(pos.sta[1,timeind], pos.sta[0,timeind], s=symsize_spacecraft, c='red', alpha=1,marker='s',lw=0,zorder=3)
 #ax.scatter(psp_phi[psp_timeind], psp_r[psp_timeind], s=symsize, c='black', marker='s', alpha=1,lw=0,zorder=3)
 ax.scatter(psp_lon_hee[psp_timeind], psp_r_hee[psp_timeind], s=symsize_spacecraft, c='black', marker='s', alpha=1,lw=0,zorder=3)
 
 if bepi_timeind > 0:
   ax.scatter(bepi_lon_hee[bepi_timeind], bepi_r_hee[bepi_timeind], s=symsize_spacecraft, c='blue', marker='s', alpha=1,lw=0,zorder=3)

 if solo_timeind > 0:
   ax.scatter(solo_lon_hee[solo_timeind], solo_r_hee[solo_timeind], s=symsize_spacecraft, c='green', marker='s', alpha=1,lw=0,zorder=3)
 
 
 fsize=10
 
 #Sun
 ax.scatter(0,0,s=100,c='yellow',alpha=0.8, edgecolors='yellow')
 #plt.figtext(0.51,0.5,'Sun', fontsize=fsize, ha='center')

 #Earth
 #plt.figtext(0.8,0.28,'Earth', fontsize=fsize, ha='center')

 #units
 plt.figtext(0.525,0.0735,'HEE longitude', fontsize=fsize, ha='left')
 #	plt.figtext(0.64,0.213,'AU', fontsize=10, ha='center')

 #----------------- legend

 plt.figtext(0.1-0.02,0.02,'Mercury', color='dimgrey', ha='center')
 plt.figtext(0.2-0.02	,0.02,'Venus', color='orange', ha='center')
 plt.figtext(0.3-0.02,0.02,'Earth', color='mediumseagreen', ha='center')
 plt.figtext(0.4-0.02,0.02,'Mars', color='orangered', ha='center')
 plt.figtext(0.5-0.02,0.02,'STEREO-A', color='red', ha='center')
 plt.figtext(0.65-0.02,0.02,'PSP', color='black', ha='center')
 plt.figtext(0.77-0.02,0.02,'Bepi Colombo', color='blue', ha='center')
 plt.figtext(0.9-0.02,0.02,'Solar Orbiter', color='green', ha='center')


 
 #set axes
 plt.thetagrids(range(0,360,45),(u'0\u00b0',u'45\u00b0',u'90\u00b0',u'135\u00b0',u'+/- 180\u00b0',u'- 135\u00b0',u'- 90\u00b0',u'- 45\u00b0'), fmt='%d')
 ax.set_theta_zero_location('S')
 plt.rgrids((0.25,0.5,0.75, 1.0,1.25, 1.5, 1.75, 2.0),('0.25','0.5','0.75','1.0','1.25 AU','1.5','1.75','2.0'), fontsize=10)
 ax.set_ylim(0, 1.35)

 #plot text for date extra so it does not move 
 #year
 plt.figtext(0.47,0.85,frame_time_str[0:4], fontsize=13, ha='center')
 #month
 plt.figtext(0.51,0.85,frame_time_str[5:7], fontsize=13, ha='center')
 #day
 plt.figtext(0.54,0.85,frame_time_str[8:10], fontsize=13, ha='center')
 #hours
 plt.figtext(0.57,0.85,frame_time_str[11:13], fontsize=13, ha='center')

 

 #save figure
 framestr = '%04i' % (k)  
 filename=outputdirectory+'/pos_anim_'+framestr+'.jpg'  
 plt.savefig(filename)

 
 plt.clf()
 
 
 
#######################
 
print('anim done')
 
 
os.system('/Users/chris/python/3DCORE/ffmpeg -r 30 -i /Users/chris/python/3DCORE/positions_animation/pos_anim_%04d.jpg -b 5000k -r 30 pos_anim.mp4 -y -loglevel quiet')


os.system('/Users/chris/python/3DCORE/ffmpeg -r 90 -i /Users/chris/python/3DCORE/positions_animation_flyby_high_res/pos_anim_%04d.jpg -b 5000k -r 90 pos_anim_flyby_high_res.mp4 -y -loglevel quiet')


print('movie done')

