# Spacecraft and planet trajectories in numpy incl. Bepi Colombo, PSP, Solar Orbiter

#Author: C. Moestl, IWF Graz, Austria
#twitter @chrisoutofspace, https://github.com/cmoestl
#December 2018 - March 2019

#needs python 3.6 with sunpy, heliopy, astropy, spiceypy, cdflib, seaborn, urllib, json, pickle, sys

#change path for ffmpeg for animation production at the very end

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


#TO DO: faded trajectories, make structured arrays for saving, get_sc_lonlat to be made into an universal function






#import scipy.io
import os
import datetime
from datetime import datetime, timedelta
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pdb
import sunpy.time
import pickle
import seaborn as sns
import sys
import heliopy.data.spice as spicedata
import heliopy.spice as spice
import astropy
import time
import numba
from numba import jit

#ignore warnings
#import warnings
#warnings.filterwarnings('ignore')

##########################################################################################
######################################### CODE START #####################################
##########################################################################################



@jit(nopython=True)
def sphere2cart(r, phi, theta):
    x = r*np.cos(theta)*np.cos(phi)
    y = r*np.cos(theta)*np.sin(phi)
    z = r*np.sin(theta)
    return (x, y, z) 

@jit(nopython=True)
def cart2sphere(x,y,z):
    r = np.sqrt(x**2+ y**2 + z**2)            # r
    theta = np.arctan2(z,np.sqrt(x**2+ y**2))     # theta
    phi = np.arctan2(y,x)                        # phi
    return (r, theta, phi)


#def get_sc_lonlat(starttime,endtime, kernel,frame):
def get_sc_lonlat_test(starttime, endtime,res_in_days):

 '''
 make spacecraft positions

 kernels: psp_pred, stereoa_pred
 frames: ECLIPJ2000 HEE
 frames https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html  Appendix. ``Built in'' Inertial Reference Frames

 '''
 spice.furnish(spicedata.get_kernel('psp_pred'))
 psp=spice.Trajectory('SPP')
 psp_time = []
 while starttime < endtime:
    psp_time.append(starttime)
    starttime += timedelta(days=res_in_days)
    
 psp.generate_positions(psp_time,'Sun','HEEQ')
 psp.change_units(astropy.units.AU)  
 psp_r, psp_lat_heeq, psp_lon_heeq=cart2sphere(psp.x,psp.y,psp.z)
 
 return psp_r, psp_lat_heeq, psp_lon_heeq


##################################################### MAIN ###############################


############# SETTINGS

start=time.time()

res_in_days=1

#1 hour res
#res_in_days=1/24.

#1min
#res_in_days=1/1440.

if res_in_days < 0.2: high_res_mode=True
else:high_res_mode=False

outputdirectory='positions_animation'
if high_res_mode:
   outputdirectory='positions_animation_flyby_high_res'

#test function
#starttime =datetime(2018, 8,13)
#endtime =datetime(2025, 8,31)
#a=get_sc_lonlat_test(starttime, endtime,0.5)
#end=time.time()
#print( 'generate position in HEEQ took time in seconds:', (end-start) )
#sys.exit()


################## MAKE TRAJECTORIES


##########################################  PSP

starttime =datetime(2018, 8,13)
endtime = datetime(2025, 8, 31)
psp_time = []
while starttime < endtime:
    psp_time.append(starttime)
    starttime += timedelta(days=res_in_days)
psp_time_num=mdates.date2num(psp_time)     

spice.furnish(spicedata.get_kernel('psp_pred'))
psp=spice.Trajectory('SPP')
psp.generate_positions(psp_time,'Sun','HEEQ')
psp.change_units(astropy.units.AU)  
[psp_r, psp_lat_heeq, psp_lon_heeq]=cart2sphere(psp.x,psp.y,psp.z)
print('PSP HEEQ done')


plt.figure(1, figsize=(12,9))
plt.plot_date(psp_time,psp_r,'-', label='R')
plt.plot_date(psp_time,psp_lat_heeq,'-',label='lat')
plt.plot_date(psp_time,psp_lon_heeq,'-',label='lon')
plt.title('PSP position HEEQ')
plt.ylabel('AU / RAD')
plt.legend()


############################################## BepiColombo

starttime =datetime(2018, 10, 21)
endtime = datetime(2025, 11, 2)
bepi_time = []
while starttime < endtime:
    bepi_time.append(starttime)
    starttime += timedelta(days=res_in_days)
bepi_time_num=mdates.date2num(bepi_time) 

spice.furnish(spicedata.get_kernel('bepi_pred'))
bepi=spice.Trajectory('BEPICOLOMBO MPO') # or BEPICOLOMBO MMO
bepi.generate_positions(bepi_time,'Sun','HEEQ')
bepi.change_units(astropy.units.AU)  
[bepi_r, bepi_lat_heeq, bepi_lon_heeq]=cart2sphere(bepi.x,bepi.y,bepi.z)
print('Bepi HEEQ done')

plt.figure(2, figsize=(12,9))
plt.plot_date(bepi_time,bepi_r,'-', label='R')
plt.plot_date(bepi_time,bepi_lat_heeq,'-',label='lat')
plt.plot_date(bepi_time,bepi_lon_heeq,'-',label='lon')
plt.title('bepi position HEEQ')
plt.ylabel('AU / RAD')
plt.legend()



#################################################### Solar Orbiter

starttime = datetime(2020, 3, 1)
endtime = datetime(2026, 1, 1)
solo_time = []
while starttime < endtime:
    solo_time.append(starttime)
    starttime += timedelta(days=res_in_days)
solo_time_num=mdates.date2num(solo_time) 

spice.furnish(spicedata.get_kernel('solo_2020'))
solo=spice.Trajectory('Solar Orbiter')
solo.generate_positions(solo_time, 'Sun', 'HEEQ')
solo.change_units(astropy.units.AU)
[solo_r, solo_lat_heeq, solo_lon_heeq]=cart2sphere(solo.x,solo.y,solo.z)
print('Solo HEEQ done')

plt.figure(3, figsize=(12,9))
plt.plot_date(solo_time,solo_r,'-', label='R')
plt.plot_date(solo_time,solo_lat_heeq,'-',label='lat')
plt.plot_date(solo_time,solo_lon_heeq,'-',label='lon')
plt.title('solo position HEEQ')
plt.ylabel('AU / RAD')
plt.legend()



######## R with all three
plt.figure(4, figsize=(16,10))
plt.plot_date(psp_time,psp.r,'-',label='PSP')
plt.plot_date(bepi_time,bepi.r,'-',label='Bepi Colombo')
plt.plot_date(solo_time,solo.r,'-',label='Solar Orbiter')
plt.legend()
plt.title('Heliocentric distance of heliospheric observatories')
plt.ylabel('AU')
plt.savefig('bepi_psp_solo_R.png')

##### Longitude all three
plt.figure(5, figsize=(16,10))
plt.plot_date(psp_time,psp_lon_heeq*180/np.pi,'-',label='PSP')
plt.plot_date(bepi_time,bepi_lon_heeq*180/np.pi,'-',label='Bepi Colombo')
plt.plot_date(solo_time,solo_lon_heeq*180/np.pi,'-',label='Solar Orbiter')
plt.legend()
plt.title('HEE Longitude')
plt.ylabel('DEG')
plt.savefig('bepi_psp_solo_longitude_HEEQ.png')






############# Earth for Mercury, Venus, STA



planet_kernel=spicedata.get_kernel('planet_trajectories')

starttime =datetime(2018, 1, 1)
endtime = datetime(2028, 12, 31)
earth_time = []
while starttime < endtime:
    earth_time.append(starttime)
    starttime += timedelta(days=res_in_days)
earth_time_num=mdates.date2num(earth_time)     

earth=spice.Trajectory('399')  #399 for Earth, not barycenter (because of moon)
earth.generate_positions(earth_time,'Sun','HEEQ')
earth.change_units(astropy.units.AU)  
[earth_r, earth_lat_heeq, earth_lon_heeq]=cart2sphere(earth.x,earth.y,earth.z)
print('Earth HEEQ')

################ mercury
merc_time_num=earth_time_num
merc=spice.Trajectory('1')  #barycenter
merc.generate_positions(earth_time,'Sun','HEEQ')  
merc.change_units(astropy.units.AU)  
[merc_r, merc_lat_heeq, merc_lon_heeq]=cart2sphere(merc.x,merc.y,merc.z)
print('merc HEEQ') 

################# venus
ven_time_num=earth_time_num
ven=spice.Trajectory('2')  
ven.generate_positions(earth_time,'Sun','HEEQ')  
ven.change_units(astropy.units.AU)  
[ven_r, ven_lat_heeq, ven_lon_heeq]=cart2sphere(ven.x,ven.y,ven.z)
print('venus HEEQ') 


############### Mars

mars_time_num=earth_time_num
mars=spice.Trajectory('4')  
mars.generate_positions(earth_time,'Sun','HEEQ')  
mars.change_units(astropy.units.AU)  
[mars_r, mars_lat_heeq, mars_lon_heeq]=cart2sphere(mars.x,mars.y,mars.z)
print('mars HEEQ') 

#############stereo-A
#https://docs.heliopy.org/en/stable/data/spice.html

sta_time_num=earth_time_num
spice.furnish(spicedata.get_kernel('stereo_a_pred'))
sta=spice.Trajectory('-234')  
sta.generate_positions(earth_time,'Sun','HEEQ')  
sta.change_units(astropy.units.AU)  
[sta_r, sta_lat_heeq, sta_lon_heeq]=cart2sphere(sta.x,sta.y,sta.z)
print('sta HEEQ') 



end=time.time()
print( 'generate position in HEEQ took time in seconds:', round((end-start),1) )




#save heeq positions 
if high_res_mode:
 pickle.dump([psp_time,psp_time_num,psp_r,psp_lon_heeq,psp_lat_heeq,bepi_time,bepi_time_num,bepi_r,bepi_lon_heeq,bepi_lat_heeq,solo_time,solo_time_num,solo_r,solo_lon_heeq,solo_lat_heeq], open( "psp_solo_bepi_highres_heeq_1min.p", "wb" ) )
 sys.exit()
else: 
 pickle.dump([psp_time,psp_time_num,psp_r,psp_lon_heeq,psp_lat_heeq,bepi_time,bepi_time_num,bepi_r,bepi_lon_heeq,bepi_lat_heeq,solo_time,solo_time_num,solo_r,solo_lon_heeq,solo_lat_heeq], open( "psp_solo_bepi_lowres_heeq_12hours.p", "wb" ) )
 
# load old
#[psp_time,psp_time_num,psp_r_hee,psp_lon_hee,psp_lat_hee,bepi_time,bepi_time_num,bepi_r_hee,bepi_lon_hee,bepi_lat_hee,solo_time,solo_time_num,solo_r_hee,solo_lon_hee,solo_lat_hee]=pickle.load( open( "psp_solo_bepi_high_res_HEE_1min.p", "rb" ) )



#########################################################################################
######################## Animation

plt.close('all')


print()
print('make animation')

#from psp start
frame_time_num=mdates.date2num(sunpy.time.parse_time('2018-Aug-1 00:00:00'))
kend=int(365/res_in_days*7.4)

#frame_time_num=mdates.date2num(sunpy.time.parse_time('2021-Apr-29 00:00:00'))
#frame_time_num=mdates.date2num(sunpy.time.parse_time('2020-Jun-03 00:00:00'))
#frame_time_num=mdates.date2num(sunpy.time.parse_time('2024-Dec-25 18:00:00'))

#high res flyby
if high_res_mode:
 frame_time_num=mdates.date2num(sunpy.time.parse_time('2020-Jan-20 00:00:00'))
 kend=500


if os.path.isdir(outputdirectory) == False: os.mkdir(outputdirectory)

sns.set_context('talk')
sns.set_style('darkgrid')


plt.figure(6, figsize=(14,10), dpi=100, facecolor='w', edgecolor='w')


fsize=10


#################################################### animation loop start


for k in np.arange(0,kend):


 ax = plt.subplot(111,projection='polar') 
 frame_time_str=str(mdates.num2date(frame_time_num+k*res_in_days))
 print( 'current frame_time_num', frame_time_str, '     ',k)

 #these have their own times
 dct=frame_time_num+k*res_in_days-psp_time_num
 psp_timeind=np.argmin(abs(dct))

 dct=frame_time_num+k*res_in_days-bepi_time_num
 bepi_timeind=np.argmin(abs(dct))

 dct=frame_time_num+k*res_in_days-solo_time_num
 solo_timeind=np.argmin(abs(dct))

 #all same times
 dct=frame_time_num+k*res_in_days-earth_time_num
 earth_timeind=np.argmin(abs(dct))

 symsize_planet=70
 symsize_spacecraft=40

 #plot all positions including text R lon lat for some 
 ax.scatter(ven_lon_heeq[earth_timeind], ven_r[earth_timeind]*np.cos(ven_lat_heeq[earth_timeind]), s=symsize_planet, c='orange', alpha=1,lw=0,zorder=3)
 ax.scatter(merc_lon_heeq[earth_timeind], merc_r[earth_timeind]*np.cos(merc_lat_heeq[earth_timeind]), s=symsize_planet, c='dimgrey', alpha=1,lw=0,zorder=3)
 ax.scatter(0, earth_r[earth_timeind]*np.cos(earth_lat_heeq[earth_timeind]), s=symsize_planet, c='mediumseagreen', alpha=1,lw=0,zorder=3)
 ax.scatter(sta_lon_heeq[earth_timeind], sta_r[earth_timeind]*np.cos(sta_lat_heeq[earth_timeind]), s=symsize_spacecraft, c='red', marker='s', alpha=1,lw=0,zorder=3)
 ax.scatter(mars_lon_heeq[earth_timeind], mars_r[earth_timeind]*np.cos(mars_lat_heeq[earth_timeind]), s=symsize_planet, c='orangered', alpha=1,lw=0,zorder=3)
 
 earth_text='Earth R / lon / lat: '+str(f'{earth_r[earth_timeind]:6.2f}')+str(f'{np.zeros(1)[0]:8.1f}')+str(f'{np.rad2deg(earth_lat_heeq[earth_timeind]):8.1f}')
 f10=plt.figtext(0.01,0.9,earth_text, fontsize=fsize+2, ha='left')
 
 mars_text='Mars  R / lon / lat: '+str(f'{mars_r[earth_timeind]:6.2f}')+str(f'{np.rad2deg(mars_lon_heeq[earth_timeind]):8.1f}')+str(f'{np.rad2deg(mars_lat_heeq[earth_timeind]):8.1f}')
 f9=plt.figtext(0.01,0.86,mars_text, fontsize=fsize+2, ha='left')

 sta_text='STA   R / lon / lat: '+str(f'{sta_r[earth_timeind]:6.2f}')+str(f'{np.rad2deg(sta_lon_heeq[earth_timeind]):8.1f}')+str(f'{np.rad2deg(sta_lat_heeq[earth_timeind]):8.1f}')
 f8=plt.figtext(0.01,0.82,sta_text, fontsize=fsize+2, ha='left')


 fadeind=int(60/res_in_days)
 
 if psp_timeind > 0:
   ax.scatter(psp_lon_heeq[psp_timeind], psp_r[psp_timeind]*np.cos(psp_lat_heeq[psp_timeind]), s=symsize_spacecraft, c='black', marker='s', alpha=1,lw=0,zorder=3)
   psp_text='PSP   R / lon / lat: '+str(f'{psp_r[psp_timeind]:6.2f}')+str(f'{np.rad2deg(psp_lon_heeq[psp_timeind]):8.1f}')+str(f'{np.rad2deg(psp_lat_heeq[psp_timeind]):8.1f}')
   f5=plt.figtext(0.01,0.78,psp_text, fontsize=fsize+2, ha='left')
   ax.plot(psp_lon_heeq[psp_timeind-fadeind:psp_timeind+fadeind], psp_r[psp_timeind-fadeind:psp_timeind+fadeind]*np.cos(psp_lat_heeq[psp_timeind-fadeind:psp_timeind+fadeind]), c='black', alpha=0.2,lw=1,zorder=3)
   #for rate in np.arange(10, -1, -1)*0.1:
   #  line.set_alpha(rate)
   #  plt.draw()
   

 if bepi_timeind > 0:
   ax.scatter(bepi_lon_heeq[bepi_timeind], bepi_r[bepi_timeind]*np.cos(bepi_lat_heeq[bepi_timeind]), s=symsize_spacecraft, c='blue', marker='s', alpha=1,lw=0,zorder=3)
   bepi_text='Bepi  R / lon / lat: '+str(f'{bepi_r[bepi_timeind]:6.2f}')+str(f'{np.rad2deg(bepi_lon_heeq[bepi_timeind]):8.1f}')+str(f'{np.rad2deg(bepi_lat_heeq[bepi_timeind]):8.1f}')
   f6=plt.figtext(0.01,0.74,bepi_text, fontsize=fsize+2, ha='left')
   ax.plot(bepi_lon_heeq[bepi_timeind-fadeind:bepi_timeind+fadeind], bepi_r[bepi_timeind-fadeind:bepi_timeind+fadeind]*np.cos(bepi_lat_heeq[bepi_timeind-fadeind:bepi_timeind+fadeind]), c='blue', alpha=0.2,lw=1,zorder=3)



 if solo_timeind > 0:
   ax.scatter(solo_lon_heeq[solo_timeind], solo_r[solo_timeind]*np.cos(solo_lat_heeq[solo_timeind]), s=symsize_spacecraft, c='green', marker='s', alpha=1,lw=0,zorder=3)
   solo_text='SolO  R / lon / lat: '+str(f'{solo_r[solo_timeind]:6.2f}')+str(f'{np.rad2deg(solo_lon_heeq[solo_timeind]):8.1f}')+str(f'{np.rad2deg(solo_lat_heeq[solo_timeind]):8.1f}')
   f7=plt.figtext(0.01,0.7,solo_text, fontsize=fsize+2, ha='left')
   ax.plot(solo_lon_heeq[solo_timeind-fadeind:solo_timeind+fadeind], solo_r[solo_timeind-fadeind:solo_timeind+fadeind]*np.cos(solo_lat_heeq[solo_timeind-fadeind:solo_timeind+fadeind]), c='green', alpha=0.2,lw=1,zorder=3)


 #plot text for date extra so it does not move 
 #year
 f1=plt.figtext(0.47,0.85,frame_time_str[0:4], fontsize=13, ha='center')
 #month
 f2=plt.figtext(0.51,0.85,frame_time_str[5:7], fontsize=13, ha='center')
 #day
 f3=plt.figtext(0.54,0.85,frame_time_str[8:10], fontsize=13, ha='center')
 #hours
 f4=plt.figtext(0.57,0.85,frame_time_str[11:13], fontsize=13, ha='center')



 plt.figtext(0.1-0.02,0.02,'Mercury', color='dimgrey', ha='center')
 plt.figtext(0.2-0.02	,0.02,'Venus', color='orange', ha='center')
 plt.figtext(0.3-0.02,0.02,'Earth', color='mediumseagreen', ha='center')
 plt.figtext(0.4-0.02,0.02,'Mars', color='orangered', ha='center')
 plt.figtext(0.5-0.02,0.02,'STEREO-A', color='red', ha='center')
 plt.figtext(0.65-0.02,0.02,'PSP', color='black', ha='center')
 plt.figtext(0.77-0.02,0.02,'Bepi Colombo', color='blue', ha='center')
 plt.figtext(0.9-0.02,0.02,'Solar Orbiter', color='green', ha='center')

 plt.suptitle('Spacecraft trajectories HEEQ 2D projection  2018-2025')	
 plt.figtext(0.53,0.0735,'HEEQ longitude', fontsize=fsize+5, ha='left')
 #signature
 plt.figtext(0.97,0.01/2,r'$C. M\ddot{o}stl$', fontsize=fsize, ha='center') 
 
 #set axes

 ax.set_theta_zero_location('S')
 plt.thetagrids(range(0,360,45),(u'0\u00b0',u'45\u00b0',u'90\u00b0',u'135\u00b0',u'+/- 180\u00b0',u'- 135\u00b0',u'- 90\u00b0',u'- 45\u00b0'), fmt='%d',fontsize=fsize+5)
 plt.rgrids((0.25,0.5,0.75, 1.0,1.25, 1.5),('0.25','0.5','0.75','1.0','1.25','1.5  AU'),angle=125, fontsize=fsize)
 ax.set_ylim(0, 1.7)
 
 #Sun
 ax.scatter(0,0,s=100,c='yellow',alpha=1, edgecolors='yellow')
 

 #save figure
 framestr = '%05i' % (k)  
 filename=outputdirectory+'/pos_anim_'+framestr+'.jpg'  
 plt.savefig(filename)


 #f1.set_visible(False)
 #f2.set_visible(False)
 #f3.set_visible(False)
 #f4.set_visible(False)
 #if psp_timeind > 0: f5.set_visible(False)
 #if bepi_timeind > 0: f6.set_visible(False)
 #if solo_timeind > 0: f7.set_visible(False)
 #f8.set_visible(False)
 #f9.set_visible(False)
 #f10.set_visible(False)
 #clear
 plt.clf()

  
########################################### loop end
 
print('anim done')
 
os.system('/Users/chris/python/3DCORE/ffmpeg -r 60 -i /Users/chris/python/3DCORE/positions_animation/pos_anim_%05d.jpg -b 5000k -r 60 pos_anim.mp4 -y -loglevel quiet')
#os.system('/Users/chris/python/3DCORE/ffmpeg -r 90 -i /Users/chris/python/3DCORE/positions_animation_flyby_high_res/pos_anim_%04d.jpg -b 5000k -r 90 pos_anim_flyby_high_res.mp4 -y -loglevel quiet')

print('movie done')

