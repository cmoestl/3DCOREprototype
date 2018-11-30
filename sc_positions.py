# Spacecraft trajectories for 3DCORE incl. Bepi Colombo, PSP, Solar Orbiter

#Author: C. Moestl, IWF Graz, Austria
#twitter @chrisoutofspace, https://github.com/cmoestl
#started December 2018

#needs python > 3.6 with sunpy, heliopy, spiceypy, cdflib, seaborn

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

#https://heliopy.readthedocs.io/en/stable/api/heliopy.spice.Trajectory.html#heliopy.spice.Trajectory.generate_positions

#load kernels if not here
spicedata.get_kernel('psp')
spicedata.get_kernel('planet_trajectories')

#Solar orbiter 'solo_2020'
#for PSP NAIF CODE is -96 (search for solar probe plus)
#-144        'SOLAR ORBITER'

#Earth 399

spice.furnish(spicedata.get_kernel('psp'))
psp=spice.Trajectory('-96')

starttime =datetime(2018, 10, 1)
endtime = datetime(2025, 1, 1)
times = []
while starttime < endtime:
    times.append(starttime)
    starttime += timedelta(days=1)

#frames https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html  Appendix. ``Built in'' Inertial Reference Frames

psp.generate_positions(times,'Sun','ECLIPJ2000')
#pspe.generate_positions(time2,'Sun','ECLIPJ2000')

psp.change_units(astropy.units.AU)  


print('PSP')

print(psp.r)
print(psp.x)
print(psp.y)
print(psp.z)  


plt.figure(1)
plt.plot_date(times,psp.r,'-')
plt.plot_date(times,psp.x,'-')
plt.plot_date(times,psp.y,'-')
plt.plot_date(times,psp.z,'-')



print()

print()







###################BepiColombo


starttime =datetime(2018, 10, 21)
endtime = datetime(2025, 1, 1)
times2 = []
while starttime < endtime:
    times2.append(starttime)
    starttime += timedelta(days=1)



#    -68         'BEPICOLOMBO MMO'

print('Bepi')


spice.furnish('/Users/chris/heliopy/data/spice/bc_mpo_fcp_00040_20181020_20251102_v01.bsp')
bepi=spice.Trajectory('BEPICOLOMBO MPO')



bepi.generate_positions(times2,'Sun','J2000')
#bepie.generate_positions(time2,'Sun','ECLIPJ2000')

bepi.change_units(astropy.units.AU)  

print(bepi.r)
print(bepi.x)
print(bepi.y)
print(bepi.z)  

plt.figure(2)
plt.plot_date(times2,bepi.r,'-')
plt.plot_date(times2,bepi.x,'-')
plt.plot_date(times2,bepi.y,'-')
plt.plot_date(times2,bepi.z,'-')

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
endtime = datetime(2025, 1, 1)
times3 = []
while starttime < endtime:
    times3.append(starttime)
    starttime += timedelta(days=1)

###############################################################################
# Generate positions
orbiter.generate_positions(times3, 'Sun', 'ECLIPJ2000')
orbiter.change_units(astropy.units.AU)




plt.figure(3, figsize=(12,8))

plt.plot_date(times,psp.r,'-',label='PSP')
plt.plot_date(times2,bepi.r,'-',label='Bepi Colombo')
plt.plot_date(times3,orbiter.r,'-',label='Solar Orbiter')
plt.legend()
plt.title('Heliocentric distance of heliospheric observatories')
plt.ylabel('AU')
plt.savefig('bepi_psp_solo.png')

