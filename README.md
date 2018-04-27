# 3DCORE

This code is concerned with data interpretation and prediction for space weather research. It is a method for forward modeling of solar storm magnetic flux ropes, called 3-Dimensional Coronal Rope Ejection (3DCORE). The code is able to produce synthetic in situ observations of the magnetic cores of solar coronal mass ejections sweeping over planets and spacecraft. Near Earth, these data are taken currently by the Wind, ACE and DSCOVR spacecraft. Other suitable spacecraft making these kind of observations carrying magnetometers in the solar wind were MESSENGER, Venus Express, MAVEN, and even Helios, and in the future this will be done by Parker Solar Probe, Solar Orbiter and possibly interplanetary CubeSats.

It is distributed with an MIT license and has this DOI: https://doi.org/10.6084/m9.figshare.5450341

If this research is used for peer-reviewed scientific publications, the paper below needs to cited and in any case, please also contact me either via twitter https://twitter.com/chrisoutofspace or at christian.moestl at oeaw.ac.at.

This code was used for producing the figures and results in the publication:
Möstl et al. (2018)
Forward Modeling of Coronal Mass Ejection Flux Ropes in the Inner Heliosphere with 3DCORE
in press at AGU Space Weather, https://arxiv.org/abs/1710.00587
doi: 10.1002/2017SW001735
Note that there has been a change in matplotlib from the publication to the current code version, so the 
3DCORE model volume is now plotted with ax.scatter instead of the model surface with ax.wireframe as in the paper.

Last update: 27 April 2018.

## Dependencies
* The code is written in python, and I run it over the command line in MAC OS X (High Sierra) in ipython.

* Install python anaconda v3 or up, I use 3.5.2. NumPy 1.12.1, SciPy 0.19.1, Matplotlib 2.2.2.

https://www.anaconda.com/download/#macos

* Add the packages sunpy (v0.8.5) and seaborn (v 0.8.1). 

http://docs.sunpy.org/en/stable/guide/installation/index.html

    $ conda config --add channels conda-forge
     
    $ conda install sunpy

    
https://seaborn.pydata.org/installing.html

    $ conda install seaborn    
    

## Running the code
* Download the repository, start an OS X command line and go to the "3DCORE_v2_public" directory
* Start ipython and run the only code file:

      $ ipython
      
      $ run 3dcore_v2.py
  
* After the packages are imported, the "inputfilename" and "outputdirectory" are specified in the code, and can be changed there.




