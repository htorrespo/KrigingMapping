#!/usr/bin/env python
# Filename: createMap.py
#
# History
# - v_1: create Kriging map wrapping R code in python 
#
#
# This code creates a 2d map using kriging from an input of spatially        #
# scattered data points. In the example, the points are some random          #
# samplings of the NGC 1023 velocity dispersion field.                       #
# A more complete description of the procedure can be found in Appendix A of #
# Pastorello+(2014) (http://adsabs.harvard.edu/abs/2014MNRAS.442.1003P).     #

#
# In order to run the code, the following python packages are requested:     #
#
# -asciidata
# -collections
# -matplotlib
# -numpy
# -os
# -pickle
# -pyfits
# -pylab
# -pyraf
# -random
# -rpy2
# -scipy
# -sys
# -time
# 
#
# Copyright: Nicola Pastorello (2015)
#
###############################################

# Create kriging map 
from Nicola import *

#Retrieve dictionary
from KrigingMapping_def_v4 import *
from createMaps__def__ import *



### TEST DATA INPUT
Xextra, Yextra, Zextra, eZextra = numpy.transpose(numpy.loadtxt('./testData.txt'))

for ii in numpy.arange(len(Xextra)):
  if not(numpy.isnan(Zextra[ii])):
    X.append(Xextra[ii])
    Y.append(Yextra[ii])
    Z.append(Zextra[ii])
    eZ.append(eZextra[ii])

genTable = transpose(numpy.array([X, Y, Z, eZ]))


#Saving new files
fileout = open('./listElements.txt', 'wb')
numpy.savetxt(fileout, genTable, delimiter='\t', header='x\ty\tz\terrz')
fileout.close()



### KRIGING PARAMETERS 

theta = 10. # Kriging range
coeff = 3. # Kriging coeff

dummy = KrigingR('./listElements.txt', visualize=False, 
         theta_r = theta, coeff_r = coeff, savePdf = True, 
         pathOutput = './', label='Kriging', sizePixelMap=300) 

dummy = KrigingMapPython('./', 'GalaxyName', genTable, label='Kriging',
                            limits = [0, 250], sizePixelMap=300)  #For the visualization

print 'DONE'