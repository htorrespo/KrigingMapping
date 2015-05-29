# createMaps__def__.py
import pyfits
from Nicola import *
from math import pi

#Retrieve dictionary
lib_path = os.path.abspath('/Users/npastorello/Desktop/Galaxies/General_studies/')
sys.path.append(lib_path)
from galaxyParametersDictionary_v5 import *

def findDell(RA, Dec, PA0, b_a, unfolded=False):
  angleRot = (numpy.pi/180.)*(PA0-90.)
  xrot, yrot = (RA *numpy.cos(angleRot) - Dec * numpy.sin(angleRot), 
                RA *numpy.sin(angleRot) + Dec * numpy.cos(angleRot))
  # 
  if unfolded:
    sign = (xrot/numpy.abs(xrot))
  else:
    sign = 1.
  Rell = (b_a*(xrot**2)+(yrot**2)/b_a)**(0.5)
  #
  return sign*Rell

def retrieve_krigPos(path, galname, offset=[0.,0.], 
					 keyword='VEL', clean=True):
  inputFits = pyfits.open(path)[1].data
  #
  RAgal = convAngCoord(CentreCoordinates[galname][0])[4]*15.
  Decgal = convAngCoord(CentreCoordinates[galname][1])[4] 
  Rell = findDell(inputFits['RA']+offset[0]/3600.-RAgal, 
  					inputFits['Dec']+offset[1]/3600.-Decgal, 
  					PA0[galname], b_a[galname])
  #
  x = (inputFits['RA']-RAgal)*3600.+offset[0] 
  y = (inputFits['Dec']-Decgal)*3600.+offset[1]
  z, ez = inputFits[keyword], inputFits['ERR'+keyword]
  #
  if clean:
    selected = numpy.nonzero(inputFits['ACCEPT'] == 'ACCEPT')
    x, y, z, ez = x[selected], y[selected], z[selected], ez[selected]
  return [x, y, z, ez]

# Minor axis mask needs conversion of positions because namely the slits have 
# the same coordinates as the Major axis mask

# I'm applying the offset after the rotation
def retrieve_krigPos_Minor(path, galname, offset=[0.,0.], 
					 keyword='VEL'):
  inputFits = pyfits.open(path)[1].data
  #
  # Fix slit position in minor axis SuperSKiMS
  skypa1 = -276.70499351
  gal_PA0 = PA0[galname]
  mask_C_RA = '02:40:22.25'
  mask_C_Dec = '+39:02:48.2'
  #
  RAgal = convAngCoord(mask_C_RA)[4]*15.
  Decgal = convAngCoord(mask_C_Dec)[4]
  #
  deltaPA = mod(90.+skypa1-gal_PA0, 360.) #Difference between mask alignment and galaxy PA0
  angleNE = mod(90.-gal_PA0, 360.) #angle between galaxy major axis and East axis
  maskPA = -mod(angleNE+deltaPA, 360) #angle between the East-West axis and the mask alignment
  #
  distRA = inputFits['RA']-RA_c    #Distance from mask CentreCoordinates (in degrees)
  distDEC = inputFits['Dec']-Dec_c #Distance
  #
  angrot = maskPA*numpy.pi/180.
  realRA = (distRA*numpy.cos(angrot)-distDEC*numpy.sin(angrot))+RA_c+offset[0]/3600.   #Coordinates given the rotation of the mask
  realDEC = (distRA*numpy.sin(angrot)+distDEC*numpy.cos(angrot))+Dec_c+offset[1]/3600. #Coordinates given the rotation of the mask
  #
  Rell = findDell(realRA-RAgal, realDEC-Decgal, PA0[galname], b_a[galname])#,unfolded=True)
  #
  x, y = (realRA-RAgal)*3600., (realDEC-Decgal)*3600.
  z, ez = SS_Minor_input[keyword], SS_Minor_input['ERR'+keyword]
  #
  return [x, y, z, ez]




def drawLayout(galname='NGC1023'):
  fig = figure(0); clf(); ax = subplot(111); grid(True)
  ax.set_xlim([40.16,40.04]); ax.set_ylim([39.03,39.11])
  ax.set_aspect('equal')
  #Galaxy centre
  RAgal = convAngCoord(CentreCoordinates[galname][0])[4]*15.
  Decgal = convAngCoord(CentreCoordinates[galname][1])[4] 
  ax.scatter(RAgal, Decgal, c='r', marker='x', s=50)
  #
  #Galaxy axes
  major = ax.plot([numpy.sin(math.radians(PA0[galname]))*2+RAgal, numpy.sin(math.radians(PA0[galname]))*-2+RAgal], 
     [numpy.cos(math.radians(PA0[galname]))*2+Decgal, numpy.cos(math.radians(PA0[galname]))*-2+Decgal], 
     'r-')
  minor = ax.plot([numpy.cos(math.radians(PA0[galname]))*2+RAgal, numpy.cos(math.radians(PA0[galname]))*-2+RAgal], 
     [numpy.sin(math.radians(PA0[galname]))*-2+Decgal, numpy.sin(math.radians(PA0[galname]))*2+Decgal],
     'r-')
  #
  # Measured offset
  ax.scatter(RAgal-6./3600., Decgal+9./3600., marker='o', c='b')
  return ax



def drawMask_design(ax, angleMask=0., offsetX=0., offsetY=0):
  #fig = figure(0)
  path = '../JacobOutputs/MajorAxis/DEIMOS_pPXF_NGC1023_v0.fits'
  inputFits = pyfits.open(path)[1].data
  #
  angrot = math.radians(angleMask)
  distRA = inputFits['RA']
  distDec = inputFits['Dec']
  RA = (distRA*numpy.cos(angrot)-distDec*numpy.sin(angrot))+offsetX/3600.
  DEC = (distRA*numpy.sin(angrot)+distDec*numpy.cos(angrot))+offsetY/3600.
  #
  ax.scatter(RA, DEC, marker=(2,0,88.3), c='k')
  plt.draw()
  return True

