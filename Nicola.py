import scipy
import numpy
from numpy import *
import asciidata
import os
import pickle
import pyfits
import matplotlib
from matplotlib import *
matplotlib.use('TkAgg')
from matplotlib.patches import Ellipse
import pylab
from pylab import *
from sys import stdout

sens=float(.1)/3600 #Errore sull'accoppiamento per posizione in decimi arcsec o sec

#FUNCTIONS
def remEmptySpaces(initialArray):
    dim= len(initialArray) 
    k=0
    for i in range(0,dim):
        if (initialArray[i] != ' '):
            k=i
            break
    cuttedArray=initialArray[k:]
    
    return cuttedArray


def getgal():
  pathdir = os.getcwd()
  pos1, pos2 = int(1e3), int(1e3)
  for ii in range(len(pathdir)-2):
    if pathdir[ii:ii+3] == 'NGC':
      pos1 = ii
    if pathdir[ii] == '/' and ii > pos1:
      pos2 = ii
      break
  return pathdir[pos1:pos2]



def convAngCoord(angString):
    dim=len(angString)
    sep=numpy.zeros(dim)
    for i in range(0,dim):  #cerca :
        if (angString[i] == ':'):
            sep[i]=1
    pos=numpy.nonzero(sep)
    
    arcsec=float(angString[(pos[0][1])+1:])
    arcmin=float(angString[(pos[0][0]+1):(pos[0][1])])
    deg=float(angString[:(pos[0][0])])
    
    arcsectot=(deg/abs(deg))*(abs(deg*3600)+(arcmin*60)+arcsec)
    deg_tot=(deg/abs(deg))*(abs(deg)+(arcmin/60)+(arcsec/3600))
    
    return deg,arcmin,arcsec,arcsectot,deg_tot

def WeightedAverage(items, value, weight):
    numerator = sum(value[i] * weight[i] for i in items)
    divisor =  sum(weight[i] for i in items)
    return (numerator / divisor) if divisor != 0 else None


def permutation_indices(data):
    return sorted(range(len(data)), key = data.__getitem__)

def confidenceCurves(x_val,y_val,n_bin): #dim(x_val) !!= dim(y_val) 
    if (len(x_val) != len(y_val)):
        print "Different arrays' sizes!"
        return 0
    #sorting dei vettori
    zipped=zip(x_val,y_val)    
    sort_arr=sorted(zipped)
    x_sorted=[point[0] for point in sort_arr]
    y_sorted=[point[1] for point in sort_arr]
    #print x_sorted 
    #binning
    n_points=len(x_sorted)
    x_conf=[]
    y_sig_1=[]
    y_sig_2=[]
    y_sig_3=[]
    for i in range(0,int(n_points)/int(n_bin)):
        if (i==0):  #primo punto e' a x=0 
            x_conf.append(0)
        x_conf.append(x_sorted[i*n_bin+n_bin/2]) #Prendo il punto x a meta' delle x nel bin
        sum_y_bin=0 #Media y nel bin
        for j in range(0,n_bin):
            sum_y_bin=sum_y_bin+y_sorted[i*n_bin+j]
        med_y_bin=sum_y_bin/n_bin
        sum_delta_y_bin=0
        for j in range(0,n_bin):
            sum_delta_y_bin=sum_delta_y_bin+((y_sorted[i*n_bin+j]-med_y_bin)**2)
        if (i==0):  #primo punto e' a y=y0 
            y_sig_1.append(sqrt(sum_delta_y_bin/n_bin))
            y_sig_2.append(sqrt(sum_delta_y_bin/n_bin)*2)
            y_sig_3.append(sqrt(sum_delta_y_bin/n_bin)*3)
        y_sig_1.append(sqrt(sum_delta_y_bin/n_bin))
        y_sig_2.append(sqrt(sum_delta_y_bin/n_bin)*2)
        y_sig_3.append(sqrt(sum_delta_y_bin/n_bin)*3)
    y_conf=zip(y_sig_1,y_sig_2,y_sig_3)
    return x_conf, y_conf # 1xn, 3xn


#CLASSES

class LinePicker:	#Per selezionare elementi in un plot
    def __init__(self,xs,ys):

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        ax.set_title('click on a line to make it thick')

        self.handles = []
        for x,y in zip(xs,ys):
            #note use of the comma after the handle
            h, = ax.plot(x,y, '-', picker=5)
            self.handles.append(h)

        fig.canvas.mpl_connect('pick_event', self.onpick)

        self.lastartist = self.handles[-1]

        pyplot.show()

    def onpick(self,event):
        self.lastartist.set_linewidth(1)
        self.lastartist = thisline = event.artist
        thisline.set_linewidth(6)

        pyplot.draw()

