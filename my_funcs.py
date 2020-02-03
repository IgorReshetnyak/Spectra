#COPYRIGHT Igor RESHETNYAK 2016-2020 

import math
import numpy
import time

def broad_gauss(data, energies, broad=0.1): 
    if len(data)!=len(energies): 
     exit("Error data and energies length differ") 
    data1=numpy.zeros(len(energies)); 
    ne=len(energies) 
    for i in range(ne): 
     norm=0.0 
     for i1 in range(ne): 
      if (abs(energies[i1]-energies[i])<5*broad): 
       data1[i]+=data[i1]*gauss(energies[i1],energies[i],broad) 
       norm+=gauss(energies[i1],energies[i],broad) 
      
     data1[i]=data1[i]/norm 
    return data1

def broad_lorentz(data, energies, broad=0.1): 
    if len(data)!=len(energies): 
     exit("Error data and energies length differ") 
    data1=numpy.zeros(len(energies)); 
    ne=len(energies) 
    for i in range(ne): 
     norm=0.0 
     for i1 in range(ne): 
      if (abs(energies[i1]-energies[i])<5*broad): 
       data1[i]+=data[i1]*lorentz(energies[i1],energies[i],broad) 
       norm+=lorentz(energies[i1],energies[i],broad) 
      
     data1[i]=data1[i]/norm 
    return data1

def broad_dynamic(data,energies,broad):
    if len(data)!=len(energies) or len(broad)!=len(energies): 
     exit("Error data and energies length differ") 
    data1=numpy.zeros(len(energies)); 
    ne=len(energies) 
    for i in range(ne): 
     norm=0.0 
     for i1 in range(ne): 
      if (abs(energies[i1]-energies[i])<5*broad[i]): 
       data1[i]+=data[i1]*gauss(energies[i1],energies[i],broad[i]) 
       norm+=gauss(energies[i1],energies[i],broad[i]) 
      
     data1[i]=data1[i]/norm 
    return data1

def gauss(x, mu, sigma):
 return (sigma * math.sqrt(2*math.pi))*math.exp(-1.0 / (2 * sigma * sigma) * (x - mu)*(x - mu))

def lorentz(x, mu, sigma):
 return (sigma /math.pi/2)*(1/((x-mu)*(x-mu)+sigma*sigma/4))

def printStage(txt):
  print txt
  print time.ctime()


def eelsToEps(spectra):
  n=len(spectra)

  if n<1: exit()

  if n%2==1: n=n-1
  im=spectra
  re=spectra 
  re=numpy.fft.fft(re,n)
  re=-2.*re.imag/n
  re[0:n/2]=-re[0:n/2]
  re=numpy.fft.fft(re,n)
  re=re.real
  re=re[1:n]
  im=im[1:n]

  e1 = re / (re ** 2 + im ** 2)
  e2 = im / (re ** 2 + im ** 2)    
  return e2
