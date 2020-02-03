import struct
import math
import numpy
from numpy.linalg import norm
import scipy.linalg as linalg

class params: #Program params
      def __init__(self):
        self.label="" #Label for files
        
        self.doDiag=True #Calculate diagonal distribution
        self.de=0.1
        self.dn=10
        
        self.doFull=True #Calculate the Full matrix-matrix product
        
        self.doSelection=True #Calculate products for given vectors
        self.relative=[0.1,0.5,0.9]
        self.colors=["red","blue","black"]
        self.absolute=[]

        self.doJDOS=True
        self.broad=0.5
        self.domega=0.05

        self.showPlots=False

class my_file: #General File
    def __init__(self, fileName, resonant=False, label=""):
         self.fileName = fileName
         self.resonant = resonant
         self.label=label
         self.n=0
         self.n1=0
    def setN(self,n,n1):
         self.n=n
         self.n1=n1


class eig_file(my_file): #Eigenvalue File
     def __init__(self, fileName, resonant=False, label="",full=False,color="red"):
         my_file.__init__(self,fileName,resonant,label)
         self.eig_pairs=[]
         self.jdos=[]
         self.full=full
         self.color=color

     def appendPair(self,val,vect):
         pair=eig_pair(val,numpy.array(vect))
         self.eig_pairs.append(pair)

     def printPair(self,i):
         n1=self.n1
         pair=self.eig_pairs[i]
         pair.printSelf(n1)

     def createFromHamFile(self,ham_f):
         self.setN(ham_f.n,ham_f.n1)


         [w,vr]=linalg.eig(ham_f.ham)

         print "Max imag part",max(w[:].imag)

         wr=w[:].real*27.21

         for i in range(len(wr)):
          self.eig_pairs.append(eig_pair(wr[i],vr[:,i]))

         self.eig_pairs.sort(key=lambda pairs: pairs.value, reverse=True)
         #print len(self.eig_pairs)

     def write(self,fileNewName):

         print "Starting to write eigenvalues and eigenvectors"

         fileName=fileNewName
         f=open(fileName,mode='wb')
         n=self.n;
         n1=self.n1;
         start=4+3*4+2*4*(n1);
         size=start+(n1)*(n1+1)*8
         #data=bytes(" "*size)

         #print len(data)

         data=struct.pack(">i",4)+struct.pack(">i",n1)+struct.pack(">i",4)#+data[12:]

         f.write(data)

         #print len(data)

         print "Number of eigenvectors",n1

         w=[]

         for i in range(n1):
          w.extend([self.eig_pairs[i].value/27.21,0])
         
         data=struct.pack(">i",n1*8)+struct.pack(">"+"ff"*n1,*w)+struct.pack(">i",n1*8)

         f.write(data)

         #print len(data)

         for i in range(n1):
           starti=start+i*(n1+1)*8+8
           endi=start+(i+1)*(n1+1)*8
           eiv=[]
           for j in range(n1):
            eiv.extend([self.eig_pairs[i].vector[j].real,self.eig_pairs[i].vector[j].imag])
           #data=data[:starti]+struct.pack(">"+"ff"*n1,*eiv)+data[endi:]
           data=struct.pack(">i",n1*8)+struct.pack(">"+"ff"*n1,*eiv)+struct.pack(">i",n1*8)
           f.write(data)
         #print len(data)

         #f.write(data)
         f.close()

         print "End write"

     def read(self):
         fileName=self.fileName
         f=open(fileName,mode='rb')
         data=f.read()
         n1=struct.unpack(">i",data[4:8])[0]

         print "Starting to read eigenvalues and eigenvectors"

         if (self.resonant):
          n=n1
         else:
          n=n1/2
         self.setN(n,n1)

         print "Number of eigenvectors",n1

         w=struct.unpack(">" + "ff"*n1 , data[(4+3*4):(4+3*4+2*4*(n1))])
         eig=[w[i]*27.21 for i in range(2*n1)]
         start=4+3*4+2*4*(n1);

         for i in range(n1):
          if eig[2*i]>0 or self.full==True:
           starti=start+i*(n1+1)*8+8
           endi=start+(i+1)*(n1+1)*8
           eivw=struct.unpack(">"+"ff"*(n1),data[starti:endi])
           eiv=[complex(eivw[2*k],eivw[2*k+1]) for k in range(n1)]
           self.appendPair(eig[2*i],eiv)
          

         f.close()
         self.eig_pairs.sort(key=lambda pairs: pairs.value, reverse=True) 
         print "End read"



class eig_pair: #Pair eig_value,eig_vector
     def __init__(self, value, vector):
         self.value=value
         self.vector=vector
     def printSelf(self,n1):
          print self.value;
          vec=self.vector
          print ["%5.3f+%5.3fi"%(round(vec[j].real,3),round(vec[j].imag,3))  for j in range(n1)]

class ham_file(my_file):
    def __init__(self, fileName, resonant=False, label="",square=False,hermitian=False):
        my_file.__init__(self,fileName,resonant,label)
        self.square=square
        self.hermitian=hermitian
        self.ham=[]

    def buildFromBlocks(self,ham_h,ham_c):
      self.setN(ham_h.n1,ham_h.n1*2)
      self.ham=build_full(ham_h,ham_c)
      self.square=True
      self.hermitian=False
      self.resonant=False

    def buildFromMerge(self,ham1,ham2,dE):
      if (ham1.n1==ham2.n1):
        print "Starting build from merge"
        self.setN(ham1.n,ham1.n1)
        self.ham=numpy.zeros(shape=(ham1.n1,ham1.n1),dtype=complex)
        for i in range(self.n1):
          for j in range(i+1):
            if abs(ham1.ham[i,i]-ham1.ham[j,j])*27.21<dE:
              self.ham[i,j]=ham1.ham[i,j]
              self.ham[j,i]=ham1.ham[j,i]
            else:
              self.ham[i,j]=ham2.ham[i,j]
              self.ham[j,i]=ham2.ham[j,i]
        print "End build via merge"
      else:
        print "Size mismatch"
        exit()

    def buildFromExtend(self,ham1,dE):
      print "Starting build from extend"
      self.setN(ham1.n,ham1.n1)
      self.ham=numpy.zeros(shape=(ham1.n1,ham1.n1),dtype=complex)
      for i in range(self.n1):
          for j in range(i+1):
            if abs(ham1.ham[i,i]-ham1.ham[j,j])*27.21<=dE:
              self.ham[i,j]=ham1.ham[i,j]
              self.ham[j,i]=ham1.ham[j,i]
            else:
              self.ham[i,j]=0.0
              self.ham[j,i]=0.0
      print "End build via extend"



    def read(self):
      self.ham=numpy.zeros(shape=(self.n1,self.n1),dtype=complex)
      if self.square==False:
        fileName=self.fileName
        f=open(fileName,mode='rb')
        data=f.read()
        print "Starting to read hamiltonian"
        print "Block size", self.n1
        start=4


        for i in range(self.n1):
          for j in range(i+1):
      
            itp=i+1
            it=j+1
            ir=itp*(itp-1)/2 + it

            starti=(ir-1)*8
            endi=(ir)*8



            eivw=struct.unpack(">"+"ff",data[starti:endi])
            eiv=complex(eivw[0],eivw[1]) 
            if self.hermitian==True:
              self.ham[j,i]=eiv
            else:
              self.ham[j,i]=eiv.conjugate()
            self.ham[i,j]=eiv.conjugate()
          if self.hermitian==True:
            self.ham[i,i]=self.ham[i,i].real

        f.close()

      print "End read"


def cproduct(v1,v2,n):
     #angle=sum([v1[k].conjugate()*v2[k] for k in range(n)])
     angle=numpy.vdot(v1[:n],v2[:n])
     rho=abs(angle)
     phi=math.atan2(angle.imag, angle.real) #For old python versions... in newer ones there is a seprate function
     return [rho,phi]

def fproduct(v1,v2,n,n1):
     #angle=sum([v1[k].conjugate()*v2[k] for k in range(n)])
     angle=numpy.vdot(v1[:n],v2[:n])-numpy.vdot(v1[n:n1],v2[n:n1])
     rho=abs(angle)
     phi=math.atan2(angle.imag, angle.real) #For old python versions... in newer ones there is a seprate function
     return [rho,phi]

def hproduct(h,eig_f,n1,with_b=False):
	ham1=numpy.zeros(shape=(n1,n1),dtype=complex)




	for j in range(n1):
		tmp1=numpy.dot(h,eig_f.eig_pairs[j].vector)
		if with_b==True:
			n=n1/2
			for i in range(n1/2,n1):
				tmp1[i]=-tmp1[i]
		for i in range(n1):
			ham1[i,j]=numpy.vdot(eig_f.eig_pairs[i].vector,tmp1)

  
	for j in range(n1):
		ham1[j,j]=0

	for j in range(n1):
		h[j,j]=0

	print norm(ham1)/n1

	print norm(h)/n1

def build_full(ham_h,ham_c):
	n1=ham_h.n1
	print "Building full H"

	full_h=numpy.zeros(shape=(n1*2,n1*2),dtype=complex)

	for i in range(n1):
		for j in range(n1):
			full_h[i,j]=ham_h.ham[i,j]
			full_h[n1+i,n1+j]=-ham_h.ham[i,j].conjugate()
			full_h[n1+i,j]=-ham_c.ham[j,i]
			full_h[i,n1+j]=ham_c.ham[i,j].conjugate()

	print "End build"

	return full_h



class ov_file(my_file): #Overlap File
     def __init__(self, fileName, resonant=False, label="",full=False):
         my_file.__init__(self,fileName,resonant,label)
         self.ov=[]


     def read(self):
         fileName=self.fileName
         f=open(fileName,mode='rb')
         data=f.read()
         n1=struct.unpack(">i",data[4:8])[0]

         print "Starting to read eigenvalues and eigenvectors"

         if (self.resonant):
          n=n1
         else:
          n=n1/2
         self.setN(n,n1)

         print "Number of eigenvectors",n1

         #w=struct.unpack(">" + "ff"*n1 , data[(4+3*4):(4+3*4+2*4*(n1))])
         #eig=[w[i]*27.21 for i in range(2*n1)]
         start=4+1*4;
         self.ov=numpy.zeros(shape=(self.n1,self.n1),dtype=float)


         for i in range(n1):
          #if eig[2*i]>0 or self.full==True:
           starti=start+i*(n1+1)*8+8
           endi=start+(i+1)*(n1+1)*8
           eivw=struct.unpack(">"+"ff"*(n1),data[starti:endi])
           for kk in range(n1):
            #print kk
            #self.ov[i,kk]=0.0
            self.ov[i,kk]=math.sqrt(eivw[2*kk]**2+eivw[2*kk+1]**2)
            if i==kk:
              self.ov[i,kk]=math.sqrt((eivw[2*kk]-1)**2+eivw[2*kk+1]**2)
           
           
           #self.appendPair(eig[2*i],eiv)

         f.close()
         #for i in range(n1):
          #self.ov[i,i]=self.ov[i,i]-1.0
         #print self.ov
         print "End read"

class fa_file(my_file): #FA File
     def __init__(self, fileName, resonant=False, label="",full=False):
         my_file.__init__(self,fileName,resonant,label)
         self.fa=[]


     def read(self):
         fileName=self.fileName
         f=open(fileName,mode='rb')
         data=f.read()
         n1=struct.unpack(">i",data[4:8])[0]

         print "Starting to read eigenvalues and eigenvectors"

         if (self.resonant):
          n=n1
         else:
          n=n1/2
         self.setN(n,n1)

         print "Number of eigenvectors",n1

         #w=struct.unpack(">" + "ff"*n1 , data[(4+3*4):(4+3*4+2*4*(n1))])
         #eig=[w[i]*27.21 for i in range(2*n1)]
         start=4+1*4;
         self.fa=numpy.zeros(shape=(self.n1),dtype=complex)


         for i in range(1):
          #if eig[2*i]>0 or self.full==True:
           #print i 
           starti=start+i*(n1+1)*8+8
           endi=start+(i+1)*(n1+1)*8
           eivw=struct.unpack(">"+"ff"*(n1),data[starti:endi])
                      #print kk
            #self.ov[i,kk]=0.0
           self.fa=[complex(eivw[2*k],eivw[2*k+1]) for k in range(n1)]
            #if i==kk:
            
           
           
           #self.appendPair(eig[2*i],eiv)

         f.close()
         #for i in range(n1):
          #self.ov[i,i]=self.ov[i,i]-1.0
         #print self.ov
         print "End read"




class spectraSet:
  def __init__(self,label=""):
    self.spectras=[]
    self.label=label
    self.cnt=0
    self.spectra_sum=[]
    self.spectra_halfsum=[]
  
  def read(self,biPairs=[[1,1]],nameMask="%1d_%2d.dat",labelMask="snapshot %1d %02d",linesTaken=[1],readRe=True,skipRead=1):
    cnt=0
    for pair in biPairs:
      b=pair[0]
      i=pair[1]
      self.spectras.append(spectraFile(nameMask%(b,i),labelMask%(b,i)))
      self.spectras[cnt].read(lines_taken=linesTaken,readRe=readRe,skip_read=skipRead)
      cnt=cnt+1
    self.cnt=cnt

  
  def calcHalfSum(self):
    spectra_sum=numpy.zeros(self.spectras[0].cnt)
    cnt_half=int(numpy.floor(self.cnt/2))
    for i in range(cnt_half):
      for e in range(self.spectras[0].cnt):
        spectra_sum[e]=self.spectras[i].spectra[e]+spectra_sum[e]
    self.spectra_halfsum=spectra_sum/cnt_half

  def calcSum(self):
    spectra_sum=numpy.zeros(self.spectras[0].cnt)
    for i in range(self.cnt):
      for e in range(self.spectras[0].cnt):
        #print e
        spectra_sum[e]=self.spectras[i].spectra[e]+spectra_sum[e]
    self.spectra_sum=spectra_sum/self.cnt





class spectraFile: #Program params
  def __init__(self,fileName,label=""):
    self.label=label #Label for files
    self.fileName=fileName
        
    self.energies=[]
    self.spectra=[]
    self.spectra_re=[]
    self.broad_spectra=[]
    self.broad_spectra_re=[]
    self.cnt=0
    self.hasRe=False

  def read(self,skip_read=0,lines_taken=[1],de=0.0,coef=1.0,readRe=True):
    f=open(self.fileName, 'r')
    cnt=0
    energies=[]
    spectra=[]
    if readRe: spectra_re=[]

    els=lines_taken


    for line in f:
      cnt=cnt+1
      if cnt>=skip_read:
        data=line.split()
        energies.append(float(data[0])+de)

        spectra_sum=0
        if readRe: spectra_sum_re=0
        for el in els:
          spectra_sum+=coef*float(data[el])
          if readRe: spectra_sum_re+=coef*float(data[el-1])
        spectra_sum=spectra_sum/len(els)
        if readRe: spectra_sum_re=spectra_sum_re/len(els)

        spectra.append(spectra_sum)
        if readRe: spectra_re.append(spectra_sum_re)

    f.close()

    self.spectra=spectra
    self.energies=energies
    self.cnt=len(energies)
    #print cnt
    #print len(spectra)
    #print len(energies)
    if readRe: self.spectra_re=spectra_re
    if readRe: self.hasRe=True


