from math import *	
from scipy import *
import scipy
from pylab import *
from matplotlib import *
import numpy.fft as nft
import scipy.optimize as spo
#from matplotlib import pyplot as plt
import pylab as plt
from matplotlib import rc
import numpy


import string
import sys
#from matplotlib import *
#from pylab import *
import os
import random
import time
#
# gamma function lives here:
#import scipy.special
from scipy.special import gamma
#from scipy.optimize import leastsq
from matplotlib import axis as aa
from threading import Thread
#
#
import datetime as dtm
import calendar
import operator
import urllib
import MySQLdb

import rbIntervals as rbi

# maping bits:
import matplotlib	# note that we've tome from ... import *. we should probably eventually get rid of that and use the matplotlib namespace.
#matplotlib.use('Agg')
#from matplotlib.toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

################################################
# Utilities
################################################

class linefit:
	# a simple line-fit class. we'll need X,Y,wt,
	# call like: lf1=m.linefit(dataSet)
	#
	#
	datas=[]		# [X, Y, Wt]
	totalVar=0			# sum( (x_i - f(x_i)**2 )
	errs=[]
	Ndof=0
	a=None
	b=None
	activeRes=None
	AIC=None
	nPrams=None	# len(fitData)=Ndof+nPrams
	
	#
	def __init__(self, inData=[]):
		self.initialize(inData)
	
	def initialize(self, inData=[]):
		# what does the data look like?
		self.activeRes=self.linRes
		dLen=len(inData)		# if dLen==2,3, then we probably have [[X], [Y], [wt]]
		#
		if dLen<=3:
			# we should probably sort the arrays...
			self.datas=inData
		if dLen>3:
			# data is probably like [[x,y,w], [x,y,w]...]
			if len(inData[0])==1:
				# 1D array; assume a plain ole sequence for X
				self.datas=[]
				self.datas+=[range(len(inData))]
				self.datas+=[map(operator.itemgetter(0), inData)]
			if len(inData[0])>=2:
				# assume the data are ordered pairs, so we can sort them on the x coordinate.
				inData.sort(key=operator.itemgetter(0))
				self.datas=[]
				self.datas+=[map(operator.itemgetter(0), inData)]
				self.datas+=[map(operator.itemgetter(1), inData)]
			if len(inData[0])>=3: self.datas+=[map(operator.itemgetter(2), inData)]
		if len(self.datas)==2:
			# add even weight.
			self.datas+=[[]]
			for x in self.datas[0]:
				self.datas[2]+=[1]
		#
	#
	def meanVar(self):
		# aka, rChiSqr
		return self.totalVar/self.Ndof
		
	#def doOmoriFit(self, p0=[None, None, None], xmin=None, xmax=None, r0=0):
	def doOmoriFit(self, p0=[None, None, None], xmin=None, xmax=None, r0=0, fitres=None):
		# the pram r0 is for the OmoriIntRes2() residual, which incorporates backtround seismicity.
		# do a linear fit:
		# a,b parameters are first guesses.
		# if they are None, guess starting prams:
		if fitres==None: fitres=self.omoriIntRes
		if p0==None: p0=[None, None, None]
		if p0[1]==None:
			p0[1]=(float(self.datas[1][-1])-float(self.datas[1][0]))/(float(self.datas[0][-1])-float(self.datas[0][0]))
		if p0[0]==None:
			p0[0]=self.datas[1][0]-p0[1]*self.datas[0][0]
		if p0[2]==None: p0[2]=1.0
		#self.a=a
		#self.b=b
		p=scipy.array(p0)
		#
		if xmin==None: xmin=self.datas[0][0]
		if xmax==None: xmax=self.datas[0][-1]
		# 
		# get X,Y,W:
		X=[]
		Y=[]
		W=[]
		for i in xrange(len(self.datas[0])):
			if self.datas[0][i]<xmin: continue
			#
			X+=[self.datas[0][i]]
			Y+=[self.datas[1][i]]
			W+=[self.datas[2][i]]
			if X[-1]>=xmax: break

		#
		# now, fit data:
		
		#print "do the fit..."
		# note: args are (y, x, wt)
		#plsq=spo.leastsq(self.linRes, p, args=(scipy.array(self.datas[1]), scipy.array(self.datas[0]), scipy.array(self.datas[2])), full_output=1)
		print "prams: %s" % str(p)
		#plsq=spo.leastsq(fitres, p, args=(scipy.array(Y), scipy.array(X), scipy.array(W), r0), full_output=1)
		plsq=spo.leastsq(fitres, p, args=(scipy.array(Y), scipy.array(X), scipy.array(W)), full_output=1)
		#print "fit done. sum error..."
		for sig in self.errs:
			self.totalVar+=sig*sig
		#self.Ndof=len(self.datas[0])-len(p)
		self.Ndof=len(X)-len(p)
		self.nPrams=len(p)
		self.AIC=self.totalVar+2*nPrams
		self.a=plsq[0][0]
		self.b=plsq[0][1]
		'''
		plsq=spo.leastsq(linRes, p, args=(ymax,x), full_output=0, maxfev=200000)
		amax=plsq[0][0]
		bmax=plsq[0][1]
		'''
		return plsq
	
	def doLogFit(self, lbase=10.0, a=None, b=None, xmin=None, xmax=None, thisdatas=None):
		# fit the log-log representation.
		if thisdatas==None: thisdatas=self.datas
		#print "datalen: %d" % len(thisdatas)
		logdatas=[[], [], []]
		#
		# get logarithms of data:
		for i in xrange(len(thisdatas[0])):
			logdatas[0]+=[math.log(thisdatas[0][i], lbase)]
			logdatas[1]+=[math.log(thisdatas[1][i], lbase)]
			wt=1
			if len(thisdatas)>=3: wt=thisdatas[2][i]
			logdatas[2]+=[wt]				
		#
		#return logdatas
		
		return self.doFit(a, b, xmin, xmax, logdatas)
		
	def doFit(self, a=None, b=None, xmin=None, xmax=None, thisdatas=None, fop=1):
		if thisdatas==None: thisdatas=self.datas
		# do a linear fit:
		# a,b parameters are first guesses.
		# if they are None, guess starting prams:
		self.errs=[]
		self.totalVar=0
		if b==None:
			b=(float(thisdatas[1][-1])-float(thisdatas[1][0]))/(float(thisdatas[0][-1])-float(thisdatas[0][0]))
		if a==None:
			a=thisdatas[1][0]-b*thisdatas[0][0]
		self.a=a
		self.b=b
		p=scipy.array([a,b])
		
		if xmin==None: xmin=thisdatas[0][0]
		if xmax==None: xmax=thisdatas[0][-1]
		# 
		# get X,Y,W:
		X=[]
		Y=[]
		W=[]
		for i in xrange(len(thisdatas[0])):
			if thisdatas[0][i]<xmin: continue
			#
			X+=[thisdatas[0][i]]
			Y+=[thisdatas[1][i]]
			W+=[thisdatas[2][i]]
			if X[-1]>=xmax: break
		#
		# now, fit data:
		
		#print "do the fit..."
		# note: args are (y, x, wt)
		#plsq=spo.leastsq(self.linRes, p, args=(scipy.array(self.datas[1]), scipy.array(self.datas[0]), scipy.array(self.datas[2])), full_output=1)
		plsq=spo.leastsq(self.activeRes, p, args=(scipy.array(Y), scipy.array(X), scipy.array(W)), full_output=fop)
		#print "fit done. sum error..."
		for sig in self.errs:
			self.totalVar+=sig*sig
		#self.Ndof=len(thisdatas[0])-len(p)
		self.AIC=self.totalVar+2.*len(p)
		self.Ndof=len(X)-len(p)
		self.dataLen=len(X)
		self.a=plsq[0][0]
		self.b=plsq[0][1]
		'''
		plsq=spo.leastsq(linRes, p, args=(ymax,x), full_output=0, maxfev=200000)
		amax=plsq[0][0]
		bmax=plsq[0][1]
		'''
		#
		# now, get error for a,b.
		# start by calculating Delta and DeltaPrime (error-partition for w!=1, w=1):
		self.delta=0
		self.deltaPrime=0
		aX=scipy.array(X)
		aXX=aX**2
		aY=scipy.array(Y)
		aYY=aY**2
		aW=scipy.array(W)
		aWW=aW**2
		
		#sigsSq=scipy.array(self.errs)**2	# i think this does not get squared (we're using linRes() )... it would not be a bad idea to just calc them here...
		#self.delta=sum
		# delta first (this might not be quite right except when w=1/sig^2, so be careful...)
		# note: we assume W=1/sigma**2
		#self.delta=sum(aW)*sum(aXX*aW) - (sum(aX*aW))**2	# check this formulation for order of operations...
		self.delta=sum(aWW)*sum(aXX*aWW) - (sum(aX*aWW))**2.	# check this formulation for order of operations...
		self.deltaPrime=float(len(X))*sum(aXX) - sum(aX)**2.
		#
		# weighted pram errors (variance):
		#thisSig=self.totalVar**.5
		#
		self.vara=(1.0/self.delta)*sum(aXX*aWW)
		self.varb=(1.0/self.delta)*sum(aWW)	# note that w=1/sig^2 in most texts, as per gaussian stats. this is more general than that, and it might make a big fat
		#										# mess when applied outside gaussian distributions.
		# w_i=1 case (note, this can be generallized to w_i=w; 1/delta -> 1/(w*delta):
		self.varaprime=(self.meanVar()/self.deltaPrime)*sum(aXX)
		self.varbprime=float(len(X))*self.meanVar()/self.deltaPrime
			
		
		return plsq		

	def tofile(self, fname='data/lfdata.dat', lfheader='#data from linefit object\n'):
		fout=open(fname, 'w')
		fout.write(lfheader)
		for i in xrange(len(self.datas[0])):
			fout.write('%f\t%f\t%f\n' % (self.datas[0][i], self.datas[1][i], self.datas[2][i]))
		fout.close()
	
	def fLin(self, x,p):
		return p[0] + x*p[1]
		
	def fPL(self, x,p):
		return (10**p[0])*(x**p[1])
	
		
	def linRes(self, p,y,x,w):
		err=y-(p[0]+x*p[1])
		self.errs=w*err
		#return w*err
		return self.errs
	
	def omoriRateRes(self, p, y, x, w):
		err=y-(1.0/(p[0]*(1+(x/p[1]))**p[2]))
		self.errs=w*err
		#return w*err
		return self.errs
	
	def omoriIntRes(self, p, y, x, w):
		err=y-(p[0]*(1+(x/p[1]))**p[2])
		werr=w*err
		self.errs=w*err
		#return werr
		return self.errs
	
	def omoriIntRes2(self, p, y, x, w, r0):
		# this is an omori like function that includes a "background rate" r0.
		#err=w*(y-(p[0]*(1+(x/p[1]))**p[2]))
		#r0=0.2222
		#err=w*(y- (1/(r0 + 1/(p[0]*(1+(x/p[1]))**p[2]) )) )
		err=w*(y- 1/(r0 + 1/(p[0]*(1+(x/p[1]))**p[2])) )
		self.errs=err
		return err
	
	def getFitPlotAry(self):
		print "from getFitPlotAry(): %s, %s, %s, %s" % (self.datas[0][0], self.datas[0][-1], self.a, self.b)
		Fx=[self.datas[0][0], self.datas[0][-1]]
		Fy=[self.datas[0][0]*self.b + self.a, self.datas[0][-1]*self.b + self.a]
		return [Fx, Fy]
		
	def quickPlot(self, toFig=True, colorNum=0):
		fitXY=self.getFitPlotAry()
		# toFig: do a figure here or assume that there is an active figure in the ether. this way, we can put many quickPlots onto a single canvas.
		
		if toFig:
			plt.figure(0)
			plt.clf()
		plt.plot(self.datas[0], self.datas[1], color=pycolor(colorNum))
		plt.plot(fitXY[0], fitXY[1], color=pycolor(colorNum), label='a=%f, b=%f, rChi^2=%f' % (self.a, self.b, self.meanVar()**.5))
		if toFig:
			plt.legend(loc='upper right')
			plt.show()

def getLogs(data=[], lbase=10):
	#output=[]
	#for rw in data:
	#	newrow=[]
	#	for elem in rw:
	#		newrow+=[math.log(elem, lbase)]
	#	output+=[newrow]
	output=[]
	for elem in data:
		if elem==0:
			output+=[None]
			continue
		output+=[math.log(elem, lbase)]
	return output
		
def pycolor(num=0):
	if num==None: num=0
	num=int(num%7)
	clrs=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
	
	return clrs[num]
def pyicon(num=0):
	if num==None: num=0
	num=int(num%20)
	# note: skip ',', the "pixel" marker
	icns=['.', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|', '_']
	return icns[num]
	
def getMonthName(monthNum=1):
	if monthNum==1:
		return "January"
	if monthNum==2:
		return "February"
	if monthNum==3:
		return 'March'
	if monthNum==4:
		return 'April'
	if monthNum==5:
		return 'May'
	if monthNum==6:
		return 'June'
	if monthNum==7:
		return 'July'
	if monthNum==8:
		return 'August'
	if monthNum==9:
		return 'September'
	if monthNum==10:
		return 'October'
	if monthNum==11:
		return 'November'
	if monthNum==12:
		return 'December'

def plotPolygons(polygons=None):
	import ygmapbits as yg
	# getPIsquareRays
	if polygons==None: polygons=yg.getReducedPolys()
	#if polygons==None: polygons=getPIsquareRays()
	#
	plt.figure(0)
	for ply in polygons:
		#if len(ply)<5: continue
		plt.fill(map(operator.itemgetter(1), ply), map(operator.itemgetter(0),ply), '.-')
	
	plt.show()
	
	return polygons
	
def printPolyLens(polygons=None):
	import ygmapbits as yg
	if polygons==None: polygons=yg.getReducedPolys()
	#
	i=0
	for ply in polygons:
		print "poly(%d): %d" % (i, len(ply))
		i+=1
	return None

##################
# write an all-python method to fetch ANSS data and insert into MySQL
# this will eventually replace {NEICANSS2sql.py}.fetchAndInsertANSS()
def getANSStoFilehandler(lon=[-125, -115], lat=[32, 45], minMag=4.92, dates0=[dtm.date(2001,01,01), dtm.date(2010, 12, 31)], Nmax=999999):
	# fetch data from ANSS; return a file handler.
	#
	# use urllib in "post" mode. an example from http://www.python.org/doc/current/library/urllib.html#urllib.FancyURLopener)
	# using "get" (aka, query-string method; note the ?%s string at the end of the URL, this is a single pram call to .urlopen):
	#
	#>>> import urllib
	#>>> params = urllib.urlencode({'spam': 1, 'eggs': 2, 'bacon': 0})
	#>>> f = urllib.urlopen("http://www.musi-cal.com/cgi-bin/query?%s" % params)
	#>>> print f.read()
	#
	# using "post" (note this is a 2 pram call):
	#>>> import urllib
	#>>> params = urllib.urlencode({'spam': 1, 'eggs': 2, 'bacon': 0})
	#>>> f = urllib.urlopen("http://www.musi-cal.com/cgi-bin/query", params)
	#>>> print f.read()
	#
	# make ANSS prams dictionary (thank james for the bash-template):
	# ANSSquery has day-resolution:
	dates=[dtm.date(dates0[0].year, dates0[0].month, dates0[0].day), dtm.date(dates0[1].year, dates0[1].month, dates0[1].day)]
	anssPrams={'format':'cnss', 'output':'readable', 'mintime':str(dates[0]).replace('-', '/'), 'maxtime':str(dates[1]).replace('-', '/'), 'minmag':str(minMag), 'minlat':lat[0], 'maxlat':lat[1], 'minlon':lon[0], 'maxlon':lon[1], 'etype':'E', 'searchlimit':Nmax}
	f = urllib.urlopen('http://www.ncedc.org/cgi-bin/catalog-search2.pl', urllib.urlencode(anssPrams))
	#
	# we might return f, a string of f, or maybe a list of lines from f. we'll work that out shortly...
	return f

def getANSSlist(lon=[-125, -115], lat=[32, 45], minMag=4.92, dates0=[dtm.date(2001,01,01), dtm.date(2010, 12, 31)], Nmax=999999, fin=None):
	#
	# note: this appears to be a bad idea for global downloads. a full catalog is ~4GB, which kills my computer.
	#
	# note: this may be repeated exactly in ygmapbits.py
	# fetch new ANSS data; return a python list object of the data.
	# fin: data file handler. if this is None, then get one from ANSS.
	dates=[dtm.date(dates0[0].year, dates0[0].month, dates0[0].day), dtm.date(dates0[1].year, dates0[1].month, dates0[1].day)]
	anssList=[]
	if fin==None:
		#print "get data from ANSS...(%s, %s, %s, %s, %s)" % (lon, lat, minMag, dates, Nmax)
		fin = getANSStoFilehandler(lon, lat, minMag, dates, Nmax)
		#fin = getANSStoFilehandler([-180, 180], [-90, 90], 0, [datetime.date(1910,01,01), datetime.date(2010, 01, 16)], 9999999)

		print "data handle fetched..."
		
	for rw in fin:
		if rw[0] in ["#", "<"] or rw[0:4] in ["Date", "date", "DATE", "----"]:
			#print "skip a row... %s " % rw[0:10]
			continue
		#anssList+=[rw[:-1]]
		# data are fixed width delimited
		# return date-time, lat, lon, depth, mag, magType, nst, gap, clo, rms, src, catEventID (because those are all the available bits)
		#print "skip a row... %s " % rw
		rwEvdt=rw[0:22].strip()
		rwLat=rw[23:31].strip()
		if rwLat=='' or isnumeric(str(rwLat))==False or rwLat==None:
			continue
			#rwLat=0.0
		else:
			rwLat=float(rwLat)
		rwLon=rw[32:41].strip()
		if rwLon=='' or isnumeric(str(rwLon))==False or rwLon==None:
			#rwLon=0.0
			continue
		else:
			rwLon=float(rwLon)
		rwDepth=rw[42:48].strip()
		if rwDepth=='' or isnumeric(str(rwDepth))==False or rwDepth==None or str(rwDepth).upper() in ['NONE', 'NULL']:
			#rwDepth=0.0
			rwDepth=None
			continue
		else:
			rwDepth=float(rwDepth)
		rwMag=rw[49:54].strip()
		if rwMag=='' or isnumeric(str(rwMag))==False or rwMag==None:
			#rwMag=0.0
			continue
		else:
			rwMag=float(rwMag)
		rwMagType=rw[55:59].strip()
		rwNst=rw[60:64].strip()
		if rwNst=='':
			rwNst=0.0
		else:
			rwNst=float(rwNst)
		rwGap=rw[65:68].strip()
		rwClo=rw[69:73].strip()
		rwrms=rw[74:78].strip()
		if rwrms=='':
			rwrms=0.0
		else:
			rwrms=float(rwrms)		
		rwsrc=rw[79:83].strip()
		rwCatEventId=rw[84:96].strip()
		
		#anssList+=[[rw[0:22].strip(), float(rw[23:31].strip()), float(rw[32:41].strip()), float(rw[42:48].strip()), float(rw[49:54].strip()), rw[55:59].strip(), float(rw[60:64].strip()), rw[65:68].strip(), rw[69:73].strip(), float(rw[74:78].strip()), rw[79:83].strip(), rw[84:96].strip()]]
		anssList+=[[rwEvdt, rwLat, rwLon, rwDepth, rwMag, rwMagType, rwNst, rwGap, rwClo, rwrms, rwsrc, rwCatEventId]]
	return anssList
#
def replaceANSS2SQL(anssList=None, catID=523, tempCatID=None):
	# just a familiar-name wrapper:
	return ANSSlist2SQL(anssList, catID, tempCatID)
	
def ANSSlist2SQL(anssList=None, catID=523, tempCatID=None):
	# this is the end-product; call it and update the WHOLE CATALOG. BUT, this replaces the whole catalog, so it takes forever.
	# we need an update system...
	#
	# insert an ANSSlist into MySQL as catalogID={catID}. use tempCatID to temporarily dump existing data to a safe place,
	# should something go wrong. reserve an option to skip this step for brevity.
	#
	sqlHost = 'localhost'
	sqlUser = 'myoder'
	sqlPassword = 'yoda'
	sqlPort = 3306
	sqlDB = 'QuakeData'
	myConn = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	myConn2 = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	strDel=""
	# if no list has been provided, get a current, complete world catalog:
	if anssList==None: anssList=getANSSlist([-180, 180], [-90, 90], 0, [dtm.date(1932,01,01), dtm.date.fromordinal(dtm.datetime.now().toordinal())], 9999999)
	print "ANSS list retrieved, len: %d" % len(anssList)

	if tempCatID==None:
		# get a save catID; use that one.
		#strGetMaxID="select min(catalogID) from Earthquakes where catalogID>%d" % catID
		strGetMaxID="select max(catalogID) from Earthquakes"
		myConn.query(strGetMaxID)
		rID=myConn.store_result()
		# there should only be one value. for now, assume this is true (which is, in general, really bad DB programming form):
		tempCatID=int(rID.fetch_row()[0][0])+1
		#rID=None
		strDel='delete from Earthquakes where catalogID=%d; update Earthquakes set catalogID=%d where catalogID=%d' % (tempCatID, tempCatID, catID)
		
	elif tempCatID==0:
		# skip this step; just delete the catalog and hope all goes well.
		strDel='delete from Earthquakes where catalogID=%d' % catID		
	#
	else:
		# do we have a real backup catalogID? update the "real" id with the tempCatID value.
		strDel='update Earthquakes set catalogID=%d where catalogID=%d' % (tempCatID, catID)
		
	myConn.query(strDel)
	
	#
	print "begin sql insert loop..."
	for rw in anssList:
		# write an insert string or command or whatever...
		# one string, or a bunch? i don't think it matters for inserts.
		# ... and this whole thing hast to be error handled (None, Null, quotes, etc.)...
		#print "rw: %s" % rw
		#insStr = "insert into QuakeData.Earthquakes (catalogID, eventDateTime, lat, lon, depth, mag, magType, nst, gap, clo, rms, src, catEventID) values (%d, '%s',%s, %s,%s, %s,%s, %s,%s, '%s',%s, '%s', '%s') " % (catID, str(rw[0]), str(rw[1]), str(rw[2]), str(rw[3]), str(rw[4]), str(rw[5]), str(rw[6]), str(rw[7]), str(rw[8]), str(rw[9]), str(rw[10]), str(rw[11]) )
		insStr = "insert into QuakeData.Earthquakes (catalogID, eventDateTime, lat, lon, depth, mag, magType, src, catEventID) values (%d, '%s', %f, %f, %f, %f, '%s', '%s', '%s')" % (catID, str(rw[0]), float(rw[1]), float(rw[2]), float(rw[3]), float(rw[4]), str(rw[5]), str(rw[10]), str(rw[11]) )
		#print insStr
		# it might be faster to bundle these statements. also, it might not be a bad idea to remove unique constraints from the Earthquakes table to facilitate faster
		# inserts. nominally, actual inserts could be done on threads to make it super speedy.
		myConn2.query(insStr)
	
	myConn.close()
	myConn2.close()


def updateANSS2SQL(catID=523):
	# this is the end-product. call this to update MySQL.QuakeData.Earthquakes with the most recent events; it will call other functions as necessary.
	# this function assumes all existing data are correct and complete.
	# here, we get the date of the most recent event and select from ANSS all events after that.
	#
	# also, ANSS allows only day-level queries, so we have to delete everything date>date0, then replace with all events date>date0
	#
	sqlHost = 'localhost'
	sqlUser = 'myoder'
	sqlPassword = 'yoda'
	sqlPort = 3306
	sqlDB = 'QuakeData'
	myConn = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	myConn2 = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	strDel=""
	# if no list has been provided, get a current, complete world catalog:
	#if anssList==None: anssList=getANSSlist([-180, 180], [-90, 90], 0, [dtm.date(1932,01,01), datetime.date.fromordinal(datetime.datetime.now().toordinal())], 9999999)
	anssList=[]
	# get a maximum date:
	#sqlMaxDate="select max(eventDateTime) from Earthquakes where catalogID=%d" % catID
	maxdates=myConn.cursor()
	maxdates.execute("select max(eventDateTime) from Earthquakes where catalogID=%d" % catID)
	maxDt=maxdates.fetchone()[0]
	print "getList prams: (%s, %s, %s, %s, %s)" % ([-180, 180], [-90, 90], 0, [maxDt, dtm.date.fromordinal(dtm.datetime.now().toordinal())], 9999999)
	anssList=getANSSlist([-180, 180], [-90, 90], 0, [maxDt, dtm.date.fromordinal(dtm.datetime.now().toordinal())], 9999999)
	
	# now, delete most recent days events (they will be replaced):
	strDel = "delete from Earthquakes where eventDateTime>='%s' and catalogID=%d" % (str(dtm.date(maxDt.year, maxDt.month, maxDt.day)), catID)
	myConn.query(strDel)
	for rw in anssList:
		# write an insert string or command or whatever...
		# one string, or a bunch? i don't think it matters for inserts.
		# ... and this whole thing hast to be error handled (None, Null, quotes, etc.)...
		#print "rw: %s" % rw
		#insStr = "insert into QuakeData.Earthquakes (catalogID, eventDateTime, lat, lon, depth, mag, magType, nst, gap, clo, rms, src, catEventID) values (%d, '%s',%s, %s,%s, %s,%s, %s,%s, '%s',%s, '%s', '%s') " % (catID, str(rw[0]), str(rw[1]), str(rw[2]), str(rw[3]), str(rw[4]), str(rw[5]), str(rw[6]), str(rw[7]), str(rw[8]), str(rw[9]), str(rw[10]), str(rw[11]) )
		insStr = "insert into QuakeData.Earthquakes (catalogID, eventDateTime, lat, lon, depth, mag, magType, src, catEventID) values (%d, '%s', %f, %f, %f, %f, '%s', '%s', '%s')" % (catID, str(rw[0]), float(rw[1]), float(rw[2]), float(rw[3]), float(rw[4]), str(rw[5]), str(rw[10]), str(rw[11]) )
		#print insStr
		# it might be faster to bundle these statements. also, it might not be a bad idea to remove unique constraints from the Earthquakes table to facilitate faster
		# inserts. nominally, actual inserts could be done on threads to make it super speedy.
		myConn2.query(insStr)
	#
	myConn.close()
	myConn2.close()
	
	return anssList

###############

def isnumeric(value):
  return str(value).replace(".", "").replace("-", "").isdigit()

def getValsAbove(inList, aboveVal=0):
	# input a 1D list. return: the value if it's above aboveVal, otherwise aboveVal:
	retList=[]
	for x in inList:
		if x>=aboveVal:
			retList+=[x]
		else:
			retList+=[aboveVal]
	return retList

def getValsBelow(inList, belowVal=0):
	# input a 1D list. return: the value if it's above aboveVal, otherwise aboveVal:
	retList=[]
	for x in inList:
		if x<=belowVal:
			retList+=[x]
		else:
			retList+=[belowVal]
	return retList

def frange(xmax=10.0, x0=0.0, dx=1.0):
	if dx==0: dx=1.0
	if (xmax<x0 and dx>0) or (xmax>x0 and dx<0): dx=-dx
	X=[float(x0)]
	while (X[-1]<xmax and xmax>x0) or (X[-1]>xmax and xmax<x0):
		X+=[float(X[-1]+dx)]
	return X

def datetimeToFloat(dtIn=dtm.datetime.now()):
	# returns a float version in units of days.
	fdt= dtIn.toordinal() + float(dtIn.hour)/24.0 + dtIn.minute/(24.0*60.0) + dtIn.second/(24.0*3600.0) + dtIn.microsecond/(24.0*3600000000.0)
	return fdt

def timeDeltaToFloat(tdelta=None, timeunit='day'):
	if tdelta==None: return None	# or maybe some default val.
	# time deltas are (days, secs., microsecs.)
	# let's, by default, convert to days. when we add other time units, we'll just multiply or divide accordingly.
	#
	if timeunit.lower() in ['day', 'days', 'dy', 'dys']:
		#
		rval=tdelta.days + (tdelta.seconds + (float(tdelta.microseconds)/10**6)/(3600*24.0) )
	else:
		rval=timeDeltaToFloat(tdelta, 'day')
	#
	return rval

#def datetimeFromFloat(fin=1000.0):
#	dt=dtm.datetime.fromordinal(int(fin))
#	dayfrac=fin%1
#	hrs=int(round(dayfrac*24))
#	mins=int(round(dayfrac*24*60/24))
#	secs=int(round(dayfrac*24*3600/(24*60)))
#	msecs=int(round(dayfrac*24*3600%1)*10**6)
#	print hrs, mins, secs, msecs
#	#
#	return dtm.datetime(dt.year, dt.month, dt.day, hrs, mins, secs, msecs)

def deg2rad(theta):
	return 2.0*pi*theta/360.0
	
def ellipseY(x, a, b):
	#print b, x, a
	return b*(1.0-x*x/(a*a))**.5
	
def rotatexy(x, y, Lat, Lon, theta):
	# x,y to transform via blah, blah.
	#
	theta=deg2rad(float(theta))
	xprime = (x-Lon)*cos(theta) - (y-Lat)*sin(theta)
	yprime = (x-Lon)*sin(theta) + (y-Lat)*cos(theta)
	
	return [xprime, yprime]

def datetimeFromString(strDtin=dtm.datetime.now().isoformat(' '), delim='-'):
	possibleDelims=['/', '-']
	if delim not in strDtin:
		for dlm in possibleDelims:
			if dlm in strDtin: delim=dlm
			
	#print strDtin
	strDt=strDtin.split(' ')[0]
	#print strDt
	if ' ' in strDtin and ':' in strDtin:
		strTm=strDtin.split(' ')[1]
	else:
		strTm="0:0:0.00"
	if '.' not in strTm: strTm=strTm+'.00'
	#
	dts=strDt.split(delim)
	tms=strTm.split(':')
	
	secFract='.'+tms[2].split('.')[1]
	if secFract=='':
		secFract=0.0
	else:
		secFract=float(secFract)
	
	yr=int(dts[0])
	mnth=int(dts[1])
	dy=int(dts[2])
	#
	hrs=int(tms[0])
	mns=int(tms[1])
	secs=tms[2].split('.')[0]
	if secs=='': secs=0
	secs=int(secs)
	msecs=int(secFract*10**6)
	#
	# and we keep seeing seconds>60...
	if secs>=60:
		secs=secs%60
		mns+=secs/60
	if mns>=60:
		mns=mns%60
		hrs+=mns/60
	retDate=dtm.datetime(yr, mnth, dy, hrs, mns, secs, msecs)
	if hrs>=24:
		hrs=hrs%24
		retDate+=dtm.timedelta(days=hrs/24)
	#
	
	
	#print yr, mnth, dy, hrs, mns, secs, msecs
	
	#return None
	#return dtm.datetime(yr, mnth, dy, hrs, mns, secs, msecs)
	return retDate
	
def averageOver(inData=[], n=1):
	# return average over n elements.
	# return 1:1 rows, note that the first n elements will be averaged over n`<n
	#import numpy
	outData=[[],[]]	# <x>, stdev
	#
	inlogs=getLogs(inData)
	N=1
	currVal=0
	for i in xrange(len(inData)):
		#outData[0]+=[sum(inData[i-N:i+1])/float(N)]
		#
		#outData[0]+=[numpy.mean(inData[i+1-N:i+1])]
		#outData[1]+=[numpy.std(inData[i+1-N:i+1])]
		#
		#outData[0]+=[(numpy.prod(inData[i+1-N:i+1]))**(1./N)]		# this is correct, but to get stdev (which has exactly what meaning in this case?), we need to get the logs anyway.
		outData[0]+=[numpy.mean(inlogs[i+1-N:i+1])]
		outData[1]+=[numpy.std(inlogs[i+1-N:i+1])]
		
		if N<n: N+=1
	return outData

def logaverageOver(inData=[], n=1):
	# return average over n elements.
	# return 1:1 rows, note that the first n elements will be averaged over n`<n
	#import numpy
	outData=[[],[]]	# <x>, stdev
	#
	N=1
	currVal=0
	for i in xrange(len(inData)):
		#outData[0]+=[sum(inData[i-N:i+1])/float(N)]
		outData[0]+=[numpy.mean(inData[i+1-N:i+1])]
		outData[1]+=[numpy.std(inData[i+1-N:i+1])]
		
		if N<n: N+=1
	return outData

def greaterof(a,b):
	if a>=b: return a
	if b>a: return b

def lesserof(a,b):
	if a<=b: return a
	if b<a: return b	

def vlinePadList(lst, minVal=0):
	# pad list vals so they plot as vertical spikes.
	newLst=[]
	for rw in lst:
		# assume: [[x,y]...]
		if rw[1]<minVal: continue
		newLst+=[[rw[0], minVal], [rw[0], rw[1]], [rw[0], minVal]]
	#
	return newLst

def getIntervals(catList, winLen):
	catLen=len(catList)
	i=(catLen-1-winLen)	# start winLen positions from the end.
	thisInterval=0
	#N=1
	intervals=[]	# [[eventDateTime, totalInterval]]
	while i>=0:
		#
		thisInterval=datetimeToFloat(catList[i+winLen][0])-datetimeToFloat(catList[i][0])
		intervals+=[[catList[i+winLen][0], thisInterval]]
		i-=1
		
	#
	#return [intervals, catList]
	return intervals
	
class eqcatalog:
	# a simple catalog class. we'll need X,Y,wt,
	# call like: lf1=m.catalog()
	#
	#
	mc=None	# catalog threshold
	cat=[]	# let's use [[evDateTime, lat, lon, mag, a, b], [...]] and we'll plot with map(), etc.
	subcats=[]
	activefig=None
	catmap=None
	#
	#
	def __init__(self, inData=[]):
		#self.cat=[]
		#self.subcats=[]
		self.initialize(inData)
	
	def initialize(self, inData=[]):
		# what does the data look like?
		self.cat=[]
		self.subcats=[]
		#
		self.cat=inData
		self.sqlport=3306
		self.sqlhost='localhost'
		inData=None
		self.catmap=None
		self.__name__='eqcatalog'
		self.mapres='l'	# basemap map resolution. at some pont, add functions to sort out nonsensical values.
		#
	#
	def writeCatToFile(self, foutname=None, cat=None):
		if foutname==None: foutname='outcat.cat'
		if cat==None: cat=self.getcat(0)
		#
		fout=open(foutname, 'w')
		fout.write("#dtm, lat, lon, mag\n")
		for rw in cat:
			fout.write("%d/%d/%d %d:%d:%d.%d\t%f\t%f\t%f\n" % (rw[0].year, rw[0].month, rw[0].day, rw[0].hour, rw[0].minute, rw[0].second, rw[0].microsecond, rw[1], rw[2], rw[3]))
		fout.close()
			
	def loadCatFromFile(self, fname=None, minmag=None):
		# standard file format: date \s time \t lat \t lon \t mag
		#
		#print "fname: %s" % fname
		self.cat=[]
		fin=open(fname)
		for rw in fin:
			if rw[0] in ('#', '\t', '\n', ' '): continue
			#rws=rw.split('\t')
			rws=rw.split()
			#rint rws
			#if len(rws)<4: continue
			if len(rws)<5: continue
			if minmag!=None:
				#if float(rws[3])<minmag: continue
				if float(rws[4])<minmag: continue
			#self.cat+=[[datetimeFromString(rws[0], float(rws[1]), float(rws[2]), float(rws[3])]]
			#if '24:' in rws[0]: print rws[0]
			#print rws[1]
			#print rws
			self.cat+=[[datetimeFromString(rws[0] + ' '+ rws[1]), float(rws[2]), float(rws[3]), float(rws[4])]]
		fin.close()
	#
	def getMainEvent(self, thiscat=None):
		# return catalog row of max magnitude (epicenter location (more or less)) event. note, by default we use ths shock-cat because it will
		# be faster than using the fullCat AND, the fullCat is likely to have other large earthquakes.
		if thiscat==None: thiscat=self.cat
		maxMag=thiscat[0][3]
		maxIndex=0
		for i in xrange(len(thiscat)):
			#print i, maxMag, maxIndex, cat[i][3]
			if thiscat[i][3]>maxMag:
				maxIndex=i
				maxMag=thiscat[i][3]
		return thiscat[maxIndex] + [maxIndex]
	
	def getIndexDtm(self, mindt=None, cat=None, datecol=0):
		if cat==None or type(cat).__name__ not in ('list', 'tuple'): cat=self.cat
		if mindt==None or type(mindt).__name__!='datetime': mindt=self.getMainEvent()[0]
		#
		for rw in cat:
			if rw[datecol]>=mindt: return rw
		return None
		#
		return thiscat[maxIndex] + [maxIndex]
	
	def getSubCat(self, catindex=0):
		if len(self.subcats)==0: return None
		#
		return self.subcats[catindex][1]
	
	def getcat(self, catindex=0):
		# more general and simpler than getSubCat. 0 -> maincat, 1, etc. -> subcats:
		# it would probably be a good idea to also restructure how catalogs are stored inte class:
		# catalogs=[[maincat], [subcat1], [subcat2], ...]
		# self.cat -> catalogs[0]
		if catindex==0: return self.cat
		if len(self.subcats)<(catindex): return None	# this cat index does not exist
		return self.subcats[catindex-1][1]
	
	def getLatLonRange(self, cat=None, latloncols=[1,2]):
		if cat==None: cat=self.cat
		if latloncols==None: latloncols=[1,2]	# latitude, lon cols of catalog (order is lat, lon).
		#
		minLat=cat[0][latloncols[0]]
		maxLat=cat[0][latloncols[0]]
		minLon=cat[0][latloncols[1]]
		maxLon=cat[0][latloncols[1]]
		#
		for rw in cat:
			thisLat=rw[latloncols[0]]
			thisLon=rw[latloncols[1]]
			#
			if thisLat>maxLat: maxLat=thisLat
			if thisLat<minLat: minLat=thisLat
			if thisLon>maxLon: maxLon=thisLon
			if thisLon<minLon: minLon=thisLon
		#
		return [[minLat, minLon], [maxLat, maxLon]]
	
	# subcat (subcatalog) functions:
	def ellipseCat(self, fullcat=None, theta=0, clat=35.9, clon=-120.5, ra=1.0, rb=1.0):
		#
		
		#print "event (start) date, catname: %s, %s, %s" % (eventDate, catFname, self.catname)
		#
		if fullcat==None: fullcat=self.cat
		#self.subcats+=[[subcatname], []]
		tempcat=[]
		
		#nEventsSinceMS=0
		for row in fullcat:
			# rotate each element into our aftershock axis, is it in the ellipse?
			newVec=rotatexy(row[2], row[1], clat, clon, theta)
			#
			# is the rotated vector in our ellipse?
			if abs(newVec[0])>ra: continue
			Y=ellipseY(newVec[0], ra, rb)
			if abs(newVec[1])>Y: continue
			# dtm, lat, lon, mag, tX, tY 		(note this is like y,x, x`, y` for the space coordinates).
			#self.subcats[-1][1]+=[[row[0], row[1], row[2], row[3], newVec[0], newVec[1]]]
			tempcat+=[[row[0], row[1], row[2], row[3], newVec[0], newVec[1]]]
		return tempcat
	
	def polycat(self, cat=None, verts=None):
		# as per james' counsel, a "knot theory" approach is much simpler. independent of right/left handedness, the sum of above/below tests is >0 for
		# points inside, 0 for points outside (like EnM).
		# start by making verts -> vectors -> f(x)
		#
		# verts are like: [[x0,y0], [x1, y1], ..., [xn, yn], [x0, y0]]; last point is optional.
		if cat==None: cat=self.cat
		if verts==None or len(verts)<3:
			# don't know. if we don't have verts, what can we do?
			# also, we need at least 3 verts, or we just have a line.
			return None
		#
		if verts[-1]!=verts[0]: verts+=[verts[0]]
		#
		vecs=[]	# like [ [[x0,y0], [x1, y1]], [[x1, y1], [x2, y2]], ...]
		#vecdirs=[]	# vector directions; -1=left, 0=none, 1=right. this determines whether we want to be over or under. the first x-direction vector definds a right/left poly.
		# get lat, lon extrema and make vectors:
		extremeVerts=[verts[0][0], verts[0][0], verts[0][1], verts[0][1]]	# [minLon, maxLon, minLat, maxLat]
		for i in xrange(len(verts)-1):
			vecs+=[[verts[i], verts[i+1]]]
			if verts[i+1][0]>extremeVerts[1]: extremeVerts[1]=verts[i+1][0]
			if verts[i+1][0]<extremeVerts[0]: extremeVerts[0]=verts[i+1][0]
			if verts[i+1][1]>extremeVerts[3]: extremeVerts[3]=verts[i+1][1]
			if verts[i+1][1]<extremeVerts[2]: extremeVerts[2]=verts[i+1][1]
			#
			# and keep a list of vector directions (right,left; do we need up, down?):
			thisdir=0	# reserve for vertical elements.
			if verts[i+1][0]>verts[i][0]: thisdir=1 #CatMap
			if verts[i+1][0]<verts[i][0]: thisdir=-1
			#vecdirs+=[thisdir]
		#
		# we don't really need the center, but it might be useful later:
		center=scipy.array([extremeVerts[0] + (extremeVerts[1]-extremeVerts[0])/2.0, extremeVerts[2] + (extremeVerts[3]-extremeVerts[2])/2.0])
		#
		# and this way, we don't need the poly-direction.
		# now we can spin through the catalog. inout=sum(x^ * above/below). inout=0 means out; inout>0 means in.
		# where x^ is {-1, 1} for left, right; above/below is {-1, 1} for point is above, below. i don't think which one is -1 and which is 1 matters
		# so long as we are consistent. also, as per old-school gaussian integrals, the number of times we cross a boundary: odd-> in , even -> out
		# applies as well.
		polycat=[]
		for iev in xrange(len(cat)):
			event=cat[iev]
			x=event[2]
			y=event[1]
			# for speed, if we are outside the extreme vertices, move on:
			if (x<extremeVerts[0] or x>extremeVerts[1] or y<extremeVerts[2] or y>extremeVerts[3]):
					#print "extreme kill (%d, %d)" % (x, y)
					#keepEvent=0
					# and we're done...
					continue
			#
			#keepEvent=1	# start by assuming we keep the event.
			inPolyTracker=0	# running up/down score. by default, do not keep the event.
			#print "*#*#*#"
			for ivec in xrange(len(vecs)):
				vec=vecs[ivec]
				# make a line (if it's not vertical):
				if vec[1][0]-vec[0][0]==0: continue	# vertical segments do not contribute, and we'll get x/0 error.
				b=(vec[1][1]-vec[0][1])/(vec[1][0]-vec[0][0])
				a=vec[0][1]-b*vec[0][0]
				y0=a+b*x
				#
				# xrange:
				if vec[0][0]>vec[1][0]:
					bigX=vec[0][0]
					smallX=vec[1][0]
				if vec[0][0]<vec[1][0]:
					bigX=vec[1][0]
					smallX=vec[0][0]
				#
				# debug:
				#if iev<40:
				#	print vec[0][0], vec[1][0], x, lookUpDown, y, y0
				#	print (x>=smallX and x<=bigX), (lookUpDown==-1 and y>y0 ), (lookUpDown==1 and y<y0)
				#
				# are we in the current xrange?
				if (x<smallX or x>bigX): continue
				# if it's on the line, keep it:
				if y==y0:
					inPolyTracker=1
					continue
				# is it inside the polygon?
				if y>y0: isUp=1							# point is above
				if y<y0: isUp=-1							# point is below
				if vec[1][0]>vec[0][0]: vecDir=1		# to the right
				if vec[1][0]<vec[0][0]: vecDir=-1	# to the left
				inPolyTracker+=(vecDir*isUp)
				#
			#
			if inPolyTracker!=0: polycat+=[event]
			
		#print extremeVerts
		return polycat
		
		
	def polycat_cp(self, cat=None, verts=None):
		# my original version of polycat using cross products to determine the direction of the polygon. there is a faster way...
		# verts are like: [[x0,y0], [x1, y1], ..., [xn, yn], [x0, y0]]; last point is optional.
		if cat==None: cat=self.cat
		if verts==None or len(verts)<3:
			# don't know. if we don't have verts, what can we do?
			# also, we need at least 3 verts, or we just have a line.
			return None
		#
		if verts[-1]!=verts[0]: verts+=[verts[0]]
		#
		vecs=[]	# like [ [[x0,y0], [x1, y1]], [[x1, y1], [x2, y2]], ...]
		vecdirs=[]	# vector directions; -1=left, 0=none, 1=right. this determines whether we want to be over or under. the first x-direction vector definds a right/left poly.
		# get lat, lon extrema and make vectors:
		extremeVerts=[verts[0][0], verts[0][0], verts[0][1], verts[0][1]]	# [minLon, maxLon, minLat, maxLat]
		for i in xrange(len(verts)-1):
			vecs+=[[verts[i], verts[i+1]]]
			if verts[i+1][0]>extremeVerts[1]: extremeVerts[1]=verts[i+1][0]
			if verts[i+1][0]<extremeVerts[0]: extremeVerts[0]=verts[i+1][0]
			if verts[i+1][1]>extremeVerts[3]: extremeVerts[3]=verts[i+1][1]
			if verts[i+1][1]<extremeVerts[2]: extremeVerts[2]=verts[i+1][1]
			#
			# and keep a list of vector directions (right,left; do we need up, down?):
			thisdir=0	# reserve for vertical elements.
			if verts[i+1][0]>verts[i][0]: thisdir=1
			if verts[i+1][0]<verts[i][0]: thisdir=-1
			vecdirs+=[thisdir]
		#
		#print vecdirs
		# now, is the poly right or left handed? from the poly center, calculate the mean r x v.
		center=scipy.array([extremeVerts[0] + (extremeVerts[1]-extremeVerts[0])/2.0, extremeVerts[2] + (extremeVerts[3]-extremeVerts[2])/2.0])
		#print "verts: %s" % str(verts)
		#print "vecs: %s" % str(vecs)
		#print "center: %s" % str(center)
		polyDir=0
		for vec in vecs:
			# get the cross product r x vec:
			thisvec=scipy.array(vec[1])-scipy.array(vec[0])
			rvec=scipy.array(vec[0])-center
			cprod=numpy.cross(rvec, thisvec)
			#print vec, thisvec, rvec, cprod, type(cprod)
			# so i guess when your vectors are coplanar, scipy.array knows to retun just a scalar for the cross-product.
			polyDir+=cprod
		if polyDir>0: polyDir=1
		if polyDir<0: polyDir=-1
		#print "polyDir: %f" % polyDir
		#
		# now we can spin through the catalog to find elements above/below poly segments, depending on the direction of the segment and right/left
		# handedness of the poly.
		polycat=[]
		for iev in xrange(len(cat)):
			event=cat[iev]
			x=event[2]
			y=event[1]
			keepEvent=1	# start by assuming we keep the event.
			#print "*#*#*#"
			for ivec in xrange(len(vecs)):
				# test the event against each polygon segment. if it falls outside one or more, don't keep it...
				vec=vecs[ivec]
				# make a line:
				if vec[1][0]-vec[0][0]==0: continue
				#
				b=(vec[1][1]-vec[0][1])/(vec[1][0]-vec[0][0])
				a=vec[0][1]-b*vec[0][0]
				#
				lookUpDown=vecdirs[ivec]*polyDir
				y0=a+b*x
				#keep criteria:
				#if (x>=vec[0][0] and vec<=vec[1][0]) and ((lookUpdown==-1 and y<=y0 ) or (lookUpdown==1 and y>=0)):
				# so discard criteria is opposite (in y):
				if vec[0][0]>vec[1][0]:
					bigX=vec[0][0]
					smallX=vec[1][0]
				if vec[0][0]<vec[1][0]:
					bigX=vec[1][0]
					smallX=vec[0][0]
				#	
				if (x<extremeVerts[0] or x>extremeVerts[1] or y<extremeVerts[2] or y>extremeVerts[3]):
					print "extreme kill (%d, %d)" % (x, y)
					keepEvent=0
				if ((x>=smallX and x<=bigX) and ((lookUpDown==-1 and y>y0 ) or (lookUpDown==1 and y<y0))) :
					keepEvent=0
					print "f(x) kill (%d, %d)" % (x, y)
					# and for efficiency:
					continue
				#
			#
			if keepEvent==1: polycat+=[event]
		
			
		print extremeVerts
		return polycat
			
		
#[0,0], [2,0], [4,4], [2,6], [0,4]			
	
	def addEllipCat(self, subcatname='newcat', fullcat=None, theta=0, clat=35.9, clon=-120.5, ra=1.0, rb=1.0):
		#
		if fullcat==None: fullcat=self.cat
		
		newcat=self.ellipseCat(fullcat, theta, clat, clon, ra, rb)
		self.subcats+=[[subcatname, newcat]]
		
	def getMagSubcat(self, fullcat=None, minmag=2.5, magcol=3):
		newcat=[]
		for rw in fullcat:
			if rw[magcol]>=minmag: newcat+=[rw]
		return newcat
	def addMagSubcat(self, subcatname='magsubcat', fullcat=None, minmag=2.5, magcol=3):
		subcatname='%s-%s' % (subcatname, str(minmag))
		self.subcats+=[[subcatname, getMagSubcat(fullcat, minmag, magcol)]]
	#
	def getxytSubcat(self, fullcat=None, dts=[], lats=[], lons=[], llcols=[1,2]):
		if type(dts).__name__!='list': dts=[]
		if type(lats).__name__!='list': lats=[]
		if type(lons).__name__!='list': lons=[]
		while len(dts)<2: dts+=[None]
		#
		newcat=self.getTimeRangeCat(fullcat, dts[0], dts[1])
		newcat=self.getLatLonSubcat(newcat, lats, lons, llcols)
		#
		return newcat
	def addxytSubcat(self, subcatname='xytsubcat', fullcat=None, dts=[], lats=[], lons=[], llcols=[1,2]):
		self.subcats+=[[subcatname, self.getxytSubcat(fullcat, dts, lats, lons, llcols)]]
	
	def getTimeRangeCat(self, fullcat=None, dtFrom=None, dtTo=None):
		if fullcat==None: fullcat=self.cat
		if dtFrom==None: dtFrom=fullcat[0][0]
		if dtTo==None: dtTo=fullcat[-1][0]
		newcat=[]
		for rw in fullcat:
			if rw[0]>=dtFrom and rw[0]<=dtTo: newcat+=[rw]
		#
		return newcat
	def addTimeRangeCat(self, subcatname='dtSubcat', fullcat=None, dtFrom=None, dtTo=None):
		self.subcats+=[[subcatname, self.getTimeRangeCat(fullcat, dtFrom, dtTo)]]
	
	def getLatLonSubcat(self, fullcat, lats=[], lons=[], llcols=[1,2]):
		# llcols: lat, lon
		llrange=None
		if lats in [[], None]:
			if llrange==None: llrange=self.getLatLonRange(fullcat, latloncols=[llcols[0], llcols[1]])
			deltaLats=llrange[1][0]-float(llrange[0][0])
			lats=[llrange[0][0]+deltaLats/2.0, llrange[1][0]-deltaLats/2.0]
		
		if lons in [[], None]:
			if llrange==None: llrange=self.getLatLonRange(fullcat, latloncols=[llcols[0], llcols[1]])
			deltaLons=llrange[1][1]-float(llrange[0][1])
			lons=[llrange[0][1]+deltaLons/2.0, llrange[1][1]-deltaLons/2.0]
			
			#lats={get min, max lat,lon from catalog}
		# and same for lons...
		#
		newcat=[]
		for rw in fullcat:
			if (rw[llcols[0]]>=lats[0] and rw[llcols[0]]<=lats[1]) and (rw[llcols[1]]>=lons[0] and rw[llcols[1]]<=lons[1]):
				newcat+=[rw]
			#
		#
		return newcat
	def addLatLonSubcat(self, subcatname='xysubcat', fullcat=None, lats=[], lons=[], llcols=[1,2]):
		self.subcats+=[[subcatname, self.getLatLonSubcat(fullcat, lats, lons, llcols)]] 
	
	def getxytmSubcat(self, fullcat=None, dts=[], lats=[], lons=[], minmag=2.5, llmcols=[1,2,3]):
		# just do the whole thing here, so it's fast:
		if type(dts).__name__!='list': dts=[]
		if type(lats).__name__!='list': lats=[]
		if type(lons).__name__!='list': lons=[]
		if type(llmcols).__name__!='list': llmcols=[1,2,3]
		while len(dts)<2: dts+=[None]
		while len(llmcols)<3: llmcols+=[llmcols[-1]+1]
		llrange=None
		#return llmcols
		if lats in [[], None] or len(lats)!=2:
			if llrange==None: llrange=self.getLatLonRange(fullcat, latloncols=[llmcols[0], llmcols[1]])
			deltaLats=llrange[1][0]-float(llrange[0][0])
			lats=[llrange[0][0]+deltaLats/2.0, llrange[1][0]-deltaLats/2.0]
		
		if lons in [[], None] or len(lons)!=2:
			if llrange==None: llrange=self.getLatLonRange(fullcat, latloncols=[llmcols[0], llmcols[1]])
			deltaLons=llrange[1][1]-float(llrange[0][1])
			lons=[llrange[0][1]+deltaLons/2.0, llrange[1][1]-deltaLons/2.0]
		#
		#newcat=self.getTimeRangeCat(fullcat, dts[0], dts[1])
		#newcat=self.getLatLonSubcat(newcat, lats, lons, llcols)
		newcat=[]
		print lats, lons, dts
		for rw in fullcat:
			if rw[llmcols[0]]>=lats[0] and rw[llmcols[0]]<=lats[1] and rw[llmcols[1]]>=lons[0] and rw[llmcols[1]]<=lons[1] and rw[llmcols[2]]>=minmag and rw[0]>=dts[0] and rw[0]<=dts[1]:
				newcat+=[rw]
		return newcat
		#
	def addxytmSubcat(self, subcatname='xytmsubcat', fullcat=None, dts=[], lats=[], lons=[], minmag=2.5, llmcols=[1,2,3]):
		#print llmcols
		self.subcats+=[[subcatname, self.getxytmSubcat(fullcat, dts, lats, lons, minmag, llmcols)]]
			
	def mapOverlay(self, catalog=None, fignum=0, dots='b.', doShow=False):
		# this does not quite work yet. the map does not rescale properly for the distinct catalogs with different lat/lon ranges.
		# it looks like a good approach might be to create a map-class, which can contain a catalog or vice-versa, or maybe one
		# could be a sub-class, but i dont' think that hierarchy is clear.
		# the basic idea: map: lat/lon range, lat/lon center, projection, etc., catalogOverlays [] (are a list of catalogs overlayed on the map. note
		# the lat/lon range will be (at least) max/min(lon/lat from any cat) overlayed onto the map). also, annotationOverlays (text, graphics, etc.),
		# other stuff too...
		if catalog==None: catalog=self.cat
		f0=plt.figure(fignum)
		#
		#set up map:
		llr=self.getLatLonRange(catalog)	# latLonRange
		llr[0][0]-=2.0
		llr[0][1]-=2.0
		llr[1][0]+=2.0
		llr[1][1]+=2.0
		
		cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
		catmap=Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution=self.mapres, projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
		canvas=FigureCanvas(f0)
		catmap.ax=f0.add_axes([0,0,1,1])
		f0.set_figsize_inches((8/catmap.aspect,8.))
		#
		catmap.drawcoastlines(color='gray')
		catmap.drawcountries(color='gray')
		catmap.fillcontinents(color='beige')
		xfull, yfull=catmap(map(operator.itemgetter(2), catalog), map(operator.itemgetter(1), catalog))
		#epx, epy=catmap(epicenter[0], epicenter[1])
		catmap.plot(xfull, yfull, dots, label='Full Catalog')
		#catmap.plot(epx, epy, 'ro')
		#canvas.print_figure(saveName)
		
		if doShow: plt.show()
		
		return None
		
	def plotCatsMap(self, catalogses=None, maincat=0, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='best', maincatname='full cat', fignum=0, doCLF=True, bigmag=6.0, padfactor=.25):
		# somehow, this is returning a skewed map - i think basically, the basemap object re-callibrates itself to the smaller catalog, so x,y=thisthing.basemapobject(lat, lon) returns something off by a bit.
		
		# same as plotCatMap, but multiple catalogs. we assume the lat/lon range comes from the first catalog.
		# maincat is the "main catalog", the subcat we care about most. we assume the primary catalog is the broadest; maincat contains the epicenter, etc.
		if catalogses==None: catalogses=[maincatname, self.cat] + self.subcats

		#catalogs=[self.cat] + map(operator.itemgetter(1), self.subcats)
		#catnames=[maincatname] + map(operator.itemgetter(0), self.subcats)
		catalogs=map(operator.itemgetter(1), catalogses)
		catnames=map(operator.itemgetter(0), catalogses)
		#return [catalogs, catnames]
		catalog=catalogs[0]
		
		if epicenter==None:
			#mainshock=self.getMainEvent(catalog)
			mainshock=self.getMainEvent(catalogs[maincat])
			epicenter=[mainshock[2], mainshock[1]]
		#
		f0=plt.figure(fignum)	
		if doCLF: plt.clf()
		#		
		#set up map:
		llr=self.getLatLonRange(catalog)	# latLonRange #return [[minLat, minLon], [maxLat, maxLon]]
		latpad=padfactor*(llr[1][0]-llr[0][0])
		lonpad=padfactor*(llr[1][1]-llr[0][1])
		llr[0][0]-= latpad	#.5
		#if llr[0][0]<90.: llr[0][0]=90.
		llr[0][1]-= lonpad	#.5
		#if llr[0][1]<-180.: llr[0][1]=-180.
		llr[1][0]+= latpad	#.5
		#if llr[1][0]>90.: llr[1][0]=90.
		llr[1][1]+= latpad	#.5
		#if llr[1][1]>180.: llr[1][1]=180.
		
		print "setting up map prams"
		
		cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
		print "create basmap object."
		catmap=Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution =self.mapres, projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
		print "bm object created..."
		canvas=FigureCanvas(f0)
		catmap.ax=f0.add_axes([0,0,1,1])
		#f0.set_figsize_inches((8/catmap.aspect,8.))
		
		#f0.set_figsize_inches((10/catmap.aspect,10.))
		#f0.set_size_inches((10/catmap.aspect,10.))
		f0.set_size_inches((10.,15.))
		#
		print "draw stuff on map..."
		catmap.drawcoastlines(color='gray', zorder=0)
		catmap.drawcountries(color='gray', zorder=0)
		catmap.fillcontinents(color='beige', zorder=0)
		#catmap.drawrivers(color='b')
		catmap.drawstates()
		catmap.drawmeridians(range(int(llr[0][1]-2.0), int(llr[1][1]+2.0)), color='k', labels=[1,1,1,1])
		catmap.drawparallels(range(int(llr[0][0]-2.0), int(llr[1][0]+2.0)), color='k', labels=[1, 1, 1, 1])
		#
		'''
		catmap.llcrnrlon=llr[0][1]+2.0
		catmap.llcrnrlat=llr[0][0]+2.0
		catmap.urcrnrlon=llr[1][1]-2.0
		catmap.urcrnrlat=llr[1][0]-2.0
		'''
		print "plot catalogs..."
		icat=0
		for ct in catalogs:
			xfull, yfull=catmap(map(operator.itemgetter(2), ct), map(operator.itemgetter(1), ct))
			catmap.plot(xfull, yfull, '.', label='%s' % catnames[icat], ms=2, zorder=1, alpha=.5)
			icat+=1
		
		# now, plot all events m>m0 from the full catalog:
		#bigmag=5.0
		for rw in catalog:
			if rw[3]<bigmag: continue
			thisx, thisy=catmap(rw[2], rw[1])
			catmap.plot(thisx, thisy, '*', label='%s, %s\n (%s, %s)' % (str(rw[3]), str(rw[0]), str(rw[2]), str(rw[1])), ms=15, zorder=2)
			
		#epx, epy=catmap(epicenter[0], epicenter[1])
		#catmap.plot(epx, epy, 'ro', label='epicenter', zorder=1)
		#
		#################
		#################
		# this is how to draw an ellipse... obviously, this does not really belong in this part of the script;
		#  it was part of the learning process...
		###############
		#canvas.print_figure(saveName)
		#from matplotlib.patches import Ellipse
		#f=plt.figure(0)
		#
		#ax1=f0.gca()
		#el = Ellipse([-120.5, 35.9], .8, .3, -40, facecolor='b', alpha=0.4)
		#Xel, Yel = catmap(el.get_verts()[:,0],el.get_verts()[:,1])
		#catmap.plot(Xel, Yel, '-r', lw=2)
		#catmap.ax.fill(Xel, Yel, ec='r', fc='r', alpha=.4)
		###
		##################
		#ax1.add_artist(el)
		#catmap.ax.add_artist(el)
		#
		#ax=plt.gca()
		#el = Ellipse((self.tLon, self.tLat), 2.0*self.tA, 2.0*self.tB, -self.tTheta, facecolor='b', alpha=0.4)
		#catmap.ax.add_artist(el)
		#ax.add_artist(el)
		#
		#plt.plot(map(operator.itemgetter(2), self.fullCat), map(operator.itemgetter(1), self.fullCat), '+')
		#plt.plot(map(operator.itemgetter(2), self.shockCat), map(operator.itemgetter(1), self.shockCat), '.')
		#plt.plot(map(operator.itemgetter(2), fcat), map(operator.itemgetter(1), fcat), '+', label='Full Catalog')
		#plt.plot(map(operator.itemgetter(2), scat), map(operator.itemgetter(1), scat), '.', label='Aftershock zone')
		#plt.plot([epicenter[0]], [epicenter[1]], 'ro', label='epicenter')
		plt.legend(loc=legendLoc, numpoints=1)
		if doSave: plt.savefig('pltsave-%s' % saveName)
		
		canvas.print_figure(saveName)
		
		if doShow: plt.show()
		
		return catmap
	
	def plotCatMap(self, catalog=None, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='upper left', doCLF=True, eqicon='b,', myaxis=None, fignum=None, padfactor=.25, plotevents=True):
		# temporary:
		padfactor=0.
		if catalog==None: catalog=self.cat
		
		if epicenter==None:
			mainshock=self.getMainEvent(catalog)
			epicenter=[mainshock[2], mainshock[1]]
		#
		if doShow>=1 and fignum==None: fnum=doShow
		if fignum!=None: fnum=fignum
		
		f0=plt.figure(int(doShow))	
		if doCLF: plt.clf()
		#		
		#set up map:
		#llr=self.getLatLonRange(catalog)	# latLonRange
		llr=self.getLatLonRange(catalog)	# latLonRange #return [[minLat, minLon], [maxLat, maxLon]]
		latpad=padfactor*(llr[1][0]-llr[0][0])
		lonpad=padfactor*(llr[1][1]-llr[0][1])
		llr[0][0]-= latpad	#.5
		#if llr[0][0]<90.: llr[0][0]=90.
		llr[0][1]-= lonpad	#.5
		#if llr[0][1]<-180.: llr[0][1]=-180.
		llr[1][0]+= latpad	#.5
		#if llr[1][0]>90.: llr[1][0]=90.
		llr[1][1]+= latpad	#.5
		#if llr[1][1]>180.: llr[1][1]=180.
		
		cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
		if self.catmap==None: self.catmap=Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution=self.mapres, projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
		catmap=self.catmap
		
		canvas=FigureCanvas(f0)
		if myaxis==None: myaxis=f0.add_axes([0,0,1,1])
		#catmap.ax=f0.add_axes([0,0,1,1])
		catmap.ax=myaxis
		#f0.set_figsize_inches((8/catmap.aspect,8.))
		#
		catmap.drawcoastlines(color='gray', zorder=0)
		catmap.drawcountries(color='gray', zorder=0)
		catmap.drawstates(color='gray', zorder=0)
		catmap.drawrivers(color='gray', zorder=0)
		catmap.fillcontinents(color='beige', zorder=0)
		
		catmap.drawmeridians(range(int(llr[0][1]-2.0), int(llr[1][1]+2.0)), color='k', labels=[1,1,1,1])
		catmap.drawparallels(range(int(llr[0][0]-2.0), int(llr[1][0]+2.0)), color='k', labels=[1, 1, 1, 1])
		
		if plotevents:
			xfull, yfull=catmap(map(operator.itemgetter(2), catalog), map(operator.itemgetter(1), catalog))
			epx, epy=catmap(epicenter[0], epicenter[1])
			#catmap.plot(xfull, yfull, 'b,', label='Full Catalog')
			catmap.plot(xfull, yfull, eqicon, label='earthquakes', alpha=.5, zorder=2)
			catmap.plot(epx, epy, 'ro', zorder=2)
		
		# if we are inclned to save:
		if doSave and saveName!=None: canvas.print_figure(saveName)
		
		#
		#ax=plt.gca()
		#el = Ellipse((self.tLon, self.tLat), 2.0*self.tA, 2.0*self.tB, -self.tTheta, facecolor='b', alpha=0.4)
		#catmap.ax.add_artist(el)
		#ax.add_artist(el)
		#
		#plt.plot(map(operator.itemgetter(2), self.fullCat), map(operator.itemgetter(1), self.fullCat), '+')
		#plt.plot(map(operator.itemgetter(2), self.shockCat), map(operator.itemgetter(1), self.shockCat), '.')
		#plt.plot(map(operator.itemgetter(2), fcat), map(operator.itemgetter(1), fcat), '+', label='Full Catalog')
		#plt.plot(map(operator.itemgetter(2), scat), map(operator.itemgetter(1), scat), '.', label='Aftershock zone')
		#plt.plot([epicenter[0]], [epicenter[1]], 'ro', label='epicenter')
		plt.legend(loc=legendLoc, numpoints=1)
		if doSave: plt.savefig(saveName)
		if doShow: plt.show()
		
		return catmap

	def rbomoriQuadPlot(self, catnum=0, mc=2.5, winlen=501, rbthresh=1.0, bigmag=5.0, fignum=0, intlist=None, rbavelen=1, thislw=1.0, mainEV=None, plotevents=False):
		# make an awesome quad plot: omor-times, rbRatios, mag-seismicity, map-catalog
		# catalog -> a yodapy.eqcatalog() object
		# thislw: linewidth
		#
		#rbavelen=1	# rb-averaging length (aka, <nrb>_thisnum
		#print "winlen: %d" % winlen
		#if catalog==None: catalog=getMexicaliCat()
		#intlist=[25, 256, 512]
		if intlist==None: intlist=[int(winlen/2), winlen, winlen*2]
		catalog=self
	
		plt.figure(fignum)
		plt.clf()
		#ax0=plt.axes([.1,.1,.85, .35])
		# define axis (subplots) boundaries:
		ydelim=.03
		xdelim=.05
		xTS=[0.05, .5]
		yTS0=0.05
		dyTS=.3
			
		myaxes=[]
		nax=0
		# mags:
		x0=xTS[0]
		y0=yTS0+dyTS*nax
		#myaxes+=[plt.axes([xTS[0], y0, xTS[1], dyTS])]
		myaxes+=[plt.axes([.1, .03, .45, .3])]
		nax+=1
		# intervals:
		#x0=xTS[0]
		y0=yTS0+dyTS*nax
		#myaxes+=[plt.axes([xTS[0], y0, xTS[1], dyTS], sharex=myaxes[0])]
		myaxes+=[plt.axes([.1, .37, .45, .27], sharex=myaxes[0])]
		nax+=1
		# ratios:
		#x0=xTS[0]
		y0=yTS0+dyTS*nax
		#myaxes+=[plt.axes([xTS[0], y0, xTS[1], dyTS], sharex=myaxes[0])]
		myaxes+=[plt.axes([.1, .68, .45, .3], sharex=myaxes[0])]
		#
		# map:
		nax+=1
		xs=[xTS[1]+xdelim, .95]
		ys=[yTS0, 1.0]
		#myaxes+=[plt.axes([xs[0], xs[1], ys[0], ys[1]])]
		#myaxes+=[plt.axes([.6, .05, .35, .90], sharex=myaxes[0])]
		myaxes+=[plt.axes([.6, .05, .35, .90])]
		#
		# get RB ratios:
		try:
			catalog.rb
		except:
			catalog.rb=rbi.intervalRecordBreaker(None)
		#ratios=self.getIntervalRatios(minmag, windowLen, cat0, deltaipos, avlen)
	#	ratios=catalog.rb.getIntervalRatios(mc, winlen, catalog.getcat(catnum), 1)
		#
		#plotIntervals(self, intervals=[10, 100, 1000], minmag=2.0, catalog=None, fignum=0, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None], thisAxes=None):
		#catalog.plotIntervals(intlist, mc, catalog.getcat(catnum), fignum, [0,1,2,3], [None, None], [myaxes[0], myaxes[1]])
		#format: plotInts(self, intervals=[10, 100, 1000], catalog=None, minmag=2.0, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None]):
		epicen=None
		if mainEV!=None: epicen=[mainEV[2], mainEV[1]]
		#
		catalog.plotInts(intervals=intlist, catalog=catalog.getcat(catnum), minmag=mc, ax=myaxes[1], thislw=thislw, legendPos='upper left')
		#myaxes[1].set_label('mean intervals $\\tau$')
		myaxes[1].set_ylabel('mean intervals $\\ < \\tau >$', size=14)
		#
		# format: #plotMags(self, catalog=None, minmag=2.0, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None]):	
		catalog.plotMags(catalog.getcat(catnum), mc, myaxes[0])
		myaxes[0].set_ylabel('mag', size=14)
		#
		# plotIntervalRatiosAx(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, thisAx=None, ratios=None, deltaipos=1, avlen=1, mainEV=None)
		catalog.rb.plotIntervalRatiosAx(minmag=mc, windowLen=winlen, cat0=catalog.getcat(catnum), hitThreshold=rbthresh, bigmag=bigmag, thisAx=myaxes[2], ratios=None, deltaipos=1, avlen=rbavelen, mainEV=mainEV)
		myaxes[2].set_ylabel('RB ratio $r(N=%d)$' % winlen, size=14)
	
		#plt.figure(fignum)
		myfsize=12
		myaxes[1].text(.1, .1, '$<\\tau>$', None, rotation='vertical', size=myfsize)
		myaxes[0].text(.1, .1, 'mags', None, rotation='vertical', size=myfsize)
		myaxes[2].text(.1, .1, '$r(1) = frac{n_{rb-large}}{n{rb-small}}$', None, rotation='vertical', size=myfsize)
	
		#
		# plotCatMap(self, catalog=None, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='upper left', doCLF=True, eqicon='b,', myaxis=None, fignum=None, padfactor=.25, plotevents=True)
		X=catalog.plotCatMap(catalog=catalog.getcat(catnum), doShow=True, doSave=False, saveName=None, epicenter=epicen, legendLoc='best', doCLF=False, eqicon='b,', myaxis=myaxes[3], fignum=0, padfactor=.15, plotevents=plotevents)
		# and large events:
		bigEq=[]
		bigeqindex=0
		#print catalog.catmap
		for rw in catalog.getcat(catnum):
			if rw[3]>=bigmag and rw[0]>=catalog.getMainEvent(catalog.getcat(catnum))[0]:
				eqx,eqy = catalog.catmap(rw[2], rw[1])
				#catalog.catmap.plot(eqx, eqy, '*', label='m=%f, %s' % (rw[3], str(rw[0])))
				#catalog.catmap.plot(eqx, eqy, '*', label='m=%.2f, %d' % (rw[3], bigeqindex), ms=15)
				catalog.catmap.plot(eqx, eqy, '*', label='m=%.2f' % (rw[3]), ms=15)
				bigeqindex+=1
				print "eq: %d, %f, %s" % (bigeqindex, rw[3], str(rw[0]))
		catalog.catmap.ax.legend(loc='best', numpoints=1)
	#	#if len(bigEqs)>0:
		#	#XX=catalog.plotCatMap(bigEqs, True, False, None, 'best', False, '*', myaxes[3], fignum)
		#	XX=catalog.plotCatMap(bigEqs, True, False, None, None, 'best', False, '*', myaxes[3], fignum)
		# sharex=ax0
		#
		#
		plt.show()
	
		# problem: the x-axis of RB-ratios is in float form, formatted as datetimes (basically using "plot" type functions). the intervals/magnitudes are date_plot
		# functions. in short, they have different x-axis variable types, so we can't share the rb-ratios x-axis with the other two. it probably makes sense to convert
		# the interval-ratio plots to the rb-ratio style (because the rb-rations uses the "fillbetween" function).
	
		catalog.axlist=myaxes
		return catalog

	
	def testMap(self):
		import cPickle
		import time
		#
		fig=plt.figure()
		#
		t1 = time.clock()
		#m = Basemap(width=920000, height=1100000, resolution='f', projection='tmerc', lon_0=-4.2, lat_0=54.6)
		#m = Basemap(llcrnrlon=-11.0, llcrnrlat=45.0, urcrnrlon=3.0, urcrnrlat=59.0, resolution='f', projection='tmerc', lon_0=-4.2, lat_0=54.6)
		
		lllon=-115
		lllat=32
		urlon=-105
		urlat=42
		lon0=lllon + (urlon-lllon)/2.0
		lat0=lllat + (urlat-urlat)/2.0
		print "center: %f, %f" % (lon0, lat0)
		m = Basemap(llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat, resolution=self.mapres, projection='tmerc', lon_0=lon0, lat_0=lat0)
		m.drawcountries()
		m.drawrivers()
		print time.clock()-t1,' secs to create original Basemap instance'

		# cPickle the class instance.
		cPickle.dump(m,open('map.pickle','wb'),-1)

		# clear the figure
		plt.clf()
		# read cPickle back in and plot it again (should be much faster).
		t1 = time.clock()
		m2 = cPickle.load(open('map.pickle','rb'))
		# draw coastlines and fill continents.
		m.drawcoastlines()
		# fill continents and lakes
		m.fillcontinents(color='coral',lake_color='aqua')
		# draw political boundaries.
		m.drawcountries(linewidth=1)
		# fill map projection region light blue (this will
		# paint ocean areas same color as lakes).
		m.drawmapboundary(fill_color='aqua')
		# draw major rivers.
		m.drawrivers(color='b')
		print time.clock()-t1,' secs to plot using using a pickled Basemap instance'
		# draw parallels
		circles = np.arange(48,65,2).tolist()
		m.drawparallels(circles,labels=[1,1,0,0])
		# draw meridians
		meridians = np.arange(-12,13,2)
		m.drawmeridians(meridians,labels=[0,0,1,1])
		plt.title("High-Res British Isles",y=1.04)
		plt.show()
	
	def setSpecialCatSQL(self, catname='parkfield'):
		#reload(yp)
		if catname in ['parkfield', 'pf', 'PF', 'park']:
			#theta=40.0, clat=35.9, clon=-120.5, ra=.4, rb=.15
			#self.setCatFromSQL(dtm.datetime(1969,1,1), dtm.datetime.now(), [34.9, 36.9], [-121.5, -119.5], 1.5, "Earthquakes", 523, 'asc')
			self.setCatFromSQL(dtm.datetime(1972,1,1), dtm.datetime.now(), [34.4, 37.4], [-121.11, -119.4], 1.5, "Earthquakes", 523, 'asc')
			self.addEllipCat('PFshock (.8 x .15)', self.cat, 40.0, 35.9, -120.5, 0.8, 0.15)
			self.addEllipCat('PFshock (.4 x .15)', self.cat, 40.0, 35.9, -120.5, 0.4, 0.15)
		if catname == 'PF5yr':
			#theta=40.0, clat=35.9, clon=-120.5, ra=.4, rb=.15
			#self.setCatFromSQL(dtm.datetime(1969,1,1), dtm.datetime.now(), [34.9, 36.9], [-121.5, -119.5], 1.5, "Earthquakes", 523, 'asc')
			self.setCatFromSQL(dtm.datetime(1999, 9, 28), dtm.datetime(2009,9,29), [34.4, 37.4], [-121.11, -119.4], 1.5, "Earthquakes", 523, 'asc')
			self.addEllipCat('PFshock (.8 x .15)', self.cat, 40.0, 35.9, -120.5, 0.8, 0.15)
			self.addEllipCat('PFshock (.4 x .15)', self.cat, 40.0, 35.9, -120.5, 0.4, 0.15)
		if catname in ['taiwan']:
			self.setCatFromSQL(dtm.datetime(1980,1,1), dtm.datetime(2010,6,1), [-90, 90], [-180, 180], 2.0, 'Earthquakes', 21)
			
	def setCatFromSQL(self, startDate=dtm.datetime(1999,9,28, 17,15,24), endDate=dtm.datetime(2009,9,28, 17,15,24), lats=[32.0, 37.0], lons=[-125.0, -115.0], minmag=3.0, catalogName='Earthquakes', catalogID=523, ordering='asc'):
		self.cat=self.getCatFromSQL(startDate, endDate, lats, lons, minmag, catalogName, catalogID, ordering)
		return None
	
	def getCatFromSQL(self, startDate=dtm.datetime(1999,9,28, 17,15,24), endDate=dtm.datetime(2009,9,28, 17,15,24), lats=[32.0, 37.0], lons=[-125.0, -115.0], minmag=2.0, catalogName='Earthquakes', catalogID=523, ordering='asc'):
		# return a catalog:
		if lats[0]>lats[1]: lats.reverse()
		if lons[0]>lons[1]: lons.reverse()
		#if yp.datetimeToFloat(startDate)>yp.datetimeToFloat(endDate):
		if startDate>endDate:
			middledate=startDate
			startDate=endDate
			endDate=middledate
			middledate=None
		if ordering not in ['asc', 'desc']: ordering='desc'
		
		import _mysql
		import MySQLdb
		#
		#sqlHost = 'localhost'
		sqlHost = self.sqlhost
		sqlUser = 'myoder'
		sqlPassword = 'yoda'
		sqlPort = self.sqlport
		sqlDB = 'QuakeData'
		con=MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
		c1=con.cursor()
		sqlstr='select eventDateTime, lat, lon, mag from %s where catalogID=%d and lat between %f and %f and lon between %f and %f and mag>=%f and eventDateTime between \'%s\' and \'%s\' order by eventDateTime %s' % (catalogName, catalogID, lats[0], lats[1], lons[0], lons[1], minmag, str(startDate), str(endDate), ordering)
		catList=[]
		#print sqlstr
		#
		c1.execute(sqlstr)
		rw=c1.fetchone()
		while rw!=None:
		#	# spin through the cursor; write a catalog. note formatting choices...
			catList+=[[rw[0], float(rw[1]), float(rw[2]), float(rw[3])]]
			rw=c1.fetchone()
		#catList=self.fetchall()
		c1.close()
		con.close()
		# now we have a catalog of the parkfield area (note, it is partially defined by our "parkfieldquakes" MySQL view.
		#
		#makeShockCat(incat, outcat)
		#makeShockCat(fullcatout, shockcatout)
		return catList
	
	def plotGRdistsFromTo(self, frmDt=None, toDt=None, catlist=None):
		# plot GRdist between two dates for a catalog.
		# assume standard catalog format.
		if catlist==None: catlist=self.cat
		if frmDt==None: frmDt=catlist[0][0]
		if toDt==None: toDt=catlist[-1][0]
		#
		mymags=[]
		for rw in catlist:
			if rw[0]>=frmDt and rw[0]<=toDt: mymags+=[rw[3]]
			# if rw[0]>toDt: break	# unless the cat is out of sequence. computers are fast; just spin through for now.
		#
		return self.plotGRdist(mymags)
		
	def plotGRdist(self, mags=None, doShow=True, fname='GRdist.png', plotTitle="Magnitude Distribution", fignum=0):
		# mags: a 1D array of magnitudes
		if mags==None: mags=map(operator.itemgetter(3), self.cat)
		# if mags rows are not scalar, assume a full standard type catalog has been passed.
		try:
			if len(mags[0])>=3: mags=map(operator.itemgetter(3), mags)
		except TypeError:
			# a list of scalars will throw a "can't get len." error. we should be able to skip without doing anything at all.
			# maybe a better approach is to test the type of mags[0] for list or tuple...
			dummyvar=None	# place-holder
		#
		mags.sort()
		# get rid of biggest event (probably a large off-GR earthquake):
		#mags.pop()
		#mags.reverse()
		#print mags
		#print len(mags)
		
		if doShow==True or fname!=None:
			# make a plot and show and/or save
			#Y=range(1, len(mags)+1)
			Y=frange(1, len(mags), -1)
			#print Y
			#print len(Y)
			plt.figure(fignum)
			plt.clf()
			plt.semilogy(mags, Y, '.-')
			plt.xlabel("Magnitude, m")
			plt.ylabel("Number of Events, n")
			plt.title(plotTitle)
			if fname!=None: plt.savefig(fname)
			if doShow: plt.show()
		
		return mags	

	def getIntervals(self, catList, winLen):
		catLen=len(catList)
		i=(catLen-1-winLen)	# start winLen positions from the end.
		thisInterval=0
		#N=1
		intervals=[]	# [[eventDateTime, totalInterval]]
		while i>=0:
			#
			thisInterval=datetimeToFloat(catList[i+winLen][0])-datetimeToFloat(catList[i][0])
			intervals+=[[catList[i+winLen][0], thisInterval]]
			i-=1
		
		#
		#return [intervals, catList]
		return intervals
	#
	def plotIntervals(self, intervals=[10, 100, 1000], minmag=2.0, catalog=None, fignum=0, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None], thisAxes=None):
		if type(plotDates).__name__!='list': plotDates=[None, None]
		while len(plotDates)<2: plotDates+=[None]
		#
		if catalog==None: catalog=self.cat
		#X = plotIntervals(intervals, minmag, catalog, fignum, dtmlatlonmagCols)
		#return X
		#zonedat=[35.9, -120.5, .4, .15, 40.0]
		cols=dtmlatlonmagCols	# for efficient notation
		#zonedat=[35.9, -120.5, .4, .05, 40.0]	# this will be done in advance of this function call, when the catalog is made.
		#minmag=2.0
		#dts=['1950-01-01', str(dtm.datetime.now())]
		#sqlcat="Earthquakes"
		#catid=523
	
		#
		plt.figure(fignum)
		if thisAxes==None:
			#plt.figure(fignum)
			plt.clf()
			ax0=plt.axes([.1,.1,.85, .35])
			plt.xlabel("time")
			plt.ylabel("mags")
			ax1=plt.axes([.1, .55, .85, .35], sharex=ax0)
			plt.ylabel("mean interval")
			plt.xlabel("")
			plt.title("Mean intervals, $m_c=%s$" % str(minmag))
		else:
			ax0=thisAxes[0]
			ax1=thisAxes[1]
	
		#dtms=map(operator.itemgetter(cols[0]), catalog)
		#lats=map(operator.itemgetter(cols[1]), catalog)
		#lons=map(operator.itemgetter(cols[2]), catalog)
		#mags=map(operator.itemgetter(cols[3]), catalog)
		mags=[]
		activecat=[]
		for rw in catalog:
			if rw[cols[3]]<minmag: continue
			mags+=[[rw[cols[0]], rw[cols[3]]]]
			activecat+=[rw]
		mags=vlinePadList(mags, minmag-abs(minmag)*.1)	# return the mags data padded for vertical line style plotting. this is just a trick to get width=1 histograms.
		#
		ax0.plot_date(map(operator.itemgetter(0), mags), map(operator.itemgetter(1), mags), '-')
		shockints=[]
		
		#print "plotdates: %s" % str(plotDates)
		for wlen in intervals:
			#print "plotting for wlen: %d" % wlen
			##
			#shockints+=[getIntervals(catalog, wlen)]
			shockints+=[self.getIntervals(activecat, wlen)]
			#
			# trim off max/min date ends for prettier plots:
			#print "mindt: %s" % str(shockints[-1][0])
			while (plotDates[1]!=None and plotDates[1]<shockints[-1][0][0]): a=shockints[-1].pop(0)
			while plotDates[0]!=None and plotDates[0]>shockints[-1][-1][0]: a=shockints[-1].pop()
			#
			#plt.plot(map(operator.itemgetter(0), shockints[-1]), scipy.array(map(operator.itemgetter(1), shockints[-1]))/float(wlen), '-', label='winLen=%d' % wlen)
			#
			X=map(operator.itemgetter(0), shockints[-1])
			# pylab.date2num(dtm)
			#XX=date2num(X)
			ax1.plot(X, scipy.array(map(operator.itemgetter(1), shockints[-1]))/float(wlen), '-', label='$N=%d$' % wlen, lw=1.0)
			#ax1.semilogy(map(operator.itemgetter(0), shockints[-1]), scipy.array(map(operator.itemgetter(1), shockints[-1]))/float(wlen), '-', label='winLen=%d' % wlen)
			# fg.autofmt_xdate()
			
		
			#plt.plot(map(operator.itemgetter(0), shockints[-1]), scipy.array(map(operator.itemgetter(1), shockints[-1])), '-', label='winLen=%d' % wlen)
		#
		plt.legend(loc='upper left')
			
		#plt.legend(loc='lower left')
	
	
		plt.show()
	
		return shockints
		
	def plotInts(self, intervals=[10, 100, 1000], catalog=None, minmag=2.0, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None], thislw=1.0, legendPos='best'):
		# plot mean intervals only. having figured out how to put a bunch of independently generated plots onto a single canvas, we
		# split up some of these plots...
		#
		if type(plotDates).__name__!='list': plotDates=[None, None]
		while len(plotDates)<2: plotDates+=[None]
		#
		if catalog==None: catalog=self.cat
		cols=dtmlatlonmagCols	# for efficient notation
		activecat=[]
		for rw in catalog:
			if rw[cols[3]]>=minmag: activecat+=[rw]	# build active catalog of events m>mc
		#
		shockints=[]
		for wlen in intervals:
			shockints+=[self.getIntervals(activecat, wlen)]
			# trim off catalog elements outside min/max date range:
			while (plotDates[1]!=None and plotDates[1]<shockints[-1][0][0]): a=shockints[-1].pop(0)
			while plotDates[0]!=None and plotDates[0]>shockints[-1][-1][0]: a=shockints[-1].pop()
			#
			ax.plot(map(operator.itemgetter(0), shockints[-1]), scipy.array(map(operator.itemgetter(1), shockints[-1]))/float(wlen), '-', label='$N=%d$' % wlen, lw=thislw)
			#ax1.semilogy(map(operator.itemgetter(0), shockints[-1]), scipy.array(map(operator.itemgetter(1), shockints[-1]))/float(wlen), '-', label='winLen=%d' % wlen)
			
		
			#plt.plot(map(operator.itemgetter(0), shockints[-1]), scipy.array(map(operator.itemgetter(1), shockints[-1])), '-', label='winLen=%d' % wlen)
		#
		ax.legend(loc=legendPos)
			
		return shockints
		
	def plotMags(self, catalog=None, minmag=2.0, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None]):
		if type(plotDates).__name__!='list': plotDates=[None, None]
		while len(plotDates)<2: plotDates+=[None]
		#
		if catalog==None: catalog=self.cat
		cols=dtmlatlonmagCols	# for efficient notation
		#
		mags=[]
		activecat=[]
		for rw in catalog:
			if rw[cols[3]]>=minmag: mags+=[[rw[cols[0]], rw[cols[3]]]]
		mags=vlinePadList(mags, minmag-abs(minmag)*.1)	# return the mags data padded for vertical line style plotting. this is just a trick to get width=1 histograms.
		#
		ax.plot(map(operator.itemgetter(0), mags), map(operator.itemgetter(1), mags), '-')
		ax.legend(loc='best')
		return mags

def yodaecho(something=None):
	# this is for loadFileToH/VList() so we can have a null conversion (return whatever is there).
	return something

def loadFileToHlist(fname=None, castFunct=yodaecho):
	# in a later version, consider castFunct=[int, int, float, str, int...] so we can make a mask to convert
	if fname==None: return None
	try:
		x=castFunct(42.42)
	except:
		castFunct=yodaecho
		
	#
	f=open(fname, 'r')
	X=[]
	for rw in f:
		if rw[0] in ['#', ' ', '\t']: continue
		rws=rw.split('\t')
		X+=[[]]
		for elem in rws:
			X[-1]+=[castFunct(elem)]
			
		#
	#
	f.close()
	#
	return X

def floatYear(thisdt):
	yr=thisdt.year
	Ndays=dtm.datetime(yr,12,31).timetuple()[7]
	T=thisdt.timetuple()
	nday=T[7]-1
	dayfrac=(T[3]*60.*60.*10**6+T[4]*60*10**6 + T[5]*10**6 + thisdt.microsecond)/((24.*60.*60. + 60*60 + 60)*10**6)
	#
	return yr + (nday+dayfrac)/Ndays
floatyear=floatYear

def decistring(floatval, ndeciplaces):
	#return a string of a float with a fixed number of decimals
	strval=str(round(floatval, ndeciplaces))
	strvals=strval.split('.')
	return strvals[0]+'.'+strvals[1][0:ndeciplaces]

'''
def logSpacedLogPoints(XY, logbase=10.0):
	# we want [X,Y]. we'll figure out how to reshape [[x,y]...] arrays later.
	# this version for power-law (or exponential) data (aka, we'd plot it logX).
	#
	outsies=[[],[]]	#output array.
	X=XY[0]
	Y=XY[1]
	logsies=[int(math.log(X[0], logbase))]
	outsies[0]+=[X[0]]
	outsies[1]+=[Y[0]]
	for i in xrange(1, len(X)):
		logx=int(math.log(X[i], logbase))
		#if logx==math.log(outsies[0][-1], logbase): continue
		if logx==logsies[-1]: continue
		logsies+=[logx]
		outsies[0]+=[X[i]]
		outsies[1]+=[Y[i]]
	#
	return outsies
'''
fitMarkerShortList=['o', '^', 's', 'p', '*', 'h', '+', 'H', 'D', 'x']

def integerSpacedPoints(XY, intFactor=10):
	# linear data come in (presumably the log of PL/exp data), so the data density
	# are probably higher at one end of the distribution.
	# intFactor: basically, how many points between base-10 integers (effectively, base for log-binning).
	outsies=[[],[]]	#output array.
	X=XY[0]
	Y=XY[1]
	intsies=[int(X[0]/intFactor)]
	outsies[0]+=[X[0]]
	outsies[1]+=[Y[0]]
	#
	for i in xrange(1, len(X)):
		intx=int(X[i]/intFactor)
		if intx==intsies[-1]: continue
		intsies+=[intx]
		outsies[0]+=[X[i]]
		outsies[1]+=[Y[i]]
	return outsies

def catfromANSS(lon=[135., 150.], lat=[30., 41.5], minMag=4.0, dates0=[dtm.date(2005,01,01), None], Nmax=999999, fout='cats/mycat.cat'):
	# get a basic catalog. then, we'll do a poly-subcat. we need a consistent catalog.
	# eventually, cut up "japancatfromANSS()", etc. to call this base function and move to yodapy.
	if dates0[1]==None:
		# i think this needs a "date" object, and datetime breaks.
		# so, make a Now() for date.
		nowdtm=dtm.datetime.now()
		dates0[1]=dtm.date(nowdtm.year, nowdtm.month, nowdtm.day)
	#	
	catlist=getANSSlist(lon, lat, minMag, dates0, Nmax, None)
	f=open(fout, 'w')
	f.write("#anss catalog\n")
	f.write("#lon=%s\tlat=%s\tm0=%f\tdates=%s\n" % (str(lon), str(lat), minMag, str(dates0)))
	for rw in catlist:
		# simple, right? except that ANSS has a habit of writind useless date-times like "2001/10/08 24:00:07.62" (hour=24), or
		# where minute=60. we could toss these. for now, assume 2001/10/8 24:00:00 -> 2001/10/9/00:00:00. change by proper time-arithmetic.
		# first, parse the date-string:
		strDt, strTm=rw[0].split()[0], rw[0].split()[1]
		if '/' in strDt: delim='/'
		if '-' in strDt: delim='-'
		strDts=strDt.split(delim)
		strTms=strTm.split(':')
		yr=int(strDts[0])
		mnth=int(strDts[1])
		dy=int(strDts[2])
		hr=int(strTms[0])
		mn=int(strTms[1])
		sc=float(strTms[2])
		microsecs=(10**6)*sc%1.
		# one approach is to start with year, month and add all the subsequent quantities using datetime.timedelta objects, which we have to
		# do once we get into callendar addition anyway...
		#so let's assume the date part is correct:
		myDt=dtm.datetime(yr, mnth, dy)
		#mytimedelta=dtm.timedelta(hours=hr)
		myDt+=dtm.timedelta(hours=hr)
		myDt+=dtm.timedelta(minutes=mn)
		myDt+=dtm.timedelta(seconds=sc)
		myDt+=dtm.timedelta(microseconds=microsecs)
		#
		myDtStr='%d/%d/%d %d:%d:%d.%d' % (myDt.year, myDt.month, myDt.day, myDt.hour, myDt.minute, myDt.second, myDt.microsecond)	
		#
		#f.write('%s\t%s\t%s\t%s\n' % (rw[0], rw[1], rw[2], rw[4]))
		f.write('%s\t%s\t%s\t%s\n' % (myDtStr, rw[1], rw[2], rw[4]))
	f.close()
	
	 
	return catlist

		
