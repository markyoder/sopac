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
from matplotlib.patches import Ellipse
import matplotlib.patches as patches

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

#import MySQLdb
#from mpl_toolkits.basemap import Basemap as Basemap
import mpl_toolkits.basemap as mpbm
#from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.backends.backend_agg as mpcm

import yodapy as yp
#import rbIntervals as rbp

# script(s) to fetch, process, plot, etc. GPS data from SoPac.

def getStationPositions(sitecodes=['ucd1', 'ucsb'], yearSelected=2009, doy=11, datums='NONE', refFrame='RERUN ITRF2005', outputFormat='csv'):
	#
	# use the sopac website to fetch and process raw GPS data.
	# spoof the query page: http://sopac.ucsd.edu/cgi-bin/sector.cgi
	# for our purposes, this page posts to: /cgi-bin/sector.cgi
	# as per: <form method="post" action="/cgi-bin/sector.cgi" enctype="application/x-www-form-urlencoded" name="basicForm" expires="1d">
	#
	# this is how we do a similar operation to fetch ANSS earthquake data:
	
	#dates=[dtm.date(dates0[0].year, dates0[0].month, dates0[0].day), dtm.date(dates0[1].year, dates0[1].month, dates0[1].day)]
	#anssPrams={'format':'cnss', 'output':'readable', 'mintime':str(dates[0]).replace('-', '/'), 'maxtime':str(dates[1]).replace('-', '/'), 'minmag':str(minMag), 'minlat':lat[0], 'maxlat':lat[1], 'minlon':lon[0], 'maxlon':lon[1], 'etype':'E', 'searchlimit':Nmax}
	#f = urllib.urlopen('http://www.ncedc.org/cgi-bin/catalog-search2.pl', urllib.urlencode(anssPrams))
	#
	popupArrays='ALL'
	popupSites='NONE'
	#popupSites='ucd1'
	textSites=getTextSiteString(sitecodes) ###
	#textSites=''
	gpsDatum=datums
	refFrame=getRefFrame(refFrame)
	#refFrame="RERUN ITRF2005"
	#doy="011"
	#
	#sopacPrams={'hiddenGetCoord':1, 'hiddenReloadSites':0, 'popupArrays':popupArrays, 'popupSites':popupSites, 'textSites':textSites, 'yearSelected':str(yearSelected), 'daySelected':doy, 'decYearSelected':'NONE', 'datumSelected':gpsDatum, 'sourceName':refFrame, 'output':outputFormat }
	# it looks like SOPAC folks changed the name of the "output" variable from "output" to "format_type"
	sopacPrams={'hiddenGetCoord':1, 'hiddenReloadSites':0, 'popupArrays':popupArrays, 'popupSites':popupSites, 'textSites':textSites, 'yearSelected':str(yearSelected), 'daySelected':doy, 'decYearSelected':'NONE', 'datumSelected':gpsDatum, 'sourceName':refFrame, 'format_type':outputFormat }
	#
	#print sopacPrams
	#
	f = urllib.urlopen('http://sopac.ucsd.edu/cgi-bin/sector.cgi', urllib.urlencode(sopacPrams))
	
	returnData=[]
	for rw in f:
		returnData+=[rw]
	f.close()
	
	#printArrayToFile(returnData)
	
	return returnData

def getStationVelocities(stations=['ucd1', 'ucsb'], startDate=dtm.datetime(2000, 1,1), endDate=dtm.datetime.now(), timeStep=dtm.timedelta(days=1)):
	# get station velocities from raw parameters (station positions are not provided).
	# 1) get all the station positions
	# 2) calculate v_x_i, v_y_i for each station, interval
	# 3) <v_x>, <v_y>, var_x, var_y
	#
	stationPositions=getStationTraces(stations, startDate, endDate, timeStep)
	return getStationVelocities0(stationPositions)
	
def getStationVelocities0(stationPositions=None):
	# get station velocities from a set of station positions.
	#
	#if stationPositions==None: stationPositions=getStationPositions()
	if stationPositions==None: stationPositions=getStationTraces()
	#
	# position data are like: [["time", "sta1", "sta2", ...], [t
	stationVelocities=[stationPositions[0]]	# and return an array in the same format as stationPositions: [t, [vx, vy, vz, sigx, sizy, sigz]_1, []_2, ...]
	for i in xrange(2, len(stationPositions)):
		dt=float(stationPositions[i][0]-stationPositions[i-1][0])
		#stationVelocities+=[[dt]]
		stationVelocities+=[[stationPositions[i][0]]]
		for ii in xrange(1, len(stationPositions[i])):
			# for each set of station data in the row:
			P1=stationPositions[i][ii]
			P0=stationPositions[i-1][ii]	# current and previous position for this station (P1,0 for shorthand)
			# so these are [t, lat, lon, z, obsErrLat, oblErrLon, obsErrZ]
			stationVelocities[-1]+=[[(P1[0]-P0[0])/dt, (P1[1]-P0[1])/dt, (P1[2]-P0[2])/dt, P1[3]**2.+P0[3]**2., P1[4]**2.+P0[4]**2.,P1[5]**2.+P0[5]**2.]]	#[vx, vy, vz, varX, varY, varZ]
			#
		#
	# return format is like: [[dt, [vx, vy, vz, sigvx, sigvy, sigvz], ]
	return stationVelocities

def getMeanStationVelocities1(stationTraces=None):
	# get station velocities from a set of station positions.
	# the previous way we did this was probaly pretty silly - which was just to average (x_i-x_{i-1})/dt.
	# in this method:  get positions; fit to a linear function: lat = lat0 + v_lat*t, lon = lon0 + v_lon*t
	#
	#if stationPositions==None: stationPositions=getStationPositions()
	times=map(operator.itemgetter(0), stationTraces[1:])
	#positions=[times]
	returnVelocities=[]
	lflon=yp.linefit()
	lflat=yp.linefit()
	for i in xrange(1, len(stationTraces[1])):
		thispositions=map(operator.itemgetter(i), stationTraces)
		thislats=map(operator.itemgetter(0), thispositions[1:])
		thislons=map(operator.itemgetter(1), thispositions[1:])
		thisStation = thispositions[0]
		#print "lens: %d, %d, %d" % (len(times), len(thislats), len(thislons))
		#
		# now, fit x,y to lines.
		lflon.datas=[times, thislons, scipy.zeros(len(thislons))+1]
		lflat.datas=[times, thislats, scipy.zeros(len(thislats))+1]
		#print lflon.datas
		lonfits=lflon.doFit()
		latfits=lflat.doFit()
		#
		#returnVelocities+=[[thisStation, lflon.a, lflat.a, lflon.b, lflat.b, lflon.varaprime, lflat.varaprime, lflon.varbprime, lflat.varbprime]]
		# note that we have returned the a parameter as well. this can be used to infer the expected starting position (aka, we do not treat
		# the first measurement as being special). obviously, use times[0] to get the starting position. be careful, since the raw lat/lon outputs
		# (a values) will look "normal" since movement is slow and t0 is relatively recent.
		#
		# but let's go ahead and return the starting position, as per the data set provided.
		lon0=lflon.a + lflon.b*times[0]
		lat0=lflat.a + lflat.b*times[0]
		returnVelocities+=[[thisStation, lon0, lat0, lflon.b, lflat.b, lflon.varaprime, lflat.varaprime, lflon.varbprime, lflat.varbprime]]
	
	return returnVelocities
	

def getMeanStationVelocities(vels=None):
	if vels==None: vels=getStationVelocities()
	#
	times=scipy.array(map(operator.itemgetter(0), vels[1:]))
	
	stationses=[]
	for i in xrange(1, len(vels[1])):
		stationses+=[map(operator.itemgetter(i), vels)]
	#
	#print "returning stationses"
	#return stationses
	meanvels=[]
	for sta in stationses:
		meanvels+=[[sta[0]]]
		for ii in xrange(len(sta[1])):
			meanvels[-1]+=[numpy.average(map(operator.itemgetter(ii), sta[1:]))]	# average of all values; 3 velocities, 3 obs-errors
		for ii in xrange(3):
			meanvels[-1]+=[numpy.std(map(operator.itemgetter(ii), sta[1:]))**2.]	# statistical errors for the first 3 cols. (velocities).

	#return stationses
	# so the return data are like: [[station, vlat, vlon, vz, obsSigvlat, obsSigvlon, obsSigvz, statSigvlat, statSigvlon, statSigvz], ...]
	# the first set of errors is the mean-Variance of the reported positional error (sig1^2+sig2^2); the second set of errors is the statistical variance in the calculated velocities.
	return meanvels

def velocityEnchilada(stations=['FARB', 'P248', 'P217', 'P271', 'UCD1', 'CRU1', 'P475', 'POTR', 'SOLI', 'CUHS', 'P286', 'DLNO', 'LEMA'], startDate=dtm.datetime(2005, 1,1), endDate=dtm.datetime.now(), timeStep=dtm.timedelta(days=10)):
	# the whole enchilada in one call...
	# ?? , 'OVRO', 'MOJA', 'FIBR', 'PEAR'
	# note: we should probably change this to use getMeanStationVelocities1(traces), which uses a better line-fit
	# method to determine velocities.
	locs=getStationTraces(stations, startDate, endDate, timeStep)
	vels=getStationVelocities0(locs)
	#
	X=plotMeanStationVelocities(vels, locs, 10**9, None)
	#
	return [locs, vels]

def plotMeanStationVelocities(vels, pos, vectorTime=1000000., llrange=None):
	# velocities and positions. really, we only need the first row for each position (we'll strip the first 2 rows for each set of station positions, so we can
	# take the whole enchilada as well).
	#
	# get mean velocities; draw vectors from the first position of each station, then draw error-ellipses.
	# how do we scale the vectors (and error-ellipses)? we'll sort out something. maybe find the max/min velocity vectors and scale the representaiton to alpha*map extents.
	#
	#  llrange -> [lrlon, lrlat, urlon, urlat]
	# vectorTime: time to determine vector length.
	#
	myheaders=pos[0]
	startPoints=[]
	for i in xrange(1, len(pos[1])):
		startPoints+=[pos[1][i]]
	#print "start points: %s" % str(startPoints)
	#
	meanVels=getMeanStationVelocities(vels)
	#print "mean vels: %s" % str(meanVels)
	#
	#endPoints=[]
	#for istation in xrange(len(startPoints)):
	#	endPoints+=[[]]
	#	for i in xrange(3):
	#		endPoints[-1]+=[startPoints[istation][i] + vectorTime*meanVels[istation][i] ]
	
	#return [meanVels, startPoints]
	#
	# and plot the vectors:
	plt.ion()
	plt.figure(0)
	plt.clf()
	
	for i in xrange(len(startPoints)):
		#print "arrow prams: %f, %f, %e, %e" % (startPoints[i][0], startPoints[i][1], meanVels[i][1]*vectorTime, meanVels[i][2]*vectorTime)
		thiscolor=yp.pycolor(i)
		#plt.arrow(startPoints[i][0], startPoints[i][1], meanVels[i][1]*vectorTime, meanVels[i][2]*vectorTime, \
		#lw=3, ec=thiscolor, fc=thiscolor)
		plt.arrow(startPoints[i][1], startPoints[i][0], meanVels[i][2]*vectorTime, meanVels[i][1]*vectorTime, \
		lw=2, ec=thiscolor, fc=thiscolor, head_width=.05, head_length=.05 )
		#width=100, head_width=500, head_length=50)
		#plt.plot([startPoints[i][0], startPoints[i][0]+meanVels[i][1]*vectorTime], [startPoints[i][1], startPoints[i][1]+meanVels[i][2]*vectorTime] )
		plt.plot([startPoints[i][1], startPoints[i][1]+meanVels[i][2]*vectorTime], [startPoints[i][0], startPoints[i][0]+meanVels[i][1]*vectorTime], color=thiscolor)
	
	plt.figure(1)
	plt.clf()
	#
	# make basemap stuff
	# we'll need max extents:
	lats=map(operator.itemgetter(0), startPoints)
	lons=map(operator.itemgetter(1), startPoints)
	vlats=map(operator.itemgetter(1), meanVels)
	vlons=map(operator.itemgetter(2), meanVels)
	
	lats+=(scipy.array(lats)+scipy.array(vlats)*vectorTime).tolist()
	lons+=(scipy.array(lons)+scipy.array(vlons)*vectorTime).tolist()
	
	print lats, lons
	llpad=.25
	llr=[[min(lats)-llpad, min(lons)-llpad], [max(lats)+llpad, max(lons)+llpad]]
	print llr
	cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
	catmap=mpbm.Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution ='l', projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
	#
	thisax=plt.gca()
	catmap.ax=thisax
	#
	catmap.drawcoastlines(color='gray')
	catmap.drawcountries(color='gray')
	catmap.drawstates(color='gray')
	catmap.drawrivers(color='gray')
	catmap.fillcontinents(color='beige')
	
	catmap.drawmeridians(range(int(llr[0][1]-2.0), int(llr[1][1]+2.0)), color='k', labels=[1,1,1,1])
	catmap.drawparallels(range(int(llr[0][0]-2.0), int(llr[1][0]+2.0)), color='k', labels=[1, 1, 1, 1])
	#
	#return None
	for i in xrange(len(startPoints)):
		#print "arrow prams: %f, %f, %f, %f" % (startPoints[i][0], startPoints[i][1], meanVels[i][1], meanVels[i][2])
		thiscolor=yp.pycolor(i)

		# the arrow is more complicated here. get the dx, dy from X,Y below.
		X, Y=catmap([startPoints[i][1], startPoints[i][1]+meanVels[i][2]*vectorTime], [startPoints[i][0], startPoints[i][0]+meanVels[i][1]*vectorTime])
		arrowHeadX, arrowHeadY=catmap([startPoints[i][1]+meanVels[i][2]*vectorTime+.05], [startPoints[i][0]+meanVels[i][1]*vectorTime+.05])
		plt.plot(X,Y, color=thiscolor, lw=2, label='%s' % meanVels[i][0])
		plt.arrow(X[0], Y[0], (X[1]-X[0]), (Y[1]-Y[0]), lw=2, ec=thiscolor, fc=thiscolor, head_width=arrowHeadX[0]-X[1], head_length=arrowHeadY[0]-Y[1], )
	#print llrange
	plt.legend(loc='lower left', numpoints=1)
	#plt.show()
	
	return None
	

def getStationTraces(stations=['ucd1', 'ucsb'], startDate=dtm.datetime(2000, 1,1), endDate=dtm.datetime.now(), timeStep=dtm.timedelta(days=1)):
	lltype=1	# 0: meters and such, 1: wgsLat/Lon, 2: nadLat/Lon
	retList=[["time"]]
	for station in stations:
		retList[0]+=[station]
	#
	# get stations data:
	thisDt=startDate
	# timeStep=dtm.timedelta(days=1)
	#
	while thisDt<=endDate:
		thisYr=thisDt.year
		thisDOY=thisDt.timetuple()[7]
		retList+=[[thisDt.toordinal()]]
		#
		gpsData=getStationPositions(stations, thisYr, thisDOY)
		# returns:
		# columns key: 'site; x; y; z [ITRF2005]; wgsLat; wgsLon; wgsHt; nadLat; nadLon; nadHt; xSig; ySig; zSig; wgsLatSig; wgsLonSig; wgsHtSig; nadLatSig; nadLonSig; nadHtSig\n'
		# then a row, according to the key, for each site.
		for i in xrange(1, len(gpsData)):
			rowData=gpsData[i][:-1].split(';')	# note: drop terminating '\n'
			if lltype==0:
				XYdata=map(float, rowData[1:4]) + map(float, rowData[11:14])	# x,y,z, dx, dy, dz
			if lltype==1:
				XYdata=map(float, rowData[4:7]) + map(float, rowData[14:17])	# x,y,z, dx, dy, dz (wgs lat/lon)
				#print XYdata
			if lltype==2:
				XYdata=map(float, rowData[7:10]) + map(float, rowData[17:20])	# x,y,z, dx, dy, dz (nad lat/lon)
			retList[-1]+=[XYdata]
			
			
		# step forward:
		thisDt=thisDt+timeStep
	
	return retList

def plotStationTraces(stations=['ucd1', 'ucsb'], startDate=dtm.datetime(2000, 1,1), endDate=dtm.datetime.now(), deltaT=dtm.timedelta(days=1),verbose=False):
	traces=getStationTraces(stations, startDate, endDate, deltaT)
	plt.figure(0)
	plt.clf()
	#x0=-2628826.
	#y0=-4247932.
	x0=0
	y0=0
	#
	#tms=map(operator.itemgetter(0), traces))
	stations=traces[0]
	traces=traces[1:]
	#return traces
	numTraces=len(traces[0])-1
	#return traces
	for i in xrange(numTraces):
		xyz=map(operator.itemgetter(i+1), traces)
		if verbose: print "plot trace %d, len: %d, type: %s" % (i+1, len(xyz), type(xyz).__name__)
		if verbose: print map(operator.itemgetter(0), xyz)
		plt.plot(map(operator.itemgetter(0), xyz), map(operator.itemgetter(1),xyz), '.-', label=stations[i+1])
	plt.legend(loc='upper right')
	plt.show()
	#
	return traces
	#return stations
#
def plotDonsTriangle(stations=['farb', 'diab', 'ohln'], startDate=dtm.datetime(2000, 1,1), endDate=dtm.datetime.now(), deltaT=None):
	# deltaT=dtm.timedelta(days=1)
	plt.ion()
	timelen=10**8.0	# length magnifier for velocity vectors.
	timelenmap=10**0.0
	if deltaT==None:
		# instead of a proper trace, just get the first and last positions (start date, end date).
		# but, the easiest way to do this programatically is to use the getTraces() function.
		thistdelta=endDate-startDate
		traces=getStationTraces(stations, startDate, endDate, thistdelta)
	else:
		traces=getStationTraces(stations, startDate, endDate, deltaT)
		# and next, figure out if we ended on endDate. if so, add it (or just add it and delete it if it's the same).
		#if traces[-1][0]!=endDate.toordinal(): traces+=getStationTraces(stations, endDate, endDate)[1]
	#
	lfVelocities=getMeanStationVelocities1(traces)	# returns [[station, lon0, lat0, v_lon, vlat, siglon, siglat, sigVel_lon, sigVel_lat], []..]
	#
	#times=map(operator.itemgetter(0, traces))
	#stations=[]
	#for i in xrange(len(traces[0])):
	#	stations+=[map(operator.itemgetter(i), traces)]
	#
	polys=[]	# in general, we want polygons; one for each time-step of stations. these will be [[time, [[x,y], [x,y], ...]] ]
	for stas in traces:	# row of stations
		if stas[0]=='time':
			# header row
			continue
		# otherwise, we get [time, [position-sta1], [position-sta2], ...]
		polys+=[[stas[0], []]]	# new polygon; start with time. so a row is like [time, [positions]]
		#for sta in stas[1]:
		for i in xrange(1, len(stas)):
			sta=stas[i]
			# each station...
			#print polys, sta
			polys[-1][1]+=[[sta[1], sta[0]]]	# lon, lat 
		# and close the loop:
		polys[-1][1]+=[polys[-1][1][0]]
		#
	vels=[]
	for i in xrange(len(polys[0][1])-1):
		#vels+=[[polys[-1][1][i][0]-polys[0][1][i][0], polys[-1][1][i][1]-polys[0][1][i][1]] ]
		vels+=[[lfVelocities[i][3], lfVelocities[i][4]]] 
		#
	#
	# plotting:
	plt.figure(1)
	plt.clf()
	plt.figure(0)
	plt.clf()
	#
	for ply in polys:
	#while iply<len(polys):
		#ply=polys[iply]
		Xs=map(operator.itemgetter(0), ply[1])
		Ys=map(operator.itemgetter(1), ply[1])
		plt.fill(Xs, Ys, alpha=.25)
		#
	#
	for iply in xrange(len(polys[0][1])-1):
		thiscolor=yp.pycolor(iply)
		#print "iply: %d" % iply
		#plt.plot([Xs[iply], Xs[iply]+vels[iply][0]*timelen], [Ys[iply], Ys[iply]+vels[iply][1]*timelen], color=thiscolor)
		#plt.arrow(startPoints[i][1], startPoints[i][0], meanVels[i][2]*vectorTime, meanVels[i][1]*vectorTime, \
#lw=2, ec=thiscolor, fc=thiscolor, head_width=.05, head_length=.05 )
		plt.arrow(Xs[iply], Ys[iply], vels[iply][0]*timelen, vels[iply][1]*timelen, \
lw=2, ec=thiscolor, fc=thiscolor, head_width=.05, head_length=.05 )
	#
	# now the map plot:
	# set up a map:
	###############
	plt.figure(1)
	lats=Ys
	lons=Xs
	vlats=map(operator.itemgetter(1), vels)
	vlons=map(operator.itemgetter(0), vels)
	
	#lats+=(scipy.array(lats)+scipy.array(vlats)*timelen).tolist()
	#lons+=(scipy.array(lons)+scipy.array(vlons)*timelen).tolist()
	
	#print lats, lons
	llpad=.25
	llr=[[min(lats)-llpad, min(lons)-llpad], [max(lats)+llpad, max(lons)+llpad]]
	print llr
	cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
	catmap=mpbm.Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution ='l', projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
	#
	thisax=plt.gca()
	catmap.ax=thisax
	#
	catmap.drawcoastlines(color='gray')
	catmap.drawcountries(color='gray')
	catmap.drawstates(color='gray')
	catmap.drawrivers(color='gray')
	catmap.fillcontinents(color='beige')
	
	catmap.drawmeridians(range(int(llr[0][1]-2.0), int(llr[1][1]+2.0)), color='k', labels=[1,1,1,1])
	catmap.drawparallels(range(int(llr[0][0]-2.0), int(llr[1][0]+2.0)), color='k', labels=[1, 1, 1, 1])
	#
	iply=0	# just a polygon index tracker.
	for ply in polys:
	#while iply<len(polys):
		#ply=polys[iply]
		iply+=1
		Xs=map(operator.itemgetter(0), ply[1])
		Ys=map(operator.itemgetter(1), ply[1])
		mapXs, mapYs=catmap(Xs, Ys)
		
		if iply==1: plt.fill(mapXs, mapYs, alpha=.25)	# here, iply is the index of the polygon.
		#
	#
	# here, iply is the index of a single vertex of one (first or last) polygon.
	for iply in xrange(len(polys[-1][1])-1):
		thiscolor=yp.pycolor(iply)
		#print "iply: %d" % iply
		#plt.plot([Xs[iply], Xs[iply]+vels[iply][0]*timelen], [Ys[iply], Ys[iply]+vels[iply][1]*timelen], color=thiscolor)
		#plt.arrow(startPoints[i][1], startPoints[i][0], meanVels[i][2]*vectorTime, meanVels[i][1]*vectorTime, \
#lw=2, ec=thiscolor, fc=thiscolor, head_width=.05, head_length=.05 )
		mapX, mapY = catmap(Xs[iply], Ys[iply])
		#plt.arrow(Xs[iply], Ys[iply], vels[iply][0]*timelen, vels[iply][1]*timelen, lw=2, ec=thiscolor, fc=thiscolor, head_width=.05, head_length=.05 )
		#
		# and arrows for the map too. first, get the end positions and projected lengths.
		aryXs=scipy.array(Xs[iply])
		aryYs=scipy.array(Ys[iply])
		mapx1, mapy1 = catmap(aryXs+scipy.array(vels[iply][0])*timelen, aryYs + scipy.array(vels[iply][1])*timelen)
		mapVelX=scipy.array(mapx1)-mapX
		mapVelY=scipy.array(mapy1)-mapY	# projected lengths.
		# and get arrow head lengths:
		arrow_head_scaling=.03
		map_head_X, map_head_Y=catmap(aryXs+arrow_head_scaling, aryYs+arrow_head_scaling)
		map_head_w=map_head_X-scipy.array(mapX)
		map_head_l=map_head_Y-scipy.array(mapY)	# so the width/length might look weird if the projections are not square, since some "lengths" are x-like
																# and other "lengths" are y-like (and vice versa of course for "width"/y.
		#plt.arrow(mapX, mapY, vels[iply][0]*timelen, vels[iply][1]*timelen, lw=2, ec=thiscolor, fc=thiscolor, head_width=.05, head_length=.05 )
		#plt.arrow(mapX, mapY, mapVelX*timelenmap, mapVelY*timelenmap, lw=2, ec=thiscolor, fc=thiscolor, head_width=.05, head_length=.05 )
		plt.arrow(mapX, mapY, mapVelX*timelenmap, mapVelY*timelenmap, lw=2, ec=thiscolor, fc=thiscolor, head_width=map_head_w, head_length=map_head_l )
	#
	# and plot the last polygon as well:
	plt.fill(mapXs, mapYs, alpha=.25)
	###############
	pangles=getpolyangles(polys)
	plotpolyangles(pangles, 'asecs')
	#
	return [polys, vels]
	#
def getpolyangles(polys):
	# get angles between polygons.
	thetas=[]
	#basepoly=polys[0][1]
	#for i in xrange(len(basepoly)):
	#	basepoly[i]=scipy.array(basepoly[i])
		
	for rw in polys:
		# get vector pairs, [v0, v1], [v1, v2], ..., [vn, v0]
		# for each vector pair, get the angle;
		# return [[t, theta1, theta2, ..., thetaN]]
		thetas+=[[rw[0]]]
		nverts=len(rw[1])-1
		for i in xrange(len(rw[1])-1):
			v1=scipy.array(rw[1][(i+1)%nverts])-scipy.array(rw[1][i%nverts])
			v2=scipy.array(rw[1][(i+2)%nverts])-scipy.array(rw[1][(i+1)%nverts])
			thetas[-1]+=[getVectorAngle(v1, v2)]
		#
	return thetas

def plotpolyangles(polyangles, degs='rads'):
	# degs: degree measurements: 'rads', 'deg', 'amins', 'asecs'
	plt.figure(2)
	plt.clf()
	X=map(operator.itemgetter(0), polyangles)
	Ys=[]
	for i in xrange(1,len(polyangles[0])):
		Ys+=[map(operator.itemgetter(i), polyangles)]
		thisY=scipy.array(Ys[-1])-Ys[-1][0]
		# now the angle measurment. the native units are radians:
		# (there is a computationally faster way to do this, but it's less intuitive. we can
		# make a units index: ui={0,1,2,3} for amins, etc. and do a for-loop... but this is more straight-forward.
		if degs in ('deg', 'degs', 'amins', 'asecs'):
			thisY*=(360./(2.0*math.pi))
		if degs in ('amins', 'asecs'):
			thisY/=60.
		if degs in ('asecs'):
			thisY/=60.
		#
		#plt.plot(X,Ys[-1])
		plt.xlabel('time')
		plt.ylabel(degs)
		plt.plot(X,thisY)
	#
	return [X,Ys]
	
def getVectorAngle(u, v):
	u=scipy.array(u)
	v=scipy.array(v)
	c=dot(u,v)/norm(u)/norm(v)
	angle=arccos(c)
	return angle
	
#############
def printArrayToFile(ary, fname='arayData.htm'):
	fout=open(fname, 'w')
	for rw in ary:
		fout.write(rw)
		if rw[-1]!='\n': fout.write('\n')
	fout.close()
	return None
	
def getTextSiteString(siteCodes):
	siteString=''
	for elem in siteCodes:
		siteString=siteString + elem + ' '
	if siteString[-1]==' ': siteString=siteString[:-1]
	while '  ' in siteString: siteString.replace('  ', ' ')
	#
	return siteString

def getRefFrame(refFrame):
	# filter out some possible mistakes. return values are:
	# 'RERURN ITRF2005', 'ITRF2000 ONLY'
	rval=''
	if refFrame in ['RERURN ITRF2005', 'RERUN ITRF2005', '2005', 'ITRF2005']:
		rval = 'RERUN ITRF2005'
	elif refFrame in ['ITRF2000 ONLY', '2000', 'ITRF2000']:
		rval = 'ITRF2000 ONLY'
	else:
		rval = 'RERUN ITRF2005'
	#
	return rval
	
	
	
	
	
