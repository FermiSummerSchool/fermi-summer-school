import os,pyfits
from ROOT import TH1F,TH2F,TColor,gStyle,SetOwnership,TCanvas,TGraph
from numpy import array,pi,log10,zeros,arccos,cos
from math import floor
#from scipy.stats import rv_continuous as rvc

d2r=pi/180.
r2d=180./pi

gStyle.SetTextFont(132)
gStyle.SetTitleFont(132,'xyz')
gStyle.SetLabelFont(132,'xyz')
gStyle.SetCanvasColor(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetFrameFillColor(0)
gStyle.SetCanvasBorderMode(0)
gStyle.SetPadBorderSize(1)
gStyle.SetLegendBorderSize(1)
gStyle.SetOptStat(0)
gStyle.SetTitleBorderSize(1)
gStyle.SetMarkerSize(0.9)
gStyle.SetFuncWidth(1)

#function to put a 1D array with len1*len2 values into a len1Xlen2 dimension array
#arguments are:
#len1 (int) -- length of first dimension
#len2 (int) -- length of second dimension
#values (array) -- 1D array with len1*len2 values
def get2Darray(len1,len2,values):
	#make an value to cycle through the second indices
	#make it -1 since the very first thing we'll do is add 1 to it
	i=-1
	#make an array of zeros with appropriate dimensions and cycle through values
	myarray=zeros([len1,len2])
	for j,v in enumerate(values):
		#if we're at the start or we've gone through an integer multiple of len1, update the second index i
		i=(i+1 if j%len1==0 else i)
		myarray[j][i]=v
	return myarray

####################################################################################
####################################################################################
#############start functions for making effective area plots########################
####################################################################################
####################################################################################

#function to make a 2D histogram of the LATeffective area for a given set of IRFs, front, back, or front+back
#the arguments are:
#IRFs (str) -- insrument response functions to plot effective area from, include event class in name
#front (bool, optional) -- flag to plot for front converting events, default is True
#back (bool, optional) -- flag to plot for back converting events, default is False
##note, if both front and back are set to true, will plot the addition of the two effective areas
#effCorr (bool, optional) -- flag to apply livetime efficiency correction, default is False
#Fl (float, optional) -- value of livetime fraction for efficiency correction, acceptable values are from 0 (non-inclusive) to 1 (inclusive), default is 0.9
#CALDB (str, optional) -- path to CALDB directory, shouldn't need to supply this but use if IRFs not in default CALDB directory or if getting CALDB environment variable fails
#plot (bool, optional) -- flag to draw the plot, default is True
#color (bool, optional) -- flag to use color or grayscale, if False, default is True
def make2DAeffPlot(IRFs,front=True,back=False,effCorr=False,Fl=0.9,CALDB=None,plot=True,color=True):
	#try to get the CALDB path
	if CALDB==None:
		try:
			CALDB=os.getenv('CALDB',False)
		except:
			print "Error, CALDB path not supplied but could not get the CALDB environment variable"
			return None
	#check if we need to do both front and back
	if front and back:
		conv="Front + Back"
		#open effecitve area FITS files and get energy and costheta bin edges
		f_front_aeff=pyfits.open(CALDB+'/data/glast/lat/bcf/ea/aeff_%s_front.fits'%IRFs)
		f_back_aeff=pyfits.open(CALDB+'/data/glast/lat/bcf/ea/aeff_%s_back.fits'%IRFs)
		el=f_front_aeff[1].data.field('energ_lo')[0]
		eh=f_front_aeff[1].data.field('energ_hi')[0]
		ctl=f_front_aeff[1].data.field('ctheta_lo')[0]
		cth=f_front_aeff[1].data.field('ctheta_hi')[0]
		#bin up the effective area values
		front_aeff=f_front_aeff[1].data.field('effarea')[0]
		back_aeff=f_back_aeff[1].data.field('effarea')[0]
		#correct for efficiency
		if effCorr:
			if Fl<=1 and Fl>0:#make sure the livetime fraction is reasonable
				#get front efficiency parameters for piecewise fits
				front_corrVals=f_aeff[3].data.field('efficiency_pars')
				#get logarithmic midpoints of energy bins
				ecent=[10**((log10(E)+log10(e))/2.) for e,E in zip(el,eh)]
				front_c0vals=corrVals[0]
				front_c1vals=corrVals[1]
				#correct the effective area values
				for i in range(len(ctl)):
					for j,e in enumerate(ecent):
						front_aeff[i][j]*=(eff_c(log10(e),front_c0vals[0],front_c0vals[1],front_c0vals[2],front_c0vals[3],front_c0vals[4],front_c0vals[5])*Fl+eff_c(log10(e),front_c1vals[0],front_c1vals[1],front_c1vals[2],front_c1vals[3],front_c1vals[4],front_c1vals[5]))
				#do the same for the back
				back_corrVals=f_aeff[3].data.field('efficiency_pars')
				back_c0vals=corrVals[0]
				back_c1vals=corrVals[1]
				for i in range(len(ctl)):
					for j,e in enumerate(ecent):
						back_aeff[i][j]*=(eff_c(log10(e),back_c0vals[0],back_c0vals[1],back_c0vals[2],back_c0vals[3],back_c0vals[4],back_c0vals[5])*Fl+eff_c(log10(e),back_c1vals[0],back_c1vals[1],back_c1vals[2],back_c1vals[3],back_c1vals[4],back_c1vals[5]))
			else:
				print "effCorr set to True but livetime fraction =",Fl,'is not in (0,1]. No correction will be made.'
		#close the effective area FITS files
		f_front_aeff.close()
		f_back_aeff.close()
	else:
		#get which conversion type we're dealing with and open the corresponding effecitve area FITS file
		conv=('front' if front else 'back')
		f_aeff=pyfits.open(CALDB+'/data/glast/lat/bcf/ea/aeff_%s_%s.fits'%(IRFs,conv))
		#get the energy and costheta bin edges
		el=f_aeff[1].data.field('energ_lo')[0]
		eh=f_aeff[1].data.field('energ_hi')[0]
		ctl=f_aeff[1].data.field('ctheta_lo')[0]
		cth=f_aeff[1].data.field('ctheta_hi')[0]
		#bin up the effective area values appropriately
		#aeff=get2Darray(len(el),len(ctl),f_aeff[1].data.field('effarea'))
		aeff=f_aeff[1].data.field('effarea')[0]
		if effCorr:
			if Fl<=1 and Fl>0:#make sure the livetime fraction is reasonable
				#get efficiency paramters for piecwise function
				corrVals=f_aeff[3].data.field('efficiency_pars')
				#get logarithmic midpoints of energy bins
				ecent=[10**((log10(E)+log10(e))/2.) for e,E in zip(el,eh)]
				c0vals=corrVals[0]
				c1vals=corrVals[1]
				#correct effective area values
				for i in range(len(ctl)):
					for j,e in enumerate(ecent):
						aeff[i][j]*=(eff_c(log10(e),c0vals[0],c0vals[1],c0vals[2],c0vals[3],c0vals[4],c0vals[5])*Fl+eff_c(log10(e),c1vals[0],c1vals[1],c1vals[2],c1vals[3],c1vals[4],c1vals[5]))
			else:
				print "effCorr set to True but livetime fraction =",Fl,'is not in (0,1]. No correction will be made.'
		#close the effective area FITS file
		f_aeff.close()
	#make arrays for energy and costheta bins
	ebins=array([e for e in el]+[eh[-1]])
	ctbins=array([c for c in ctl]+[cth[-1]])
	#make 2D histogram
	aehist=TH2F('aehist','',len(el),ebins,len(ctl),ctbins)
	#fill histogram
	for i in range(len(el)):
		for j in range(len(ctl)):
			if front and back:#if front and back, add the two effective areas
				aehist.SetBinContent(i+1,j+1,front_aeff[j][i]+back_aeff[j][i])
			else:
				aehist.SetBinContent(i+1,j+1,aeff[j][i])
	#set axis titles
	aehist.GetXaxis().CenterTitle()
	aehist.GetYaxis().CenterTitle()
	aehist.SetXTitle('Energy (MeV)')
	aehist.SetYTitle('cos(#theta)')
	aehist.SetContour(99)
	if plot:
		#if plotting, take care of some formatting issues
		if color:
			BPalette()
		else:
			GPalette()
		aeffcan=TCanvas('aeffcan','%s %s Effective Area'%(IRFs,conv),1000,800)
		aeffcan.SetTicks(1,1)
		aeffcan.SetLogx(1)
		aeffcan.SetFrameFillColor(1)
		aeffcan.SetFrameFillStyle(1001)
		aehist.GetXaxis().SetAxisColor(0)
		aehist.GetYaxis().SetAxisColor(0)
		aehist.Draw('cont4z')
		SetOwnership(aeffcan,False)
		SetOwnership(aehist,False)
	#return the 2D hist
	return aehist

#function to make a 1D histogram of the LATeffective area vs theta or energy for a given set of IRFs, front, back, or front+back
#the arguments are:
#IRFs (str) -- insrument response functions to plot effective area from, include event class in name
#ptype (str) -- either "energy" or "theta", tells the function what to plot Aeff vs.
#value (float) -- either the energy, in MeV, or the costheta value for which you want to plot Aeff, if ptype is "energy" this is the costheta value, if ptype is "costheta" this is the energy value
#NOTE, avoid known bin edges, the values in the IRFs FITS files have weird extra digits 8 or 9 places after the decimal point
#npts (int,optional) -- number of points to plot, if not specified will plot for the number of bins in the effective area FITS file, note that this is more of a minimum number of points, if the requested value is not an integer divisor of the number of bins in the FITS file you will have at least this many points
#front (bool, optional) -- flag to plot for front converting events, default is True
#back (bool, optional) -- flag to plot for back converting events, default is False
##note, if both front and back are set to true, will plot the addition of the two effective areas
#effCorr (bool, optional) -- flag to apply livetime efficiency correction, default is False
#Fl (float, optional) -- value of livetime fraction for efficiency correction, acceptable values are from 0 (non-inclusive) to 1 (inclusive), default is 0.9
#CALDB (str, optional) -- path to CALDB directory, shouldn't need to supply this but use if IRFs not in default CALDB directory or if getting CALDB environment variable fails
#plot (bool, optional) -- flag to draw the plot, default is True
def make1DAeffPlot(IRFs,ptype,value,npts=None,front=True,back=False,effCorr=False,Fl=0.9,CALDB=None,plot=True):
	#try to get the CALDB path
	if CALDB==None:
		try:
			CALDB=os.getenv('CALDB',False)
		except:
			print "Error, CALDB path not supplied but could not get the CALDB environment variable"
			return None
	#check ptype value
	if ptype not in ['energy','theta']:
		print "Error, input for ptype of",ptype,'is invalid, only acceptable inputs are "energy" or "costheta".'
		return None
	#check if we need to do both front and back
	if front and back:
		conv="Front + Back"
		#open effecitve area FITS files and get energy and costheta bin edges
		f_front_aeff=pyfits.open(CALDB+'/data/glast/lat/bcf/ea/aeff_%s_front.fits'%IRFs)
		f_back_aeff=pyfits.open(CALDB+'/data/glast/lat/bcf/ea/aeff_%s_back.fits'%IRFs)
		el=f_front_aeff[1].data.field('energ_lo')[0]
		eh=f_front_aeff[1].data.field('energ_hi')[0]
		ctl=f_front_aeff[1].data.field('ctheta_lo')[0]
		cth=f_front_aeff[1].data.field('ctheta_hi')[0]
		#bin up the effective area values
		front_aeff=f_front_aeff[1].data.field('effarea')[0]
		back_aeff=f_back_aeff[1].data.field('effarea')[0]
		#correct for efficiency
		if effCorr:
			if Fl<=1 and Fl>0:#make sure the livetime fraction is reasonable
				#get front efficiency parameters for piecewise fits
				front_corrVals=f_aeff[3].data.field('efficiency_pars')
				#get logarithmic midpoints of energy bins
				ecent=[10**((log10(E)+log10(e))/2.) for e,E in zip(el,eh)]
				front_c0vals=corrVals[0]
				front_c1vals=corrVals[1]
				#correct the effective area values
				for i in range(len(ctl)):
					for j,e in enumerate(ecent):
						front_aeff[i][j]*=(eff_c(log10(e),front_c0vals[0],front_c0vals[1],front_c0vals[2],front_c0vals[3],front_c0vals[4],front_c0vals[5])*Fl+eff_c(log10(e),front_c1vals[0],front_c1vals[1],front_c1vals[2],front_c1vals[3],front_c1vals[4],front_c1vals[5]))
				#do the same for the back
				back_corrVals=f_aeff[3].data.field('efficiency_pars')
				back_c0vals=corrVals[0]
				back_c1vals=corrVals[1]
				for i in range(len(ctl)):
					for j,e in enumerate(ecent):
						back_aeff[i][j]*=(eff_c(log10(e),back_c0vals[0],back_c0vals[1],back_c0vals[2],back_c0vals[3],back_c0vals[4],back_c0vals[5])*Fl+eff_c(log10(e),back_c1vals[0],back_c1vals[1],back_c1vals[2],back_c1vals[3],back_c1vals[4],back_c1vals[5]))
			else:
				print "effCorr set to True but livetime fraction =",Fl,'is not in (0,1]. No correction will be made.'
		if ptype=='energy':
			#need to find correct costheta bin
			ctidx=-1
			for j in range(len(ctl)):
				if value>=ctl[j] and value<cth[j]:
					ctidx=j
					break
				if j==len(ctl)-1 and ctidx==-1:
					ctidx=len(ctl)-1
			print 'Using Aeff values from bin with costheta [%f,%f].'%(ctl[ctidx],cth[ctidx])
			#now check npts value
			npts=(len(el) if npts==None else npts)
			if npts>len(el):
				print "Requested number of points (%i) is greater than number of energy bins in IRFs (%i), will produce plot with %i points."%(npts,len(el),len(el))
				npts=len(el)
			#average the effective area values down into the desired number of points (if possible) from the corresponding costheta bin
			step=int(len(el)/npts)
			if step>1:
				aeff_new=[(sum([front_aeff[ctidx][i+step*j] for i in range(step)])+sum([back_aeff[ctidx][i+step*j] for i in range(step)]))/float(step) for j in range(len(el)/step)]
				bins=array([el[i*step] for i in range(len(el)/step)]+[eh[-1]])
			else:#if step is one then we're not doing any rebinning, just need to get the values for the desired costheta
				aeff_new=[front_aeff[ctidx][i]+back_aeff[ctidx][i] for i in range(len(el))]
				bins=array([e for e in el]+[eh[-1]])
		else:
			#same as above but picking one energy and looking at the effective area for different costheta values
			#need to find the correct energy bin
			eidx=-1
			for i in range(len(el)):
				if value>=el[i] and value<eh[i]:
					eidx=i
					break
				if i==len(el)-1 and eidx==-1:
					eidx=len(el)-1
			print 'Using Aeff values from bin with energy [%f,%f] MeV.'%(el[eidx],eh[eidx])
			npts=(len(ctl) if npts==None else npts)
			if npts>len(ctl):
				print "Requested number of points (%i) is greater than number of costheta bins in IRFs (%i), will produce plot with %i points."%(npts,len(ctl),len(ctl))
				npts=len(ctl)
			step=int(len(ctl)/npts)
			if step>1:
				aeff_new=[(sum([front_aeff[i+step*j][eidx] for i in range(step)])+sum([back_aeff[i+step*j][eidx] for i in range(step)]))/float(step) for j in range(len(ctl)/step)]
				aeff_new=[aeff_new[len(aeff_new)-1-i] for i in range(len(aeff_new))]
				bins=array([ctl[i*step] for i in range(len(ctl)/step)]+[cth[-1]])
				bins=array([r2d*arccos(bins[len(bins)-1-i]) for i in range(len(bins))])
			else:
				aeff_new=[front_aeff[i][eidx]+back_aeff[i][eidx] for i in range(len(ctl))]
				aeff_new=[aeff_new[len(ctl)-1-i] for i in range(len(ctl))]
				bins=array([c for c in ctl]+[cth[-1]])
				bins=array([r2d*arccos(bins[len(bins)-1-i]) for i in range(len(bins))])
		#close the effective area FITS files
		f_front_aeff.close()
		f_back_aeff.close()
	else:
		#get which conversion type we're dealing with and open the corresponding effecitve area FITS file
		conv=('front' if front else 'back')
		f_aeff=pyfits.open(CALDB+'/data/glast/lat/bcf/ea/aeff_%s_%s.fits'%(IRFs,conv))
		#get the energy and costheta bin edges
		el=f_aeff[1].data.field('energ_lo')[0]
		eh=f_aeff[1].data.field('energ_hi')[0]
		ctl=f_aeff[1].data.field('ctheta_lo')[0]
		cth=f_aeff[1].data.field('ctheta_hi')[0]
		#bin up the effective area values appropriately
		aeff=f_aeff[1].data.field('effarea')[0]
		if effCorr:
			if Fl<=1 and Fl>0:#make sure the livetime fraction is reasonable
				#get efficiency paramters for piecwise function
				corrVals=f_aeff[3].data.field('efficiency_pars')
				#get logarithmic midpoints of energy bins
				ecent=[10**((log10(E)+log10(e))/2.) for e,E in zip(el,eh)]
				c0vals=corrVals[0]
				c1vals=corrVals[1]
				#correct effective area values
				for i in range(len(ctl)):
					for j,e in enumerate(ecent):
						aeff[i][j]*=(eff_c(log10(e),c0vals[0],c0vals[1],c0vals[2],c0vals[3],c0vals[4],c0vals[5])*Fl+eff_c(log10(e),c1vals[0],c1vals[1],c1vals[2],c1vals[3],c1vals[4],c1vals[5]))
			else:
				print "effCorr set to True but livetime fraction =",Fl,'is not in (0,1]. No correction will be made.'
		if ptype=='energy':
			#need to find correct costheta bin
			ctidx=-1
			for j in range(len(ctl)):
				if value>=ctl[j] and value<cth[j]:
					ctidx=j
					break
				if j==len(ctl)-1 and ctidx==-1:
					ctidx=len(ctl)-1
			print 'Using Aeff values from bin with costheta [%f,%f].'%(ctl[ctidx],cth[ctidx])
			#check the requested number of points
			npts=(len(el) if npts==None else npts)
			if npts>len(el):
				print "Requested number of points (%i) is greater than number of energy bins in IRFs (%i), will produce plot with %i points."%(npts,len(el),len(el))
				npts=len(el)
			#now average the effective area values for some energies and a particular costheta value
			step=int(len(el)/npts)
			if step>1:
				aeff_new=[sum([aeff[ctidx][i+step*j] for i in range(step)])/float(step) for j in range(len(el)/step)]
				bins=array([el[i*step] for i in range(len(el)/step)]+[eh[-1]])
			else:#if step==1 then we're just grabbing the effective area values at different energies for a specific costheta
				aeff_new=[aeff[ctidx][i] for i in range(len(el))]
				bins=array([e for e in el]+[eh[-1]])
		else:
			#now do the same but this time for a particular energy and varying costheta
			#need to find the correct energy bin
			eidx=-1
			for i in range(len(el)):
				if value>=el[i] and value<eh[i]:
					eidx=i
					break
				if i==len(el)-1 and eidx==-1:
					eidx=len(el)-1
			print 'Using Aeff values from bin with energy [%f,%f] MeV.'%(el[eidx],eh[eidx])
			npts=(len(ctl) if npts==None else npts)
			if npts>len(ctl):
				print "Requested number of points (%i) is greater than number of costheta bins in IRFs (%i), will produce plot with %i points."%(npts,len(ctl),len(ctl))
				npts=len(ctl)
			step=int(len(ctl)/npts)
			if step>1:
				aeff_new=[sum([aeff[i+step*j][eidx] for i in range(step)])/float(step) for j in range(len(ctl)/step)]
				aeff_new=[aeff_new[len(aeff_new)-1-i] for i in range(len(aeff_new))]
				bins=array([ctl[i*step] for i in range(len(ctl)/step)]+[cth[-1]])
				bins=array([r2d*arccos(bins[len(bins)-1-i]) for i in range(len(bins))])
			else:
				aeff_new=[aeff[len(ctl)-1-i][eidx] for i in range(len(ctl))]
				bins=array([c for c in ctl]+[cth[-1]])
				bins=array([r2d*arccos(bins[len(bins)-1-i]) for i in range(len(bins))])
		#close the effective area FITS file
		f_aeff.close()
	#make 1D histogram
	hist=TH1F('hist','',len(aeff_new),bins)
	#fill histogram
	for i,a in enumerate(aeff_new):
		hist.SetBinContent(i+1,a)
	#set axis titles and do other style stuff
	hist.GetXaxis().CenterTitle()
	hist.GetYaxis().CenterTitle()
	xtitle=('Energy (MeV)' if ptype=='energy' else 'cos(#theta)')
	hist.SetXTitle(xtitle)
	hist.SetYTitle('A_{eff} (m^{2} )')
	hist.SetLineColor(1)
	hist.SetMarkerColor(1)
	hist.SetMarkerStyle(20)
	if plot:
		can=TCanvas('can','%s %s Effective Area vs. %s'%(IRFs,conv,ptype),1000,800)
		can.SetTicks(1,1)
		if ptype=='energy':
			can.SetLogx(1)
		hist.Draw('pl')
		SetOwnership(can,False)
		SetOwnership(hist,False)
	#return the 1D hist
	return hist

#function to make a 1D histogram of the LATeffective area for a given set of IRFs, front, back, or front+back versus instrument azimuth for a particular energy and costheta bin
#the arguments are:
#IRFs (str) -- insrument response functions to plot effective area from, include event class in name
#energy (float) -- energy in MeV for effective area value to apply phi correction to
#costheta (float) -- costheta for effective are value to apply phi correction to
#nphi (int, optional) -- number of points from 0 to 90 degrees for plotting (plot will be drawn from 0 to 360 degrees), default is 12
#front (bool, optional) -- flag to plot for front converting events, default is True
#back (bool, optional) -- flag to plot for back converting events, default is False
##note, if both front and back are set to true, will plot the addition of the two effective areas
#effCorr (bool, optional) -- flag to apply livetime efficiency correction, default is False
#Fl (float, optional) -- value of livetime fraction for efficiency correction, acceptable values are from 0 (non-inclusive) to 1 (inclusive), default is 0.9
#CALDB (str, optional) -- path to CALDB directory, shouldn't need to supply this but use if IRFs not in default CALDB directory or if getting CALDB environment variable fails
#plot (bool, optional) -- flag to draw the plot, default is True
def makephiAeffPlot(IRFs,energy,costheta,nphi=12,front=True,back=False,effCorr=False,Fl=0.9,CALDB=None,plot=True):
	#try to get the CALDB path
	if CALDB==None:
		try:
			CALDB=os.getenv('CALDB',False)
		except:
			print "Error, CALDB path not supplied but could not get the CALDB environment variable"
			return None
	#check if we need to do both front and back
	if front and back:
		conv='Front + Back'
		#open effecitve area FITS files and get energy and costheta bin edges
		f_front_aeff=pyfits.open(CALDB+'/data/glast/lat/bcf/ea/aeff_%s_front.fits'%IRFs)
		f_back_aeff=pyfits.open(CALDB+'/data/glast/lat/bcf/ea/aeff_%s_back.fits'%IRFs)
		el=f_front_aeff[1].data.field('energ_lo')[0]
		eh=f_front_aeff[1].data.field('energ_hi')[0]
		ctl=f_front_aeff[1].data.field('ctheta_lo')[0]
		cth=f_front_aeff[1].data.field('ctheta_hi')[0]
		#bin up the effective area values
		front_aeff=f_front_aeff[1].data.field('effarea')[0]
		back_aeff=f_back_aeff[1].data.field('effarea')[0]
		#correct for efficiency
		if effCorr:
			if Fl<=1 and Fl>0:#make sure the livetime fraction is reasonable
				#get front efficiency parameters for piecewise fits
				front_corrVals=f_aeff[3].data.field('efficiency_pars')
				#get logarithmic midpoints of energy bins
				ecent=[10**((log10(E)+log10(e))/2.) for e,E in zip(el,eh)]
				front_c0vals=corrVals[0]
				front_c1vals=corrVals[1]
				#correct the effective area values
				for i in range(len(ctl)):
					for j,e in enumerate(ecent):
						front_aeff[i][j]*=(eff_c(log10(e),front_c0vals[0],front_c0vals[1],front_c0vals[2],front_c0vals[3],front_c0vals[4],front_c0vals[5])*Fl+eff_c(log10(e),front_c1vals[0],front_c1vals[1],front_c1vals[2],front_c1vals[3],front_c1vals[4],front_c1vals[5]))
				#do the same for the back
				back_corrVals=f_aeff[3].data.field('efficiency_pars')
				back_c0vals=corrVals[0]
				back_c1vals=corrVals[1]
				for i in range(len(ctl)):
					for j,e in enumerate(ecent):
						back_aeff[i][j]*=(eff_c(log10(e),back_c0vals[0],back_c0vals[1],back_c0vals[2],back_c0vals[3],back_c0vals[4],back_c0vals[5])*Fl+eff_c(log10(e),back_c1vals[0],back_c1vals[1],back_c1vals[2],back_c1vals[3],back_c1vals[4],back_c1vals[5]))
			else:
				print "effCorr set to True but livetime fraction =",Fl,'is not in (0,1]. No correction will be made.'
		#now get the phi-dependence parameters
		#need to get the energy and costheta bin edges again as they have different dimensions than in extension 1
		pel=f_front_aeff[2].data.field('energ_lo')[0]
		peh=f_front_aeff[2].data.field('energ_hi')[0]
		pctl=f_front_aeff[2].data.field('ctheta_lo')[0]
		pcth=f_front_aeff[2].data.field('ctheta_hi')[0]
		#get the phi-dependence parameters binned up, these are q0 and q1 from equation 16 in section 5.2.3 of the pass7 performance paper
		front_q0vals=f_front_aeff[2].data.field('phidep0')[0]
		front_q1vals=f_front_aeff[2].data.field('phidep1')[0]
		back_q0vals=f_back_aeff[2].data.field('phidep0')[0]
		back_q1vals=f_back_aeff[2].data.field('phidep1')[0]
		#figure out which energy and costheta bins we want to grab the correct q0 and q1
		ctidx=-1
		eidx=-1
		for i in range(len(pel)):
			if energy>=pel[i] and energy<peh[i]:
				eidx=i
				break
			if i==len(pel)-1 and eidx==-1:
				eidx=len(pel)-1
		for j in range(len(pctl)-1):
			if costheta>=pctl[j] and costheta<pcth[j]:
				ctidx=j
				break
			if j==len(pctl) and ctidx==-1:
				ctidx=len(pctl)-1
		print 'Using q0 and q1 values from bin with energy [%f,%f] MeV and costheta [%f,%f]'%(pel[eidx],peh[eidx],pctl[ctidx],pcth[ctidx])
		#get the phis and then map to the xi parameter
		#just do for 0 to 90 degrees, same values for other 4 quadrants
		phis=[0+i*(90./(nphi-1.)) for i in range(nphi)]
		xis=[xi(p*d2r) for p in phis]
		#get the correction for each xi
		front_f_xi=[1.+front_q0vals[ctidx][eidx]*x**front_q1vals[ctidx][eidx] for x in xis]
		back_f_xi=[1.+back_q0vals[ctidx][eidx]*x**back_q1vals[ctidx][eidx] for x in xis]
		#find the energy and costheta bin for the effective area to be modified
		ctidx=-1
		eidx=-1
		for i in range(len(el)):
			if energy>=el[i] and energy<eh[i]:
				eidx=i
				break
			if i==len(el)-1 and eidx==-1:
				eidx=len(el)-1
		for j in range(len(ctl)-1):
			if costheta>=ctl[j] and costheta<cth[j]:
				ctidx=j
				break
			if j==len(ctl) and ctidx==-1:
				ctidx=len(ctl)-1
		print 'Using Aeff value from bin with energy [%f,%f] MeV and costheta [%f,%f]'%(el[eidx],eh[eidx],ctl[ctidx],cth[ctidx])
		#apply the correction for each xi
		scaledAeff=[front_aeff[ctidx][eidx]*ff+back_aeff[ctidx][eidx]*bf for ff,bf in zip(front_f_xi,back_f_xi)]
		#close the effective area FITS files
		f_front_aeff.close()
		f_back_aeff.close()
	else:
		#get which conversion type we're dealing with and open the corresponding effecitve area FITS file
		conv=('front' if front else 'back')
		f_aeff=pyfits.open(CALDB+'/data/glast/lat/bcf/ea/aeff_%s_%s.fits'%(IRFs,conv))
		#get the energy and costheta bin edges
		el=f_aeff[1].data.field('energ_lo')[0]
		eh=f_aeff[1].data.field('energ_hi')[0]
		ctl=f_aeff[1].data.field('ctheta_lo')[0]
		cth=f_aeff[1].data.field('ctheta_hi')[0]
		#bin up the effective area values appropriately
		#aeff=get2Darray(len(el),len(ctl),f_aeff[1].data.field('effarea'))
		aeff=f_aeff[1].data.field('effarea')[0]
		if effCorr:
			if Fl<=1 and Fl>0:#make sure the livetime fraction is reasonable
				#get efficiency paramters for piecwise function
				corrVals=f_aeff[3].data.field('efficiency_pars')
				#get logarithmic midpoints of energy bins
				ecent=[10**((log10(E)+log10(e))/2.) for e,E in zip(el,eh)]
				c0vals=corrVals[0]
				c1vals=corrVals[1]
				#correct effective area values
				for i in range(len(ctl)):
					for j,e in enumerate(ecent):
						aeff[j][i]*=(eff_c(log10(e),c0vals[0],c0vals[1],c0vals[2],c0vals[3],c0vals[4],c0vals[5])*Fl+eff_c(log10(e),c1vals[0],c1vals[1],c1vals[2],c1vals[3],c1vals[4],c1vals[5]))
			else:
				print "effCorr set to True but livetime fraction =",Fl,'is not in (0,1]. No correction will be made.'
		#now get the phi-dependence parameters
		#need to get the energy and costheta bin edges again as they have different dimensions than in extension 1
		pel=f_aeff[2].data.field('energ_lo')[0]
		peh=f_aeff[2].data.field('energ_hi')[0]
		pctl=f_aeff[2].data.field('ctheta_lo')[0]
		pcth=f_aeff[2].data.field('ctheta_hi')[0]
		#get the phi-dependence parameters binned up, these are q0 and q1 from equation 16 in section 5.2.3 of the pass7 performance paper
		q0vals=f_aeff[2].data.field('phidep0')[0]
		q1vals=f_aeff[2].data.field('phidep1')[0]
		#figure out which energy and costheta bins we want to grab the correct q0 and q1
		ctidx=-1
		eidx=-1
		for i in range(len(pel)):
			if energy>=pel[i] and energy<peh[i]:
				eidx=i
				break
			if i==len(pel)-1 and eidx==-1:
				eidx=len(pel)-1
		for j in range(len(pctl)-1):
			if costheta>=pctl[j] and costheta<pcth[j]:
				ctidx=j
				break
			if j==len(pctl) and ctidx==-1:
				ctidx=len(pctl)-1
		print 'Using q0 and q1 values from bin with energy [%f,%f] MeV and costheta [%f,%f]'%(pel[eidx],peh[eidx],pctl[ctidx],pcth[ctidx])
		#get the phis and then map to the xi parameter
		#just do for 0 to 90 degrees, same values for other 4 quadrants
		phis=[0+i*(90./(nphi-1.)) for i in range(nphi)]
		xis=[xi(p*d2r) for p in phis]
		#get the correction for each xi
		f_xi=[1.+q0vals[ctidx][eidx]*x**q1vals[ctidx][eidx] for x in xis]
		#find the energy and costheta bin for the effective area to be modified
		ctidx=-1
		eidx=-1
		for i in range(len(el)):
			if energy>=el[i] and energy<eh[i]:
				eidx=i
				break
			if i==len(el)-1 and eidx==-1:
				eidx=len(el)-1
		for j in range(len(ctl)-1):
			if costheta>=ctl[j] and costheta<cth[j]:
				ctidx=j
				break
			if j==len(ctl) and ctidx==-1:
				ctidx=len(ctl)-1
		print 'Using Aeff value from bin with energy [%f,%f] MeV and costheta [%f,%f]'%(el[eidx],eh[eidx],ctl[ctidx],cth[ctidx])
		#apply the correction for each xi
		scaledAeff=[aeff[ctidx][eidx]*f for f in f_xi]
		#close the effective area FITS file
		f_aeff.close()
	#make phi bins from 0 to 360 degrees and build and duplicate the effective area values accordingly
	phibins=array(phis+[p+90. for p in phis[1:]]+[p+180. for p in phis[1:]]+[p+270. for p in phis[1:]])
	AeffValues=array(scaledAeff+scaledAeff[1:]+scaledAeff[1:]+scaledAeff[1:])
	#make a graph of the scaled effective area vs. phi
	phigraph=TGraph(len(phibins),phibins,AeffValues)
	phigraph.SetMarkerStyle(20)
	phigraph.SetMarkerColor(1)
	phigraph.SetLineColor(1)
	if plot:
		#if plot set to True, draw and make some style choices
		pcanv=TCanvas('pcanv','%s %s A_{eff} vs. phi for E = %.1f MeV and theta = %.1f deg'%(IRFs,conv,energy,r2d*arccos(costheta)),1000,800)
		pcanv.SetTicks(1,1)
		pdummy=TH1F('pdummy','',1000,0,360)
		pdummy.SetXTitle('#phi (#circ)')
		pdummy.GetXaxis().CenterTitle()
		pdummy.SetYTitle('A_{eff} (m^{2} )')
		pdummy.GetYaxis().CenterTitle()
		pdummy.Draw()
		phigraph.Draw('plsame')
		SetOwnership(pcanv,False)
		SetOwnership(phigraph,False)
		SetOwnership(pdummy,False)
	return phigraph

#function to calculate the c0 and c1 constants for the effective area efficiency correction from the piecwise fit values vs. log energy
#see section 5.2.2 of the pass7 performance paper
#arguments are:
#logene (float) -- log10 of the energy at which the desired constant should be calculated
#a0 (float) -- first slope of the piecewise linear fit
#b0 (float) -- first intercept of the piecewise linear fit
#a1 (float) -- second slope of the piecewise linear fit
#logE1 (float) -- log10 of the first break energy of the piecewise linear fit
#a2 (float) -- third slope of the piecewise linear fit
#logE2 (float) -- log10 of the second break energy of the piecewise linear fit
def eff_c(logene,a0,b0,a1,logE1,a2,logE2):
	if logene<logE1:
		c=b0+a0*logene
	elif logene<logE2:
		c=b0+(a0-a1)*logE1+a2*logene
	else:
		c=b0+(a0-a1)*logE1+(a1-a2)*logE2+a2*logene
	return c

#function to map phi angles to the interval [0,1] as give in Eq. 15 from section 5.2.3 in pass7 performance paper
#arguments are:
#phi (flot) -- phi angle in radians
def xi(phi):
	return (4./pi)*abs((phi%(pi/2.))-(pi/4.))

####################################################################################
####################################################################################
#############end functions for making effective area plots##########################
####################################################################################
####################################################################################


####################################################################################
####################################################################################
#############start functions for making psf plots###################################
####################################################################################
####################################################################################

#function to make a 1D histogram of the LAT PSF containment angle (for a specified percent containment) or ratio of containment radii vs theta or energy for a given set of IRFs, front, back, or front+back
#the arguments are:
#IRFs (str) -- insrument response functions to plot effective area from, include event class in name
#perc (float) -- percentage containment desired, expressed as a fraction (e.g., 0.68 not 68 for 68% containment)
#ptype (str) -- either "energy" or "theta", tells the function what to plot Aeff vs.
#value (float) -- either the energy, in MeV, or the costheta value for which you want to plot Aeff, if ptype is "energy" this is the costheta value, if ptype is "costheta" this is the energy value
#npts (int,optional) -- number of points to plot, if not specified will plot for the number of bins in the effective area FITS file, note that this is more of a minimum number of points, if the requested value is not an integer divisor of the number of bins in the FITS file you will have at least this many points
#perc2 (float, optional) -- if specified, this tells the function to make a containment angle ratio plot instead, it will plot the ratio of (perc2 containment angle)/(perc containment angle)
#front (bool, optional) -- flag to plot for front converting events, default is True
#back (bool, optional) -- flag to plot for back converting events, default is False
##note, if both front and back are set to true, will plot the average of the two containment radii????
##should be some sort of weighted average but I can't recall the ratios...based on expected percentages of each event type
#CALDB (str, optional) -- path to CALDB directory, shouldn't need to supply this but use if IRFs not in default CALDB directory or if getting CALDB environment variable fails
#plot (bool, optional) -- flag to draw the plot, default is True
def make1DPSFPlot(IRFs,perc,ptype,value,npts=None,perc2=None,front=True,back=False,CALDB=None,plot=True):
	#try to get the CALDB path
	if CALDB==None:
		try:
			CALDB=os.getenv('CALDB',False)
		except:
			print "Error, CALDB path not supplied but could not get the CALDB environment variable"
			return None
	#check ptype value
	if ptype not in ['energy','theta']:
		print "Error, input for ptype of",ptype,'is invalid, only acceptable inputs are "energy" or "theta".'
		return None
	#check if we need to do both front and back
	if front and back:
		print "Calculating containment ratios for FRONT+BACK as the average"
		print "this is approximately correct near 685","containment but underestimates for larger containment percentages"
		conv="Front + Back"
		f_front_psf=pyfits.open(CALDB+'/data/glast/lat/bcf/psf/psf_%s_front.fits'%IRFs)
		f_back_psf=pyfits.open(CALDB+'/data/glast/lat/bcf/psf/psf_%s_back.fits'%IRFs)
		#get the energy and costheta bin edges
		el=f_front_psf[1].data.field('energ_lo')[0]
		eh=f_front_psf[1].data.field('energ_hi')[0]
		ctl=f_front_psf[1].data.field('ctheta_lo')[0]
		cth=f_front_psf[1].data.field('ctheta_hi')[0]
		#get front values
		fNTAIL=f_front_psf[1].data.field('NTAIL')[0]
		fSCORE=f_front_psf[1].data.field('SCORE')[0]
		fSTAIL=f_front_psf[1].data.field('STAIL')[0]
		fGCORE=f_front_psf[1].data.field('GCORE')[0]
		fGTAIL=f_front_psf[1].data.field('GTAIL')[0]
		#get back values
		bNTAIL=f_back_psf[1].data.field('NTAIL')[0]
		bSCORE=f_back_psf[1].data.field('SCORE')[0]
		bSTAIL=f_back_psf[1].data.field('STAIL')[0]
		bGCORE=f_back_psf[1].data.field('GCORE')[0]
		bGTAIL=f_back_psf[1].data.field('GTAIL')[0]
		psfpars=f_back_psf[2].data.field('PSFSCALE')
		fc0=psfpars[0][0]
		fc1=psfpars[0][1]
		bc0=psfpars[0][2]
		bc1=psfpars[0][3]
		beta=psfpars[0][4]
		if ptype=='energy':
			#need to find correct costheta bin
			ctidx=-1
			for j in range(len(ctl)):
				if value>=ctl[j] and value<cth[j]:
					ctidx=j
					break
				if j==len(ctl)-1 and ctidx==-1:
					ctidx=len(ctl)-1
			print 'Using PSF values from bin with costheta [%f,%f].'%(ctl[ctidx],cth[ctidx])
			#now get the requisite perecentage containment radii
			ecent=[10**((log10(E)+log10(e))/2.) for e,E in zip(el,eh)]
			fContRad=[getPercContRadius(perc,ecent[i],fc0,fc1,beta,fNTAIL[ctidx][i],fSCORE[ctidx][i],fSTAIL[ctidx][i],fGCORE[ctidx][i],fGTAIL[ctidx][i]) for i in range(len(ecent))]
			bContRad=[getPercContRadius(perc,ecent[i],bc0,bc1,beta,bNTAIL[ctidx][i],bSCORE[ctidx][i],bSTAIL[ctidx][i],bGCORE[ctidx][i],bGTAIL[ctidx][i]) for i in range(len(ecent))]
			#for now we'll just assume they contribute equally, i.e., roughly equal number of front and back events
			ContRad=[(B+F)/2. for F,B in zip(fContRad,bContRad)]
			if perc2!=None:
				fContRad2=[getPercContRadius(perc2,ecent[i],fc0,fc1,beta,fNTAIL[ctidx][i],fSCORE[ctidx][i],fSTAIL[ctidx][i],fGCORE[ctidx][i],fGTAIL[ctidx][i]) for i in range(len(ecent))]
				bContRad2=[getPercContRadius(perc2,ecent[i],bc0,bc1,beta,bNTAIL[ctidx][i],bSCORE[ctidx][i],bSTAIL[ctidx][i],bGCORE[ctidx][i],bGTAIL[ctidx][i]) for i in range(len(ecent))]
				ContRad2=[(B+F)/2. for F,B in zip(fContRad2,bContRad2)]
			#check the requested number of points
			npts=(len(el) if npts==None else npts)
			if npts>len(el):
				print "Requested number of points (%i) is greater than number of energy bins in IRFs (%i), will produce plot with %i points."%(npts,len(el),len(el))
				npts=len(el)
			step=int(len(el)/npts)
			if step>1:
				avecont=[sum([ContRad[i+step*j] for i in range(step)])/float(step) for j in range(len(el)/step)]
				bins=array([el[i*step] for i in range(len(el)/step)]+[eh[-1]])
				if perc!=None:
					avecont2=[sum([ContRad2[i+step*j] for i in range(step)])/float(step) for j in range(len(el)/step)]
					avecont=[c/C for c,C in zip(avecont2,avecont)]
			else:
				if perc2==None:
					avecont=ContRad
				else:
					avecont=[c/C for c,C in zip(ContRad2,ContRad)]
				bins=array([e for e in el]+[eh[-1]])
			
		else:
			#now do the same but this time for a particular energy and varying costheta
			#need to find the correct energy bin
			eidx=-1
			for i in range(len(el)):
				if value>=el[i] and value<eh[i]:
					eidx=i
					break
				if i==len(el)-1 and eidx==-1:
					eidx=len(el)-1
			print 'Using PSF values from bin with energy [%f,%f] MeV.'%(el[eidx],eh[eidx])
			#now get the requisite perecentage containment radii
			ecent=10**((log10(el[eidx])+log10(eh[eidx]))/2.)
			fContRad=[getPercContRadius(perc,ecent,fc0,fc1,beta,fNTAIL[i][eidx],fSCORE[i][eidx],fSTAIL[i][eidx],fGCORE[i][eidx],fGTAIL[i][eidx]) for i in range(len(ctl))]
			bContRad=[getPercContRadius(perc,ecent,bc0,bc1,beta,bNTAIL[i][eidx],bSCORE[i][eidx],bSTAIL[i][eidx],bGCORE[i][eidx],bGTAIL[i][eidx]) for i in range(len(ctl))]
			ContRad=[(B+F)/2. for F,B in zip(fContRad,bContRad)]
			if perc2!=None:
				fContRad2=[getPercContRadius(perc2,ecent,fc0,fc1,beta,fNTAIL[i][eidx],fSCORE[i][eidx],fSTAIL[i][eidx],fGCORE[i][eidx],fGTAIL[i][eidx]) for i in range(len(ctl))]
				bContRad2=[getPercContRadius(perc2,ecent,bc0,bc1,beta,bNTAIL[i][eidx],bSCORE[i][eidx],bSTAIL[i][eidx],bGCORE[i][eidx],bGTAIL[i][eidx]) for i in range(len(ctl))]
				ContRad2=[(B+F)/2. for F,B in zip(fContRad2,bContRad2)]
			#check the requested number of points
			npts=(len(ctl) if npts==None else npts)
			if npts>len(ctl):
				print "Requested number of points (%i) is greater than number of costheta bins in IRFs (%i), will produce plot with %i points."%(npts,len(ctl),len(ctl))
				npts=len(ctl)
			step=int(len(ctl)/npts)
			if step>1:
				avecont=[sum([ContRad[i+step*j] for i in range(step)])/float(step) for j in range(len(ctl)/step)]
				bins=array([ctl[i*step] for i in range(len(ctl)/step)]+[cth[-1]])
				bins=array([r2d*arccos(bins[len(bins)-1-i]) for i in range(len(bins))])
				if perc2!=None:
					avecont2=[sum([ContRad2[i+step*j] for i in range(step)])/float(step) for j in range(len(ctl)/step)]
					avecont=[c/C for c,C in zip(avecont2,avecont)]
				avecont=[avecont[len(avecont)-1-i] for i in range(len(avecont))]
			else:
				if perc2==None:
					avecont=ContRad
				else:
					avecont=[c/C for c,C in zip(ContRad2,ContRad)]
				avecont=[avecont[len(avecont)-1-i] for i in range(len(avecont))]
				bins=array([c for c in ctl]+[cth[-1]])
				bins=array([r2d*arccos(bins[len(bins)-1-i]) for i in range(len(bins))])
	else:
		conv=('front' if front else 'back')
		f_psf=pyfits.open(CALDB+'/data/glast/lat/bcf/psf/psf_%s_%s.fits'%(IRFs,conv))
		#get the energy and costheta bin edges
		el=f_psf[1].data.field('energ_lo')[0]
		eh=f_psf[1].data.field('energ_hi')[0]
		ctl=f_psf[1].data.field('ctheta_lo')[0]
		cth=f_psf[1].data.field('ctheta_hi')[0]
		NTAIL=f_psf[1].data.field('NTAIL')[0]
		SCORE=f_psf[1].data.field('SCORE')[0]
		STAIL=f_psf[1].data.field('STAIL')[0]
		GCORE=f_psf[1].data.field('GCORE')[0]
		GTAIL=f_psf[1].data.field('GTAIL')[0]
		psfpars=f_psf[2].data.field('PSFSCALE')
		if front:
			c0=psfpars[0][0]
			c1=psfpars[0][1]
		else:
			c0=psfpars[0][2]
			c1=psfpars[0][3]
		beta=psfpars[0][4]
		if ptype=='energy':
			#need to find correct costheta bin
			ctidx=-1
			for j in range(len(ctl)-1):
				if value>=ctl[j] and value<cth[j]:
					ctidx=j
					break
				if j==len(ctl) and ctidx==-1:
					ctidx=len(ctl)-1
			print 'Using PSF values from bin with costheta [%f,%f].'%(ctl[ctidx],cth[ctidx])
			#now get the requisite perecentage containment radii
			ecent=[10**((log10(E)+log10(e))/2.) for e,E in zip(el,eh)]
			ContRad=[getPercContRadius(perc,ecent[i],c0,c1,beta,NTAIL[ctidx][i],SCORE[ctidx][i],STAIL[ctidx][i],GCORE[ctidx][i],GTAIL[ctidx][i]) for i in range(len(ecent))]
			if perc2!=None:
				ContRad2=[getPercContRadius(perc2,ecent[i],c0,c1,beta,NTAIL[ctidx][i],SCORE[ctidx][i],STAIL[ctidx][i],GCORE[ctidx][i],GTAIL[ctidx][i]) for i in range(len(ecent))]
			#check the requested number of points
			npts=(len(el) if npts==None else npts)
			if npts>len(el):
				print "Requested number of points (%i) is greater than number of energy bins in IRFs (%i), will produce plot with %i points."%(npts,len(el),len(el))
				npts=len(el)
			step=int(len(el)/npts)
			if step>1:
				avecont=[sum([ContRad[i+step*j] for i in range(step)])/float(step) for j in range(len(el)/step)]
				bins=array([el[i*step] for i in range(len(el)/step)]+[eh[-1]])
				if perc2!=None:
					avecont2=[sum([ContRad2[i+step*j] for i in range(step)])/float(step) for j in range(len(el)/step)]
					avecont=[c/C for c,C in zip(avecont2,avecont)]
			else:
				if perc2==None:
					avecont=ContRad
				else:
					avecont=[c/C for c,C in zip(ContRad2,ContRad)]
				bins=array([e for e in el]+[eh[-1]])
			
		else:
			#now do the same but this time for a particular energy and varying costheta
			#need to find the correct energy bin
			eidx=-1
			for i in range(len(el)):
				if value>=el[i] and value<eh[i]:
					eidx=i
					break
				if i==len(el)-1 and eidx==-1:
					eidx=len(el)-1
			print 'Using PSF values from bin with energy [%f,%f] MeV.'%(el[eidx],eh[eidx])
			#now get the requisite perecentage containment radii
			ecent=10**((log10(el[eidx])+log10(eh[eidx]))/2.)
			ContRad=[getPercContRadius(perc,ecent,c0,c1,beta,NTAIL[i][eidx],SCORE[i][eidx],STAIL[i][eidx],GCORE[i][eidx],GTAIL[i][eidx]) for i in range(len(ctl))]
			if perc2!=None:
				ContRad2=[getPercContRadius(perc2,ecent,c0,c1,beta,NTAIL[i][eidx],SCORE[i][eidx],STAIL[i][eidx],GCORE[i][eidx],GTAIL[i][eidx]) for i in range(len(ctl))]
			#check the requested number of points
			npts=(len(ctl) if npts==None else npts)
			if npts>len(ctl):
				print "Requested number of points (%i) is greater than number of costheta bins in IRFs (%i), will produce plot with %i points."%(npts,len(ctl),len(ctl))
				npts=len(ctl)
			step=int(len(ctl)/npts)
			if step>1:
				avecont=[sum([ContRad[i+step*j] for i in range(step)])/float(step) for j in range(len(ctl)/step)]
				bins=array([ctl[i*step] for i in range(len(ctl)/step)]+[cth[-1]])
				bins=array([r2d*arccos(bins[len(bins)-1-i]) for i in range(len(bins))])
				if perc2!=None:
					avecont2=[sum([ContRad2[i+step*j] for i in range(step)])/float(step) for j in range(len(ctl)/step)]
					avecont=[c/C for c,C in zip(avecont2,avecont)]
				avecont=[avecont[len(avecont)-1-i] for i in range(len(avecont))]
			else:
				if perc2==None:
					avecont=ContRad
				else:
					avecont=[c/C for c,C in zip(ContRad2,ContRad)]
				avecont=[avecont[len(avecont)-1-i] for i in range(len(avecont))]
				bins=array([c for c in ctl]+[cth[-1]])
				bins=array([r2d*arccos(bins[len(bins)-1-i]) for i in range(len(bins))])
	#now make histogram
	hist=TH1F('hist','',len(avecont),bins)
	for i,c in enumerate(avecont):
		hist.SetBinContent(i+1,c)
	#set axis titles and do other style stuff
	hist.GetXaxis().CenterTitle()
	hist.GetYaxis().CenterTitle()
	xtitle=('Energy (MeV)' if ptype=='energy' else 'cos(#theta)')
	hist.SetXTitle(xtitle)
	ytitle=('%.1f'%(perc*100)+'%' +'Containment Angle (#circ)' if perc2==None else '%.1f'%(perc2*100)+'%'+'/%.1f'%(perc*100)+'%'+'Containment Angle Ratio')
	hist.SetYTitle(ytitle)
	hist.SetLineColor(1)
	hist.SetMarkerColor(1)
	hist.SetMarkerStyle(20)
	if plot:
		if perc2==None:
			can=TCanvas('can','%s %s %.1f'%(IRFs,conv,perc*100) +'% Containment Angle'+' vs. %s'%ptype,1000,800)
			can.SetLogy(1)
		else:
			can=TCanvas('can','%s %s %.1f'%(IRFs,conv,perc2*100) +'%'+'/%.1f'%(perc*100)+'% Containment Angle Ratio'+' vs. %s'%ptype,1000,800)
		can.SetTicks(1,1)
		if ptype=='energy':
			can.SetLogx(1)
		hist.Draw('pl')
		SetOwnership(can,False)
		SetOwnership(hist,False)
	return hist


#function to get the angular separation in degrees for a specified containment percentage with given PSF parameters
#note however that the percentage supplied is in fractional form
#arguments are:
#perc (float) -- desired containment percentage expressed as a fraction (e.g., 0.68 not 68 for 68% containment)
#energy (float) -- energy in MeV at which the requested contaiment radius is being asked for
#c0 (float) -- c0 parameter for scaling function
#c1 (float) -- c1 parameter for scaling function
#beta (float) -- beta parameter for scaling function
#ntail (float) -- normalization of tail King function from PSF probability function fit
#score (float) -- sigma value of core King function from PSF probability function fit
#stail (float) -- sigma value of tail King function from PSF probability function fit
#gcore (float) -- gamma value of core King function from PSF probability function fit
#gtail (float) -- gamma value of tail King function from PSF probability function fit
#step (float, optional) -- step size in scaled angular deviation, default is 0.001
def getPercContRadius(perc,energy,c0,c1,beta,ntail,score,stail,gcore,gtail,step=1e-3):
	#make an index to iterate on and a success flag
	i=0
	done=False
	while not done:
		#start at zero and slowly increase
		x=i*step
		#calculate the cdf of the global psf distribution with these particular parameters
		contain=psfcdf(x,ntail,score,stail,gcore,gtail)
		#if we're above the desired percentage, call it good, else increase i and try again
		#since we're starting from zero with small steps this should be fine, might need to fine tune step size though
		done=contain>perc#for now just do this, consider more careful/clever stop parameter/way to find best x later
		i+=1
	#return the containment back as an angular deviation in degrees
	#from solving for deltanu in Eq. 35 of section 6.1.1 in the pass7 performance paper
	#note that the scaling function is the only reason we need to energy
	#could have returned x and let the user rescale, but likely safer this way
	return x*S_P(energy,c0,c1,beta)*r2d

#make a class for the PSF probability function
#the _pdf is the core and tail King functions from the pass7 performance paper, Eq. 38 in section 6.1.1
#the _cdf is the integral of that from 0 to x, don't forget the factor of 2*pi*x from the integral over solid angle
#arg0=ntail
#arg1=score
#arg2=stail
#arg3=gcore
#arg4=gtail
#class PSF(rvc):#this should work, but with FSSC STs something wasn't right with scipy.stats
	#def _pdf(self,x,arg0,arg1,arg2,arg3,arg4):
		#fc=fcore(arg0,arg2,arg1)
		#return fc*King(x,arg1,arg3)+(1.-fc)*King(x,arg2,arg4)
	#def _cdf(self,x,arg0,arg1,arg2,arg3,arg4):
		#fc=fcore(arg0,arg2,arg1)
		#return fc*(1.-(1.+x**2/(2.*arg3*arg1**2))**(1.-arg3))+(1.-fc)*(1.-(1.+x**2/(2.*arg4*arg2**2))**(1.-arg4))

def psfpdf(x,ntail,score,stail,gcore,gtail):
	fc=fcore(ntail,stail,score)
	return fc*King(x,score,gcore)+(1.-fc)*King(x,stail,gtail)

def psfcdf(x,ntail,score,stail,gcore,gtail):
	fc=fcore(ntail,stail,score)
	return fc*(1.-(1.+x**2/(2.*gcore*score**2))**(1.-gcore))+(1.-fc)*(1.-(1.+x**2/(2.*gtail*stail**2))**(1.-gtail))


#make an instance of the PSF distribution, define it globally to make things easier as to use it you have to specify
#the King function parameters, don't give it one automatically or anything
#psf=PSF(name='psf',a=0.,xa=0)

#King function as given in pass7 performance paper, Eq. 36 in section 6.1.1
def King(x,sigma,gamma):
	return (1./(2*pi*sigma**2))*(1.-1./gamma)*((1.+(1./(2.*gamma))*(x**2/sigma**2))**-gamma)
	#return ((2*pi*sigma**2)**-1)*(1.-gamma**-1)*((1.+((2.*gamma)**-1)*(x**2/sigma**2))**-gamma)

#LAT PSF scaling function as given in pass7 performance paper, Eq. 34 in section 6.1.1
def S_P(energy,c0,c1,beta):
	beta=(-beta if beta<0 else beta)
	return ((c0*(energy/100.)**-beta)**2+c1**2)**0.5

#function to get fcore value for joint King function for LAT PSF from pass7 performance paper, Eq. 39 section 6.1.1
def fcore(NTAIL,STAIL,SCORE):
	return 1./(1.+NTAIL*(STAIL**2/SCORE**2))
	#return (1.+NTAIL*(STAIL**2/SCORE**2))**-1



####################################################################################
####################################################################################
#############end functions for making psf plots#####################################
####################################################################################
####################################################################################


####################################################################################
####################################################################################
#############start functions for making edisp plots#################################
####################################################################################
####################################################################################

#ummm...I think I'm actually not going to do this part...too hard...unless I can use the pyIrfLoader

#energy dispersion scaling function from the pass7 performance paper, Eq. 48 in section 7.1.1
#easier to have it defined one to avoid errors typing the damn thing out in one of several places
def S_D(logene,costheta,c0,c1,c2,c3,c4,c5):
	return c0*logene**2+c1*costheta**2+c2*logene+c3*costheta+c4*logene*costheta+c5



####################################################################################
####################################################################################
#############end functions for making edisp plots###################################
####################################################################################
####################################################################################


#################################
###general use style functions###
#################################

def BPalette():
	r=array([0.,0.0,1.0,1.0,1.0])
	b=array([0., 1.0, 0.0, 0.0, 1.0])
	g=array([0., 0.0, 0.0, 1.0, 1.0])
	stop=array([0.,.25,.50,.75,1.0])
	TColor.CreateGradientColorTable(5,stop,r,g,b,100)
	return
	
def GrayPalette():
	R=array([0.,1.])
	G=array([0.,1.])
	B=array([0.,1.])
	Stop=array([0.,1.])
	TColor.CreateGradientColorTable(2,Stop,R,G,B,100)
	return