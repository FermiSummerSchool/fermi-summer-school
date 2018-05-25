from ROOT import TH1F,TH2F,TColor,gStyle,SetOwnership,TCanvas,TGraph
from numpy import array,pi,log10,zeros,arange,arccos,cos
import pyIrfLoader

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

#function to create npts points from xmin to xmax (inclusive) with equal spacing in log10 space
def log10_array(npts,xmin,xmax):
	xstep = log10(xmax/xmin)/(npts - 1)
	return xmin*10**(arange(npts,dtype=float)*xstep)


####################################################################################
####################################################################################
#############start functions for making effective area plots########################
####################################################################################
####################################################################################

#function to make a 2D plot of effective area with y axis = cos(theta) and x axis = energy (in MeV)
#can do for FRONT, BACK, or FRONT+BACK for given IRFs
#arguments are:
#IRFs (str) -- IRFs you want to plot (e.g., P7REP_SOURCE_V15)
#nebins (int, optional) -- number of energy bins for the plot
#nctbins (int, optional) -- number of cos(theta) bins for the plot
#emin (float, optional) -- minimum energy in MeV
#emax (float, optional) -- maximum energy in MeV
#ctmin (float, optional) -- minimum cos(theta)
#ctmax (float, optional) -- maximum cos(theta)
#front (bool, optional) -- flag to get effective area for front section, default is true
#back (bool, optional) -- flag to get effective area for back section, default is false
#note, if both back and front are true, plot will be for back+front effective area
#plot (bool, optional) -- flag to draw the 2D histogram
#color (bool, optional) -- flag to plot in color or grayscale, only important if plot==True
def AeffPlot2D(IRFs,nebins=64,nctbins=32,emin=1.e2,emax=1.e5,ctmin=0.2,ctmax=1.0,front=True,back=False,plot=True,color=True):
	#should probably put some reasonability test on min and max values...
	ebins=log10_array(nebins+1,emin,emax)
	#get the logarithmic midpoints for energy bins
	ecents=[10**((log10(e)+log10(E))/2.) for e,E in zip(ebins[:-1],ebins[1:])]
	ctbins=[ctmin+i*(ctmax-ctmin)/nctbins for i in range(nctbins+1)]
	#get theta values in degrees of centers of cos(theta) bins
	thetas=[r2d*arccos((c+C)/2.) for c,C in zip(ctbins[:-1],ctbins[1:])]
	#get the available IRFs
	pyIrfLoader.Loader_go()
	if front and back:
		conv='FRONT + BACK'
		#get front and back IRFs
		irfs_front=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::FRONT")
		irfs_back=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::BACK")
		#get effective area for both sections and then add
		Aeff_front=getAeff_nophi(irfs_front.aeff(),ecents,thetas)
		Aeff_back=getAeff_nophi(irfs_back.aeff(),ecents,thetas)
		Aeff=Aeff_front+Aeff_back
	else:
		#get irfs for proper section and get effective area
		conv=('FRONT' if front else 'BACK')
		irfs=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::"+conv)
		Aeff=getAeff_nophi(irfs.aeff(),ecents,thetas)
	#make and fill hist with some style options
	aehist=TH2F('aehist','',nebins,array(ebins),nctbins,array(ctbins))
	for c in range(nctbins):
		for e in range(nebins):
			aehist.SetBinContent(e+1,c+1,Aeff[e][c])
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

#function to make a 1D plot of effective area versus either energy or theta
#can do for FRONT, BACK, or FRONT+BACK for given IRFs
#arguments are:
#IRFs (str) -- IRFs you want to plot (e.g., P7REP_SOURCE_V15)
#nbins (int) -- number of bins
#ptype (str) -- which variable to plot against, only "energy" or "theta" accepted
#pmin (float) -- minimum value for x axis, note this is the low edge of the first bin
#pmax (float) -- maximum value for x axis, note this is the high edge of the first bin
#value (float) -- value of other variable you're not plotting against (e.g., if ptype=='energy' then this is the theta value)
#front (bool, optional) -- flag to get effective area for front section, default is true
#back (bool, optional) -- flag to get effective area for back section, default is false
#note, if both back and front are true, plot will be for back+front effective area
#plot (bool, optional) -- flag to draw the 2D histogram
def AeffPlot1D(IRFs,nbins,ptype,pmin,pmax,value,front=True,back=False,plot=True):
	#check pytpe input
	if ptype not in ['energy','theta']:
		print 'Error, ptype input must be either "energy" or "theta".'
		return None
	#bin accordingly depending on choice
	if ptype=='energy':
		bins=log10_array(nbins+1,pmin,pmax)
		cents=[10**((log10(e)+log10(E))/2.) for e,E in zip(bins[:-1],bins[1:])]
	else:
		bins=[float(pmin+i*(pmax-pmin)/float(nbins)) for i in range(nbins+1)]
		cents=[(c+C)/2. for c,C in zip(bins[:-1],bins[1:])]
	#access the available IRFS
	pyIrfLoader.Loader_go()
	if front and back:
		#get FRONT and BACK IRFs
		conv='FRONT + BACK'
		irfs_front=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::FRONT")
		irfs_back=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::BACK")
		#get effective area in bins, either at one theta or at one energy
		if ptype=='energy':
			Aeff_front=getAeff_nophi(irfs_front.aeff(),cents,[value])
			Aeff_back=getAeff_nophi(irfs_back.aeff(),cents,[value])
		else:
			Aeff_front=getAeff_nophi(irfs_front.aeff(),[value],cents)
			Aeff_back=getAeff_nophi(irfs_back.aeff(),[value],cents)
		#add effective areas
		Aeff=Aeff_front+Aeff_back
	else:
		#get appropriate IRFs
		conv=('FRONT' if front else "BACK")
		irfs=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::"+conv)
		#get effective area for either one theta or one energy
		if ptype=='energy':
			Aeff=getAeff_nophi(irfs.aeff(),cents,[value])
		else:
			Aeff=getAeff_nophi(irfs.aeff(),[value],cents)
	#make histogram and fill
	aehist=TH1F('aehist','',nbins,array(bins))
	for i,a in enumerate(Aeff):
		aehist.SetBinContent(i+1,a)
	#set axis titles and do other style stuff
	aehist.GetXaxis().CenterTitle()
	aehist.GetYaxis().CenterTitle()
	xtitle=('Energy (MeV)' if ptype=='energy' else '#theta (#circ)')
	aehist.SetXTitle(xtitle)
	aehist.SetYTitle('A_{eff} (m^{2} )')
	aehist.SetLineColor(1)
	aehist.SetMarkerColor(1)
	aehist.SetMarkerStyle(20)
	if plot:
		can=TCanvas('can','%s %s Effective Area vs. %s'%(IRFs,conv,ptype),1000,800)
		can.SetTicks(1,1)
		if ptype=='energy':
			can.SetLogx(1)
		aehist.Draw('pl')
		SetOwnership(can,False)
		SetOwnership(aehist,False)
	return aehist

#function to make a 1D plot of effective area versus phi for a given energy and theta
#can do for FRONT, BACK, or FRONT+BACK for given IRFs
#arguments are:
#IRFs (str) -- IRFs you want to plot (e.g., P7REP_SOURCE_V15)
#nph (int) -- number of phi points from 0 deg to 90 deg (inclusive)
#energy (float) -- energy for effective area calculation in MeV
#theta (float) -- theta for effective area calculation in degrees
#front (bool, optional) -- flag to get effective area for front section, default is true
#back (bool, optional) -- flag to get effective area for back section, default is false
#note, if both back and front are true, plot will be for back+front effective area
#plot (bool, optional) -- flag to draw the 2D histogram
def AeffPhiDepPlot(IRFs,nph,energy,theta,front=True,back=False,plot=True):
	#get phi points
	phis=[0+i*(90./(nph-1.)) for i in range(nph)]
	#access available IRFs
	pyIrfLoader.Loader_go()
	if front and back:
		#get front and back IRFs, get effective area for phis for given energy and theta
		conv='FRONT + BACK'
		irfs_front=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::FRONT")
		irfs_back=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::BACK")
		Aeff_front=getAeff_phi(irfs_front.aeff(),energy,theta,phis)
		Aeff_back=getAeff_phi(irfs_back.aeff(),energy,theta,phis)
		Aeff=Aeff_front+Aeff_back
	else:
		#access the necessary IRFs, get effective area for phis for given energy and theta
		conv=('FRONT' if front else "BACK")
		irfs=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::"+conv)
		Aeff=getAeff_phi(irfs.aeff(),energy,theta,phis)
	#expand the phi values out to 360 degrees, copy effective area values too
	phibins=array(phis+[p+90. for p in phis[1:]]+[p+180. for p in phis[1:]]+[p+270. for p in phis[1:]])
	AeffValues=array([a for a  in Aeff]+[a for a  in Aeff[1:]]+[a for a  in Aeff[1:]]+[a for a  in Aeff[1:]])
	#make and fill phi graph
	phigraph=TGraph(len(phibins),phibins,AeffValues)
	phigraph.SetMarkerStyle(20)
	phigraph.SetMarkerColor(1)
	phigraph.SetLineColor(1)
	if plot:
		#if plot set to True, draw and make some style choices
		phcan=TCanvas('phcan','%s %s Effective Area vs. phi for E = %.1f MeV and theta = %.1f deg'%(IRFs,conv,energy,theta),1000,800)
		phcan.SetTicks(1,1)
		pdummy=TH1F('pdummy','',1000,0,360)
		pdummy.SetXTitle('#phi (#circ)')
		pdummy.GetXaxis().CenterTitle()
		pdummy.SetYTitle('A_{eff} (m^{2} )')
		pdummy.GetYaxis().CenterTitle()
		pdummy.GetYaxis().SetRangeUser(0.,1.)
		pdummy.Draw()
		phigraph.Draw('plsame')
		SetOwnership(phigraph,False)
		SetOwnership(phcan,False)
		SetOwnership(pdummy,False)
	return phigraph

#function to get effective area with no phi dependence for given energy and theta value(s)
#can produce either a 2D or 1D array
#arguments are:
#ae (IAeff instance) -- this assumes you have already made an IRFs instance, need to pass in the aeff() function
#ene (float array) -- array of energy values in MeV
#th (float array) -- array of theta values in degrees
def getAeff_nophi(ae,ene,th):
	#make sure we're not using phi dependence
	ae.setPhiDependence(0)
	#if both ene and th have more than one entry, we're making a 2D array
	if len(ene)>1 and len(th)>1:
		#make 2D array of zeros with appropriate dimensions and fill
		values=zeros([len(ene),len(th)])
		for t in range(len(th)):
			for e in range(len(ene)):
				values[e][t]=ae.value(ene[e],th[t],0.)/10000.#divide by 10000. to get units of m^2
	#failed last if statment, so one of these guys only has 1 entry
	#if ene has more than one energy, then we're making a 1D array for a particular theta but a range of energies
	elif len(ene)>1:
		values=array([ae.value(e,th[0],0.)/10000. for e in ene])
	#otherwise we're making a 1D array at a particular energy for a range of thetas
	else:
		values=array([ae.value(ene[0],t,0.)/10000. for t in th])
	return values

#function to return an array of effective area values for a set of phi values at a particular energy and theta
#could end up being an array with one entry...just in case
#arguments are:
#ae (IAeff instance) -- this assumes you have already made an IRFs instance, need to pass in the aeff() function
#ene (float) -- energy in MeV
#th (float) -- theta in degrees
#ph (float array) -- array of phi values
def getAeff_phi(ae,ene,th,phi):
	#make sure phi dependence is turned on
	ae.setPhiDependence(1)
	#if we have more than one phi value, make a 1D array of effective area values
	if len(phi)>1:
		values=array([ae.value(ene,th,p)/10000. for p in phi])#divide by 10000. to get units of m^2
		return values
	#otherwise we're just returning one value...not sure why anyone would call this function just for that...
	else:
		return ae.value(ene,th,phi[0])/10000.

####################################################################################
####################################################################################
#############end functions for making effective area plots##########################
####################################################################################
####################################################################################

####################################################################################
####################################################################################
#############start functions for making point-spread function plots#################
####################################################################################
####################################################################################

#function to produce a 1D histogram of containment angle for specified IRFs and specified containment percentage
#can do FRONT, BACK, or FRONT & BACK.  Can also make it a plot of containment angle ratio for specified containment percentages
#arguments are:
#IRFs (str) -- IRFs you want to plot (e.g., P7REP_SOURCE_V15)
#perc (float) -- fractional containment
#ptype (str) -- which variable to plot against, only "energy" or "theta" accepted
#pmin (float) -- minimum value for x axis, note this is the low edge of the first bin
#pmax (float) -- maximum value for x axis, note this is the high edge of the first bin
#value (float) -- value of other variable you're not plotting against (e.g., if ptype=='energy' then this is the theta value)
#nbins (int) -- number of bins/points to plot
#perc2 (float) -- if specified, this is the second fractional containment for a ratio plot
#value plotted will be perc2 containment angle/perc containment angle
#front (bool, optional) -- flag to get effective area for front section, default is true
#back (bool, optional) -- flag to get effective area for back section, default is false
#note, if both back and front are true, plot will be for back+front effective area
#plot (bool, optional) -- flag to draw the 2D histogram
def PSFPlot(IRFs,perc,ptype,pmin,pmax,value,nbins,perc2=None,front=True,back=False,plot=True):
	#check choice of ptype
	if ptype not in ['energy','theta']:
		print 'Error, ptype input must be either "energy" or "theta".'
		return None
	#make bins accordingly
	if ptype=='energy':
		bins=log10_array(nbins+1,pmin,pmax)
		cents=[10**((log10(e)+log10(E))/2.) for e,E in zip(bins[:-1],bins[1:])]
	else:
		bins=[pmin+i*(pmax-pmin)/nctbins for i in range(nbins+1)]
		cents=[(c+C)/2. for c,C in zip(bins[:-1],bins[1:])]
	#access the available IRFs
	pyIrfLoader.Loader_go()
	if front and back:
		print "Calculating containment ratios for FRONT+BACK as the average"
		print "this is approximately correct near 685","containment but underestimates for larger containment percentages"
		#get front and back IRFs
		conv='FRONT + BACK'
		irfs_front=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::FRONT")
		irfs_back=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::BACK")
		#get the front and back containment angles, either for one theta and range of energies or vice versa
		#then average the 2, makes the assumption that we'll have about the same number of events in both sections
		if ptype=='energy':
			cont_front=[psfContAng(irfs_front.psf(),perc,e,value) for e in cents]
			cont_back=[psfContAng(irfs_back.psf(),perc,e,value) for e in cents]
			cont=array([(f+b)/2. for f,b in zip(cont_front,cont_back)])
			#if we have a value for perc2, calculate those containment angles
			#average these values for front and back before calculating the ratio
			if perc2!=None:
				cont2_front=array([psfContAng(irfs_front.psf(),perc2,e,value) for e in cents])
				cont2_back=array([psfContAng(irfs_back.psf(),perc2,e,value) for e in cents])
				cont2=[(f+b)/2. for f,b in zip(cont2_front,cont2_back)]
				cont=array([c/C for c,C in zip(cont2,cont)])
		else:
			cont_front=[psfContAng(irfs_front.psf(),perc,value,th) for th in cents]
			cont_back=[psfContAng(irfs_back.psf(),perc,value,th) for th in cents]
			cont=array([(f+b)/2. for f,b in zip(cont_front,cont_back)])
			if perc2!=None:
				cont2_front=[psfContAng(irfs_front.psf(),perc2,value,th) for th in cents]
				cont2_back=[psfContAng(irfs_back.psf(),perc2,value,th) for th in cents]
				cont2=[(f+b)/2. for f,b in zip(cont2_front,cont2_back)]
				cont=array([c/C for c,C in zip(cont2,cont)])
	else:
		#get appropriate IRFs
		conv=('FRONT' if front else "BACK")
		irfs=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::"+conv)
		#get containtment angles
		if ptype=='energy':
			cont=array([psfContAng(irfs.psf(),perc,e,value) for e in cents])
			if perc2!=None:
				#if we actually have a value for perc2, get second set of containment angles then calculate ratio
				cont2=[psfContAng(irfs.psf(),perc2,e,value) for e in cents]
				cont=array([c/C for c,C in zip(cont2,cont)])
		else:
			cont=array([psfContAng(irfs.psf(),perc,value,th) for th in cents])
			if perc2!=None:
				cont2=[psfContAng(irfs.psf(),perc2,value,th) for th in cents]
				cont=array([c/C for c,C in zip(cont2,cont)])
	#make the histogram and fill
	psfhist=TH1F('psfhist','',nbins,bins)
	for i,c in enumerate(cont):
		psfhist.SetBinContent(i+1,c)
	#set axis titles and do other style stuff
	psfhist.GetXaxis().CenterTitle()
	psfhist.GetYaxis().CenterTitle()
	xtitle=('Energy (MeV)' if ptype=='energy' else '#theta (#circ)')
	psfhist.SetXTitle(xtitle)
	ytitle=('%.1f'%(perc*100)+'% Containment Angle (#circ)' if perc2==None else '%.1f'%(perc2*100)+'%'+'/%.1f'%(perc*100)+'%'+'Containment Angle Ratio')
	psfhist.SetYTitle(ytitle)
	psfhist.SetLineColor(1)
	psfhist.SetMarkerColor(1)
	psfhist.SetMarkerStyle(20)
	if plot:
		if perc2==None:
			can=TCanvas('can','%s %s %.1f'%(IRFs,conv,perc*100) +'% Containment Angle'+' vs. %s'%ptype,1000,800)
			can.SetLogy(1)
		else:
			can=TCanvas('can','%s %s %.1f'%(IRFs,conv,perc2*100) +'%'+'/%.1f'%(perc*100)+'% Containment Angle Ratio'+' vs. %s'%ptype,1000,800)
		can.SetTicks(1,1)
		if ptype=='energy':
			can.SetLogx(1)
		psfhist.Draw('pl')
		SetOwnership(can,False)
		SetOwnership(psfhist,False)
	return psfhist

#this function adapted from Luca Baldini's IRFPlotter.py tool
#I say adapted but it didn't really change much
#arguments are:
#psf (IPSF instance) -- this assumes you have already made an IRFs instance, need to pass in the psf() function
#perc (float) -- fractional containment you want to find the angle for
#ene (float) -- energy in MeV
#th (float) -- theta in degrees
#threshold (float, optional) -- threshold for stopping search
def psfContAng(psf,perc,ene,th,threshold=1e-6):
        #comments below taken directly from Luca's code, good info to propagate
        """ Return the PSF containment radius for a given quantile at a given
        energy, theta and phi.
        Since we know the derivative of the integral (just the psf value)
        we use the Newton method to find the quantile:
        http://en.wikipedia.org/wiki/Newton%27s_method
        """
        r=1.0
        dr=0.5
        while abs(dr)>threshold:
            dr=(psf.angularIntegral(ene,th,0.,r)-perc)/psf.value(r,ene,th,0.)
            r-=dr
            dr/= r
            if r<0:
                r=0.1
                dr=0.5
        return r


####################################################################################
####################################################################################
#############end functions for making point-spread function plots###################
####################################################################################
####################################################################################



####################################################################################
####################################################################################
#############start functions for making energy dispersion plots#####################
####################################################################################
####################################################################################

#function to make and plot 1D histogram of energy resolution vs either energy or theta
#can do FRONT or BACK but currenty not FRONT & BACK.
#arguments are:
#IRFs (str) -- IRFs you want to plot (e.g., P7REP_SOURCE_V15)
#perc (float) -- fractional containment
#ptype (str) -- which variable to plot against, only "energy" or "theta" accepted
#pmin (float) -- minimum value for x axis, note this is the low edge of the first bin
#pmax (float) -- maximum value for x axis, note this is the high edge of the first bin
#value (float) -- value of other variable you're not plotting against (e.g., if ptype=='energy' then this is the theta value)
#nbins (int) -- number of bins/points to plot
#front (bool, optional) -- flag to get effective area for front section, default is true
#back (bool, optional) -- flag to get effective area for back section, default is false
#note, if both back and front are true nothing will happen, that isn't implemented for EDISP in this code yet
#plot (bool, optional) -- flag to draw the 2D histogram
def EdispPlot(IRFs,perc,ptype,pmin,pmax,value,nbins,front=True,back=False,plot=True):
	#check ptype input value
	if ptype not in ['energy','theta']:
		print 'Error, ptype input must be either "energy" or "costheta".'
		return None
	#make bins accordingly
	if ptype=='energy':
		bins=log10_array(nbins+1,pmin,pmax)
		cents=[10**((log10(e)+log10(E))/2.) for e,E in zip(bins[:-1],bins[1:])]
	else:
		bins=array([float(pmin+i*(pmax-pmin)/float(nbins)) for i in range(nbins+1)])
		cents=[(c+C)/2. for c,C in zip(bins[:-1],bins[1:])]
	#access the available IRFs
	pyIrfLoader.Loader_go()
	if front and back:
		#get front and back IRFs, get the energy resolution for both, then average them
		print "oops, can't do FRONT and BACK right now, sorry, try again"
		return None
		#conv='FRONT + BACK'
		#irfs_front=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::FRONT")
		#irfs_back=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::BACK")
		##get for either one theta and range of energies or vice versa
		#if ptype=='energy':
			#res_front=[energyRes(irfs_front.edisp(),perc,e,value) for e in cents]
			#res_back=[energyRes(irfs_back.edisp(),perc,e,value) for e in cents]
			#res=[(f+b)/2. for f,b in zip(res_front,res_back)]
		#else:
			#res_front=[energyRes(irfs_front.edisp(),perc,value,th) for th in cents]
			#res_back=[energyRes(irfs_back.edisp(),perc,value,th) for th in cents]
			#res=[(f+b)/2. for f,b in zip(res_front,res_back)]
	else:
		#get the appropriate IRFs, get energy resolution based on ptype
		conv=('FRONT' if front else "BACK")
		irfs=pyIrfLoader.IrfsFactory.instance().create(IRFs+"::"+conv)
		if ptype=='energy':
			res=[energyRes(irfs.edisp(),perc,e,value) for e in cents]
		else:
			res=[energyRes(irfs.edisp(),perc,value,th) for th in cents]
	#make histogram and fill
	edhist=TH1F('edhist','',nbins,bins)
	for i,e in enumerate(res):
		edhist.SetBinContent(i+1,e)
	#set axis titles and do other style stuff
	edhist.GetXaxis().CenterTitle()
	edhist.GetYaxis().CenterTitle()
	xtitle=('Energy (MeV)' if ptype=='energy' else '#theta (#circ)')
	ytitle=('#frac{#Delta E}{E} (%.1f'%(perc*100)+'%'+' containtment for #theta = %.1f^{#circ} )'%value if ptype=='energy' else '#frac{#Delta E}{E} (%.1f'%(perc*100)+'%'+ 'containment for E = %.1f MeV)'%value)
	edhist.SetYTitle(ytitle)
	edhist.SetLineColor(1)
	edhist.SetMarkerColor(1)
	edhist.SetMarkerStyle(20)
	if plot:
		can=TCanvas('can','%s %s %.1f'%(IRFs,conv,perc*100) +'% Energy Resolution'+' vs. %s'%ptype,1000,800)
		can.SetTicks(1,1)
		if ptype=='energy':
			can.SetLogx(1)
		edhist.Draw('pl')
		SetOwnership(can,False)
		SetOwnership(edhist,False)
	return edhist

#this function adapted from Luca Baldini's IRFPlotter.py tool
#but tweaked to do more than just find energy of one quantile
#gets peak of edisp distribution, then finds full width of perc containment
#then returns half of that width as the energy resolution
#arguments are:
#edisp (IEdisp instance) -- this assumes you have already made an IRFs instance, need to pass in the edisp() function
#perc (float) -- fractional containment you want to find the angle for
#trueEne (float) -- assumed true energy of even in MeV
#th (float) -- theta in degrees
#threshold (float, optional) -- threshold for stopping search
#maxSteps (int, optional) -- maximum number of steps in search
def energyRes(edisp,perc,trueEne,th,threshold=1e-8,maxSteps=50):
	#comments below taken directly from Luca's code, good info to propagate
	""" Return the energy value for which the energy dispersion
	integral value is equal to a given quantile.
	
	Since we know the derivative of the integral (just the edisp value)
	we use the Newton method to find the quantile:
	http://en.wikipedia.org/wiki/Newton%27s_method
	"""
	#peak=edisp.meanAppEnergy(trueEne,th,0.)
	#ppeak=edisp.integral(0.,peak,trueEne,th,0.)
	energies=[trueEne/2.+i*(trueEne*2.-trueEne/2.)/1000 for i in range(1001)]
	values=[edisp.value(e,trueEne,th,.0) for e in energies]
	maxval=max(values)
	for i in range(len(values)):
		if values[i]==maxval:
			myval=energies[i]
			break
	#print myval,edisp.meanAppEnergy(trueEne,th,0.)
	ppeak=edisp.integral(0.,myval,trueEne,th,0.)
	pmin=ppeak-perc/2.
	pmax=ppeak+perc/2.
	emax=trueEne
	demax=1
	steps=0
	while abs(demax)>threshold:            
		demax=(edisp.integral(0,emax,trueEne,th,0.)-pmax)/edisp.value(emax,trueEne,th,0.)
		emax-=demax
		demax/=trueEne
		steps+=1
		if emax<0:
			emax=trueEne/10.
			demax=1
		if steps>maxSteps:
			break
	emin=trueEne
	demin=1
	steps=0
	while abs(demin)>threshold:            
		demin=(edisp.integral(0,emin,trueEne,th,0.)-pmin)/edisp.value(emin,trueEne,th,0.)
		emin-=demin
		demin/=trueEne
		steps+=1
		if emin<0:
			emin=trueEne/10.
			demin=1
		if steps>maxSteps:
			break
	return (emax-emin)/2./trueEne



####################################################################################
####################################################################################
#############end functions for making energy dispersion plots#######################
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