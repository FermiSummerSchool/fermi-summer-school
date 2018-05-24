#!/usr/bin/env python

class srcList:
	#arguments are:
	#sources (string, filename of LAT source list fits file in catalog format)
	#ft1 (string, filename of event file for which the xml will be used, only used to extract ROI info)
	#out (string, name of output xml file, defaults to mymodel.xml)
	def __init__(self,sources,ft1,out='mymodel.xml'):
		if not fileCheck(sources): #check that file exists
			print "Error:  %s not found." %sources
			return
		if fileCheck(out):
			print 'Warning: %s already exists, file will be overwritten if you proceed with makeModel.' %out
		self.srcs=sources
		self.out=out
		self.roi=getPos(ft1)
	
	#define a quick print function to make sure everything looks irght
	def Print(self):
		print 'Source list file: ',self.srcs
		print 'Output file name: ',self.out
		print 'Selecting %s degrees around (ra,dec)=(%s,%s)' %(self.roi[2],self.ra,self.dec)
	
	#make the xml file
	#arguments are:
	#GDfile (str) -- optional, location and name of Galactic diffuse model to use
	#GDname (str) -- optional, name of Galactic diffuse component to use in xml model
	#ISOfile (str) -- optional, location and name of Isotropic diffuse template to use
	#ISOname (str) -- optional, name of Isotropic diffuse component to use in xml model
	#normsOnly (bool) -- optional, flag to only set normalizations parameters free
	#extDir (str) -- optional, directory with extended source templates
	#radLim (float) -- optional, radius in degrees from center of ROI beyond which source parameters are fixed
	#maxRad (float) -- optional, absolute maximum radius beyond which sources are fixed, this may be necessary when doing binned analysis and a variable source beyond radLim would be set free but this source is beyond the boundaries of the square region used for the binned likelihood
	#ExtraRad (float) -- optional, radius beyond ROI radius in event file out to which sources will be included with fixed parameters, defaul tof 10 is good for analyses starting around 100 MeV, but for higher energy fits this can be decreased
	#sigFree (float) -- optional, significance below which source parameters are fixed, even if within radLim
	#varFree (float) -- optional, variability index above which source parameters are free, if beyond radLim and/or below sigFree only the normalization parameter is set free, currently not implemented for building from xml catalog
	#psForce (bool) -- optional, flag to force extended sources to be point sources
	#makeRegion (bool) -- optional, flag to also generate ds9 region file
	#GIndexFree (bool) -- optional, the Galactic diffuse is given a power-law spectral shape but the by default the index is frozen, setting this flag to True allows that to be free for additional freedom in diffuse fit
	#ApplyEDisp (boo) -- optional, flag to apply energy dispersion to free sources (except diffuse backgrounds) default is False.
	def makeModel(self,GDfile="$(FERMI_DIR)/refdata/fermi/galdiffuse/gll_iem_v06.fits",GDname='gll_iem_v06',ISOfile="$(FERMI_DIR)/refdata/fermi/galdiffuse/iso_P8R2_SOURCE_V6_v06.txt",ISOname='iso_P8R2_SOURCE_V6_v06',normsOnly=False,extDir='',radLim=-1,maxRad=None,ExtraRad=10,sigFree=5,varFree=True,psForce=False,makeRegion=True,GIndexFree=False,ApplyEDisp=False,wd='',oldNames=False):
		self.radLim=(self.roi[2] if radLim<=0 else radLim)
		self.maxRad=(self.radLim if maxRad==None else maxRad)
		if self.maxRad<self.radLim:
			print "NOTE: maxRad (%.1f deg) is less than radLim (%.1f deg), meaning maxRad parameter is useless"%(self.maxRad,self.radLim)
		self.var=varFree
		self.psF=psForce
		self.nO=normsOnly
		extDir=(extDir if extDir!='' else '$(LATEXTDIR)/Templates/')#make sure the default for FSSC STs is correct
		self.extD=(extDir if extDir[-1]=='/' else extDir+'/')
		self.ER=ExtraRad
		self.sig=sigFree
		self.reg=makeRegion
		self.GIF=GIndexFree
		self.ed=('true' if ApplyEDisp else 'false')
		if makeRegion:
			rhold=self.out.split('.')[:-1]
			rhold=rhold[0].split('/')[-1]#just in case the user has specified the full path for the output XML model file
			wd=(os.getcwd() if wd=='' else wd)
			self.regFile=wd+'/ROI_'
			for r in rhold:
				self.regFile+=r
			self.regFile+='.reg'
		print 'Creating file and adding sources from 3FGL'
		#want ability to use either the FITS or xml versions of the catalog
		#need to tweak FITS version to have new functionality and then work out xml version
		if self.srcs.split('.')[-1]=='xml':
			addSrcsXML(self,GDfile,GDname,ISOfile,ISOname,oldNames)
		else:
			addSrcsFITS(self,GDfile,GDname,ISOfile,ISOname,oldNames)
	
import pyfits
import os
from xml.dom import minidom
from xml.dom.minidom import parseString as pS
from numpy import floor,log10,cos,sin,arccos,pi,array,log,exp
acos=arccos
#import ROOT #note that this is only done to turn tab completion on for functions and filenames
print "This is make3FGLxml version 01r0."
print "The default diffuse model files and names are for pass 8 and assume you have v10r00p05 of the Fermi Science Tools or higher."
#print "For use with the gll_psc_v02.fit and gll_psc_v05.fit and later LAT catalog files."
#print "NOTE: You must have run gtselect on the event file you use as input."
d2r=pi/180.

def addSrcsXML(sL,GD,GDn,ISO,ISOn,oldNames=False):
	inputXml=minidom.parse(sL.srcs)
	outputXml=minidom.getDOMImplementation().createDocument(None,'source_library',None)
	outputXml.documentElement.setAttribute('title','source library')
	catalog=inputXml.getElementsByTagName('source')
	Sources={}
	ptSrcNum=0
	extSrcNum=0
	ed=sL.ed
	#normNames=['Prefactor','Integral','norm']
	#freePars=['Index','Index1','Cutoff','alpha','beta']
	for src in catalog:
		if src.getAttribute('type')=='PointSource':
			for p in src.getElementsByTagName('spatialModel')[0].getElementsByTagName('parameter'):
				if p.getAttribute('name')=='RA':
					srcRA=float(p.getAttribute('value'))
				if p.getAttribute('name')=='DEC':
					srcDEC=float(p.getAttribute('value'))
		else:
			srcDEC=float(src.getAttribute('DEC'))
			srcRA=float(src.getAttribute('RA'))
		dist=angsep(sL.roi[0],sL.roi[1],srcRA,srcDEC) #check that source is within ROI radius + 10 degress of ROI center
		if srcRA==sL.roi[0] and srcDEC==sL.roi[1]:
			dist=0.0
		if dist<=sL.roi[2]+sL.ER:
			spec=src.getElementsByTagName('spectrum')
			specType=spec[0].getAttribute('type')
			specPars=spec[0].getElementsByTagName('parameter')
			Ext=(True if (src.getAttribute('type')=='DiffuseSource' and not sL.psF) else False)
			sname=src.getAttribute('name')
			fixAll=(True if str(sname) in ['3FGL J0534.5+2201i','3FGL J0833.1-4511e','3FGL J1514.0-5915e','3FGL J2021.0+4031e','3FGL J2028.6+4110e'] else False)#account for sources held fixed in 3FGL analysis
			if oldNames:#if you want the same naming convention as in make1FGLxml.py and make2FGLxml.py, e.g., preceeded by an underscore and no spaces
				sn='_'
				for N in str(sname).split(' '):
					sn+=N
			varIdx=float(src.getAttribute('Variability_Index'))
			Sources[sname]={'ra':srcRA,'dec':srcDEC,'E':Ext,'stype':str(specType)}
			specOut=outputXml.createElement('spectrum')
			specOut.setAttribute('type',specType)
			spatialOut=outputXml.createElement('spatialModel')
			srcOut=outputXml.createElement('source')
			srcOut.setAttribute('name',sname)
			srcOut.setAttribute('ROI_Center_Distance',"%.2f"%dist)
			#if dist>=sL.roi[2] or dist>=sL.maxRad:
			if dist>=sL.roi[2] or dist>=sL.maxRad or fixAll:
				Sources[sname]['free']=False
				specOut.setAttribute('apply_edisp','false')#source is fixed, so never apply edisp
				for p in specPars:
					specOut.appendChild(parameter_element("0","%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
			elif dist>sL.radLim:
				if sL.var and varIdx>=72.44:
					Sources[sname]['free']=True
					specOut.setAttribute('apply_edisp',ed)
					for p in specPars:
						FreeFlag=("1" if p.getAttribute('name')==spec[0].getAttribute('normPar') else "0")
						specOut.appendChild(parameter_element("%s"%FreeFlag,"%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
				else:
					Sources[sname]['free']=False
					specOut.setAttribute('apply_edisp','false')
					for p in specPars:
						specOut.appendChild(parameter_element("0","%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
			elif float(src.getAttribute('TS_value'))>=sL.sig:
				Sources[sname]['free']=True
				specOut.setAttribute('apply_edisp',ed)
				for p in specPars:
					freeFlag=("1" if p.getAttribute('name')==spec[0].getAttribute('normPar') or (not sL.nO and p.getAttribute('free')=="1") else "0")
					#if str(p.getAttribute('name')) in normNames or (not sL.nO and str(p.getAttribute('name')) in freePars):
					specOut.appendChild(parameter_element("%s"%freeFlag,"%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
					#else:
						#specOut.appendChild(parameter_element("0","%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
			else:
				if sL.var and varIdx>=72.44:
					Sources[sname]['free']=True
					specOut.setAttribute('apply_edisp',ed)
					for p in specPars:
						FreeFlag=("1" if p.getAttribute('name')==spec[0].getAttribute('normPar') else "0")
						specOut.appendChild(parameter_element("%s"%FreeFlag,"%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
				else:
					Sources[sname]['free']=False
					specOut.setAttribute('apply_edisp','false')
					for p in specPars:
						specOut.appendChild(parameter_element("0","%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
			if Ext:
				spatial=src.getElementsByTagName('spatialModel')
				spatialOut.setAttribute('type','SpatialMap')
				spatialOut.setAttribute('map_based_integral','true')
				efile=sL.extD+spatial[0].getAttribute('file')
				spatialOut.setAttribute('file',efile)
				srcOut.setAttribute('type','DiffuseSource')
				extSrcNum+=1
				print 'Extended source %s in ROI, make sure %s is the correct path to the extended template.'%(sname,efile)
			else:
				spatialOut.setAttribute('type','SkyDirFunction')
				spatialOut.appendChild(parameter_element("0","RA","360.0","-360.0","1.0","%.4f"%srcRA))
				spatialOut.appendChild(parameter_element("0","DEC","360.0","-360.0","1.0","%.4f"%srcDEC))
				srcOut.setAttribute('type','PointSource')
				ptSrcNum+=1
			srcOut.appendChild(specOut)
			srcOut.appendChild(spatialOut)
			outputXml.documentElement.appendChild(srcOut)
	gal=outputXml.createElement('source')
	gal.setAttribute('name',GDn)
	gal.setAttribute('type','DiffuseSource')
	galspec=outputXml.createElement('spectrum')
	galspec.setAttribute('type','PowerLaw')
	galspec.setAttribute('apply_edisp','false')
	galspec.appendChild(parameter_element("1","Prefactor","10","0","1","1"))
	if sL.GIF:
		galspec.appendChild(parameter_element("1","Index","1","-1","1","0"))
	else:
		galspec.appendChild(parameter_element("0","Index","1","-1","1","0"))
	galspec.appendChild(parameter_element("0","Scale","1e6","2e1","1","100"))
	galspatial=outputXml.createElement('spatialModel')
	galspatial.setAttribute('type','MapCubeFunction')
	galspatial.setAttribute('file',GD)
	galspatial.appendChild(parameter_element("0","Normalization","1e3","1e-3","1","1"))
	gal.appendChild(galspec)
	gal.appendChild(galspatial)
	outputXml.documentElement.appendChild(gal)
	iso=outputXml.createElement('source')
	iso.setAttribute('name',ISOn)
	iso.setAttribute('type','DiffuseSource')
	isospec=outputXml.createElement('spectrum')
	isospec.setAttribute('type','FileFunction')
	isospec.setAttribute('file',ISO)
	isospec.setAttribute('apply_edisp','false')
	isospec.appendChild(parameter_element("1","Normalization","10","0.01","1","1"))
	isospatial=outputXml.createElement('spatialModel')
	isospatial.setAttribute('type','ConstantValue')
	isospatial.appendChild(parameter_element("0","Value","10","0","1","1"))
	iso.appendChild(isospec)
	iso.appendChild(isospatial)
	outputXml.documentElement.appendChild(iso)
	xmlStr=outputXml.toprettyxml(' ').splitlines(True)
	outStr=filter(lambda xmlStr: len(xmlStr) and not xmlStr.isspace(),xmlStr)
	outfile=open(sL.out,'w')
	outfile.write(''.join(outStr))
	outfile.close()
	if not sL.psF:
		print 'Added %i point sources and %i extended sources'%(ptSrcNum,extSrcNum)
		if extSrcNum>0:
			print 'If using unbinned likelihood you will need to rerun gtdiffrsp for the extended sources or rerun the makeModel function with optional argument psForce=True'
	else:
		print 'Added %i point sources, note that any extended sources in ROI were modeled as point sources becaue psForce option was set to True'%ptSrcNum
	if sL.reg:
		BuildRegion(sL,Sources)
	return

#function to cycle through the source list and add point source entries
def addSrcsFITS(sL,GD,GDn,ISO,ISOn,oldNames):
	model=open(sL.out,'w') #open file in write mode, overwrites other files of same name
	file=pyfits.open(sL.srcs) #open source list file and access necessary fields, requires LAT source catalog definitions and names
	#mask=file[1].data.field('Signif_Avg')>=sL.sig
	#data=file[1].data[mask]
	data=file['LAT_Point_Source_Catalog'].data
	extendedinfo=file['ExtendedSources'].data
	extName=extendedinfo.field('Source_Name')	
	extFile=extendedinfo.field('Spatial_Filename')
	name=data.field('Source_Name')
	Sigvals=data.field('Signif_Avg')
	VarIdx=data.field('Variability_Index')
	EName=data.field('Extended_Source_Name')
	ra=data.field('RAJ2000')
	dec=data.field('DEJ2000')
	flux=data.field('Flux_Density')
	pivot=data.field('Pivot_Energy')
	index=data.field('Spectral_Index')
	cutoff=data.field('Cutoff')
	expIndex=data.field('Exp_Index')
	spectype=data.field('SpectrumType')
	beta=data.field('Beta')
	model.write('<?xml version="1.0" ?>\n')
	model.write('<source_library title="source library">\n')
	model.write('\n<!-- Point Sources -->\n')
	step=(sL.roi[2]+sL.ER)/5. #divide ROI radius plus ExtraRadius degrees into 5 steps for ordering of sources
	i=1
	radii=[]
	ptSrcNum=0
	extSrcNum=0
	Sources={}#dictionary for sources, useful for creating region file later.
	while i<6:
		if i*step<=sL.roi[2]+sL.ER:
			radii+=[step*i]
		else:
			radii+=[sL.roi[2]+sL.ER] #just in case of rounding errors
		i+=1
	for x in radii:
		if x==sL.roi[2]+sL.ER:
			model.write('\n<!-- Sources between [%s,%s] degrees of ROI center -->\n' %(x-step,x))
		else:
			model.write('\n<!-- Sources between [%s,%s) degrees of ROI center -->\n' %(x-step,x))
		for n,f,i,r,d,p,c,t,b,TS,ei,vi,En in zip(name,flux,index,ra,dec,pivot,cutoff,spectype,beta,Sigvals,expIndex,VarIdx,EName):
			E=(True if n[-1]=='e' else False)
			dist=angsep(sL.roi[0],sL.roi[1],r,d) #check that source is within ROI radius + 10 degress of ROI center
			if r==sL.roi[0] and d==sL.roi[1]:
				dist=0.0
			if (dist<x and dist>=x-step) or (x==sL.roi[2]+10. and dist==x):
				if E and not sL.psF:
					Sources[En]={'ra':r,'dec':d,'stype':t,'E':E}
					extSrcNum+=1
					Name='<source ROI_Center_Distance="%.3f" name="%s" type="DiffuseSource">\n' %(dist,En)
				else:
					if E:#even if forcing all to point sources, use extended name
						Sources[En]={'ra':r,'dec':d,'stype':t,'E':E}
						Name='<source ROI_Center_Distance="%.3f" name="%s" type="PointSource">\n' %(dist,En)
					else:
						Sources[n]={'ra':r,'dec':d,'stype':t,'E':E}
						if oldNames:
							srcname='_'
							for N in n.split(' '):
								srcname+=N
							Name='<source ROI_Center_Distance="%.3f" name="%s" type="PointSource">\n' %(dist,srcname)
						else:
							Name='<source ROI_Center_Distance="%.3f" name="%s" type="PointSource">\n' %(dist,n)
					ptSrcNum+=1
				if t[:8]=='PowerLaw':
					fixAll=(True if n=='3FGL J0534.5+2201i' or En in ['Cygnus Cocoon','Vela X','MSH 15-52','gamma Cygni'] else False)
					spec,free=PLspec(sL,f,i,p,dist,TS,vi,fixAll)
				elif t[:9]=='PowerLaw2':#no value for flux from 100 MeV to 100 GeV in fits file
					if i!=1.:#so calculate it by integrating PowerLaw spectral model
						F=f*p**i/(-i+1.)*(1.e5**(-i+1.)-1.e2**(-i+1.))
					else:
						F=f*p*log(1.e3)
					spec,free=PL2spec(sL,F,i,dist,TS,vi)
					#spec,free=PL2spec(sL,f100,i,dist,TS,vi)
				elif t[:11]=='LogParabola':
					spec,free=LPspec(sL,f,i,p,b,dist,TS,vi)
				else:
					spec,free=COspec(sL,f,i,p,c,ei,dist,TS,vi)
				if E:
					Sources[En]['free']=free
				else:
					Sources[n]['free']=free
				if E and not sL.psF:
					#need to fix this
					efile=None
					for EXTNAME,EXTFILE in zip(extName,extFile):
						if En==EXTNAME:
							efile=sL.extD+EXTFILE
					if efile==None:
						print 'could not find a match for',En,'in the list:'
						print extName
						efile=''
					skydir='\t<spatialModel file="%s" map_based_integral="true" type="SpatialMap">\n'%(efile)
					print 'Extended source %s in ROI, make sure %s is the correct path to the extended template.'%(En,efile)
					skydir+='\t\t<parameter free="0" max="1000" min="0.001" name="Prefactor" scale="1" value="1"/>\n'
					skydir+='\t</spatialModel>\n'
				else:
					skydir='\t<spatialModel type="SkyDirFunction">\n'
					skydir+='\t\t<parameter free="0" max="360.0" min="-360.0" name="RA" scale="1.0" value="%s"/>\n' %r
					skydir+='\t\t<parameter free="0" max="90" min="-90" name="DEC" scale="1.0" value="%s"/>\n' %d
					skydir+='\t</spatialModel>\n'
				skydir+='</source>'
				(src,)=(Name+spec+skydir,)
				ptsrc=pS(src).getElementsByTagName('source')[0]
				ptsrc.writexml(model)
				model.write('\n')
	file.close() #close file
	if not sL.psF:
		print 'Added %i point sources and %i extended sources'%(ptSrcNum,extSrcNum)
		if extSrcNum>0:
			print 'If using unbinned likelihood you will need to rerun gtdiffrsp for the extended sources or rerun the makeModel function with optional argument psForce=True'
	else:
		print 'Added %i point sources, note that any extended sources in ROI were modeled as point sources becaue psForce option was set to True'%ptSrcNum
	#add galactic diffuse with PL spectrum, fix index to zero for general use, those who want it to be free can unfreeze parameter manually
	model.write('\n<!-- Diffuse Sources -->\n')
	Name='\n<source name="%s" type="DiffuseSource">\n' %GDn
	spec='\t<spectrum type="PowerLaw" apply_edisp="false">\n'
	spec+='\t\t<parameter free="1" max="10" min="0" name="Prefactor" scale="1" value="1"/>\n'
	if sL.GIF:
		spec+='\t\t<parameter free="1" max="1" min="-1" name="Index" scale="1.0" value="0"/>\n'
	else:
		spec+='\t\t<parameter free="0" max="1" min="-1" name="Index" scale="1.0" value="0"/>\n'
	spec+='\t\t<parameter free="0" max="2e2" min="5e1" name="Scale" scale="1.0" value="1e2"/>\n'
	spec+='\t</spectrum>\n'
	skydir='\t<spatialModel file="%s" type="MapCubeFunction">\n' %GD
	skydir+='\t\t<parameter free="0" max="1e3" min="1e-3" name="Normalization" scale="1.0" value="1.0"/>\n'
	skydir+='\t</spatialModel>\n'
	skydir+='</source>'
	(src,)=(Name+spec+skydir,)
	galdiff=pS(src).getElementsByTagName('source')[0]
	galdiff.writexml(model)
	model.write('\n')
	Name='<source name="%s" type="DiffuseSource">\n' %ISOn
	spec='\t<spectrum type="FileFunction" file="%s"  apply_edisp="false">\n' %ISO
	spec+='\t\t<parameter free="1" max="10" min="1e-2" name="Normalization" scale="1" value="1"/>\n'
	spec+='\t</spectrum>\n'
	skydir='\t<spatialModel type="ConstantValue">\n'
	skydir+='\t\t<parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>\n'
	skydir+='\t</spatialModel>\n'
	skydir+='</source>'
	(src,)=(Name+spec+skydir,)
	iso=pS(src).getElementsByTagName('source')[0]
	iso.writexml(model)	
	model.write('\n</source_library>')
	model.close()
	if sL.reg:
		BuildRegion(sL,Sources)
	return

def BuildRegion(sL,Sources):
	myreg=open(sL.regFile,'w')#note that this will overwrite previous region files of the same name
	myreg.write('# Region File format: DS9 version ?')#I don't actually know, but I think it's one of the later ones, need to verify
	myreg.write('\n# Created by make3FGLxml.py')
	myreg.write('\nglobal font="roman 10 normal" move =0')
	for k in Sources.keys():
		src=Sources[k]
		#get color based on if the source is free or not
		color=('green' if src['free'] else 'magenta')
		if src['E']:#if the source is extended, always have the point be a "big" box
			myreg.write('\nJ2000;point(%.3f,%.3f) # point = box 18 color = %s text={%s}'%(src['ra'],src['dec'],color,k))
		else:#if the source is a point source, choose the point type based on spectral model
			ptype=('cross' if src['stype']=='PLSuperExpCutoff' else 'diamond' if src['stype']=='LogParabola' else 'circle')
			myreg.write('\nJ2000;point(%.3f,%.3f) # point = %s 15 color = %s text={%s}'%(src['ra'],src['dec'],ptype,color,k))
	myreg.close()
	return

def PLspec(sL,f,i,p,dist,TS,vi,fixAll):
	fscale=int(floor(log10(f)))
	if (dist>sL.roi[2] or fixAll) or (dist>sL.maxRad) or (dist>sL.radLim and (vi<72.44 or not sL.var)) or (TS<sL.sig and (vi<72.44 or not sL.var)):
		spec='\t<spectrum type="PowerLaw" apply_edisp="false">\n'
	else:
		spec='\t<spectrum type="PowerLaw" apply_edisp="%s">\n'%sL.ed
	spec+='\t<!-- Source is %s degrees away from ROI center -->\n' %dist
	if dist>sL.roi[2] or fixAll: #if beyond ROI, shouldn't attempt to fit parameters
		if fixAll:
			spec+='\t<!-- Source parameters were held fixed in 3FGL analysis, free at your own discretion -->\n'
		else:
			spec+='\t<!-- Source is outside ROI, all parameters should remain fixed -->\n'
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
		free=False
	elif(dist>sL.radLim):
		if dist>sL.maxRad:
			spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
			spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=False
		elif vi<72.44 or not sL.var:
			spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
			spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=False
		else:
			spec+='\t<!-- Source is outside specified radius limit of %s but variability index %.2f is greater than 72.44 and varFree set to True-->\n'%(sL.radLim,vi)
			free=True
			spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
	elif(TS<sL.sig):
		if vi<72.44 or not sL.var:
			spec+='\t<!-- Source signficance %.1f is less than specified minimum for a free source of %s -->\n'%(TS,sL.sig)
			spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=False
		else:
			spec+='\t<!-- Source significane %.1f is less than specified minimum for a free source of %s but variability index %.2f is greater than 72.44 and varFree is set to True -->\n'%(TS,sL.sig,vi)
			spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=True
		spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
	else:
		spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		free=True
		if sL.nO:
			spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
		else:
			spec+='\t\t<parameter free="1" max="10.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="Scale" scale="1.0" value="%f"/>\n' %p
	spec+='\t</spectrum>\n'
	return spec,free

def PL2spec(sL,F,i,dist,TS,vi):
	fscale=int(floor(log10(F)))
	if dist>sL.roi[2] or (dist>sL.maxRad) or (dist>sL.radLim and (vi<72.44 or not sL.var)) or (TS<sL.sig and (vi<72.44 or not sL.var)):
		spec='\t<spectrum type="PowerLaw2" apply_edisp="false">\n'
	else:
		spec='\t<spectrum type="PowerLaw2" apply_edisp="%s">\n'%sL.ed
	spec+='\t<!-- Source is %s degrees away from ROI center -->\n' %dist
	if(dist>sL.roi[2]): #if beyond ROI, shouldn't attempt to fit parameters
		spec+='\t<!-- Source is outside ROI, all parameters should remain fixed -->\n'
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Integral" scale="1e%i" value="%s"/>\n'%(fscale,F/10**fscale)
		spec+='\t\t<parameter free="0" max="10" min="0" name="Index" scale="-1.0" value="%s"/>\n' %i
		free=False
	elif(dist>sL.radLim):
		if dist>sL.maxRad:
			spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
			spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Integral" scale="1e%i" value="%s"/>\n'%(fscale,F/10**fscale)
			free=False
		elif vi<72.44 or not sL.var:
			spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
			spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Integral" scale="1e%i" value="%s"/>\n'%(fscale,F/10**fscale)
			free=False
		else:
			spec+='\t<!-- Source is outside specified radius limit of %s but variability index %.2f is greater than 72.44 and varFree is set to True -->\n'%(sL.radLim,vi)
			spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Integral" scale="1e%i" value="%s"/>\n'%(fscale,F/10**fscale)
			free=True
		spec+='\t\t<parameter free="0" max="10" min="0" name="Index" scale="-1.0" value="%s"/>\n' %i
	elif(TS<sL.sig):
		if vi<72.44 or not sL.var:
			spec+='\t<!-- Source significance %.1f is less than specified minimum for a free source of %s -->\n'%(TS,sL.sig)
			spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Integral" scale="1e%i" value="%s"/>\n'%(fscale,F/10**fscale)
			free=False
		else:
			spec+='\t<!-- Source significance %.1f is less than specified minimum for a free source of %s but variability index %.2f is greater than 72.44 and varFree is set to True -->\n'%(TS,sL.sig,vi)
			spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Integral" scale="1e%i" value="%s"/>\n'%(fscale,F/10**fscale)
			free=True
		spec+='\t\t<parameter free="0" max="10" min="0" name="Index" scale="-1.0" value="%s"/>\n' %i
	else:
		free=True
		spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Integral" scale="1e%i" value="%s"/>\n'%(fscale,F/10**fscale)
		if sL.nO:
			spec+='\t\t<parameter free="0" max="10" min="0" name="Index" scale="-1.0" value="%s"/>\n' %i
		else:
			spec+='\t\t<parameter free="1" max="10" min="0" name="Index" scale="-1.0" value="%s"/>\n' %i
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="LowerLimit" scale="1" value="1e2"/>\n'
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="UpperLimit" scale="1" value="1e5"/>\n'
	spec+='\t</spectrum>\n'
	return spec,free

def COspec(sL,f,i,p,c,ei,dist,TS,vi):
	f*=exp((p/c)**ei)
	fscale=int(floor(log10(f)))
	if dist>sL.roi[2] or (dist>sL.maxRad) or (dist>sL.radLim and (vi<72.44 or not sL.var)) or (TS<sL.sig and (vi<72.44 or not sL.var)):
		spec='\t<spectrum type="PLSuperExpCutoff" apply_edisp="false">\n'
	else:
		spec='\t<spectrum type="PLSuperExpCutoff" apply_edisp="%s">\n'%sL.ed
	spec+='\t<!-- Source is %s degrees away from ROI center -->\n' %dist
	i=(i if i>=0 else 2.)#some pulsars with index1 < 0 assuming standard convention, means rising counts spectrum at low E, odd
	if(dist>sL.roi[2]): #if beyond ROI, shouldn't attempt to fit parameters
		spec+='\t<!-- Source is outside ROI, all parameters should remain fixed -->\n'
		free=False
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="Index1" scale="-1.0" value="%s"/>\n' %i
		if c<=1e5:
			spec+='\t\t<parameter free="0" max="1e5" min="1e1" name="Cutoff" scale="1.0" value="%f"/>\n'%c
		else:
			spec+='\t\t<parameter free="0" max="%.2e" min="1e1" name="Cutoff" scale="1.0" value="%f"/>\n'%(2.*c,c)
	elif(dist>sL.radLim):
		if dist>sL.maxRad:
			spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
			spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=False
		elif vi<72.44 or not sL.var:
			spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
			spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=False
		else:
			spec+='\t<!-- Source is outside specified radius limit of %s but variability index %.2f is greater than 72.44 and varFree is set to True -->\n'%(sL.radLim,vi)
			spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=True
		spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="Index1" scale="-1.0" value="%s"/>\n' %i
		if c<=1e5:
			spec+='\t\t<parameter free="0" max="1e5" min="1e1" name="Cutoff" scale="1.0" value="%f"/>\n'%c
		else:
			spec+='\t\t<parameter free="0" max="%.2e" min="1e1" name="Cutoff" scale="1.0" value="%f"/>\n'%(2.*c,c)
	elif(TS<sL.sig):
		if vi<72.44 or not sL.var:
			spec+='\t<!-- Source significance %.1f is less than specified minimum for a free source of %s -->\n'%(TS,sL.sig)
			spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=False
		else:
			spec+='\t<!-- Source significance %.1f is less than specified minimum for a free source of %s but variability index %.2f is greater than 72.44 and varFree is set to True -->\n'%(TS,sL.sig,vi)
			spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=True
		spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="Index1" scale="-1.0" value="%s"/>\n' %i
		if c<=1e5:
			spec+='\t\t<parameter free="0" max="1e5" min="1e1" name="Cutoff" scale="1.0" value="%f"/>\n'%c
		else:
			spec+='\t\t<parameter free="0" max="%.2e" min="1e1" name="Cutoff" scale="1.0" value="%f"/>\n'%(2.*c,c)
	else:
		free=True
		spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		if sL.nO:
			spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="Index1" scale="-1.0" value="%s"/>\n' %i
			if c<=1e5:
				spec+='\t\t<parameter free="0" max="1e5" min="1e1" name="Cutoff" scale="1.0" value="%f"/>\n'%c
			else:
				spec+='\t\t<parameter free="0" max="%.2e" min="1e1" name="Cutoff" scale="1.0" value="%f"/>\n'%(2.*c,c)
		else:
			spec+='\t\t<parameter free="1" max="10.0" min="0.0" name="Index1" scale="-1.0" value="%s"/>\n' %i
			if c<=1e5:
				spec+='\t\t<parameter free="1" max="1e5" min="1e1" name="Cutoff" scale="1.0" value="%f"/>\n'%c
			else:
				spec+='\t\t<parameter free="0" max="%.2e" min="1e1" name="Cutoff" scale="1.0" value="%f"/>\n'%(2.*c,c)
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="Scale" scale="1.0" value="%f"/>\n'%p
	spec+='\t\t<parameter free="0" max="5" min="0" name="Index2" scale="1.0" value="%f"/>\n'%ei
	spec+='\t</spectrum>\n'
	return spec,free

def LPspec(sL,f,i,p,b,dist,TS,vi):
	fscale=int(floor(log10(f)))
	if dist>sL.roi[2] or (dist>sL.maxRad) or (dist>sL.radLim and (vi<72.44 or not sL.var)) or (TS<sL.sig and (vi<72.44 or not sL.var)):
		spec='\t<spectrum type="LogParabola" apply_edisp="false">\n'
	else:
		spec='\t<spectrum type="LogParabola" apply_edisp="%s">\n'%sL.ed
	spec+='\t<!-- Source is %s degrees away from ROI center -->\n' %dist
	if(dist>sL.roi[2]): #if beyond ROI, shouldn't attempt to fit parameters
		spec+='\t<!-- Source is outside ROI, all parameters should remain fixed -->\n'
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="norm" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="alpha" scale="1.0" value="%s"/>\n' %i
		spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="beta" scale="1.0" value="%s"/>\n'%b
		free=False
	elif(dist>sL.radLim):
		if dist>sL.maxRad:
			spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
			spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="norm" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=False
		elif vi<72.44 or not sL.var:
			spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
			spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="norm" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=False
		else:
			spec+='\t<!-- Source is outside specified radius limit of %s but variability index %.2f is greater than 72.44 and varFree is set to True -->\n'%(sL.radLim,vi)
			spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="norm" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=True
		spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="alpha" scale="1.0" value="%s"/>\n' %i
		spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="beta" scale="1.0" value="%s"/>\n'%b
	elif(TS<sL.sig):
		if vi<72.44 or not sL.var:
			spec+='\t<!-- Source significance %.1f is less than specified minimum for a free source of %s -->\n'%(TS,sL.sig)
			spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="norm" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=False
		else:
			spec+='\t<!-- Source significance %.1f is less than specified minimum for a free source of %s but variability index %.2f is greater than 72.44 and varFree is set to True -->\n'%(TS,sL.sig,vi)
			spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="norm" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
			free=True
		spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="alpha" scale="1.0" value="%s"/>\n' %i
		spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="beta" scale="1.0" value="%s"/>\n'%b
	else:
		free=True
		spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="norm" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		if sL.nO:
			spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="alpha" scale="1.0" value="%s"/>\n' %i
			spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="beta" scale="1.0" value="%s"/>\n'%b
		else:
			spec+='\t\t<parameter free="1" max="5.0" min="0.0" name="alpha" scale="1.0" value="%s"/>\n' %i
			spec+='\t\t<parameter free="1" max="10.0" min="0.0" name="beta" scale="1.0" value="%s"/>\n'%b
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="Eb" scale="1.0" value="%s"/>\n'%p
	spec+='\t</spectrum>\n'
	return spec,free

#this function searches the header of the ft1 to find the Position keyword and extract the ra and dec values
def getPos(ft1):
	file=pyfits.open(ft1)
	num=file[1].header['NDSKEYS']
	header=file[1].header
	right='POS(RA,DEC)'
	i=1
	keynum=0
	while i<=num:  #this step is necessary since it is not clear that the POS key word will have the same number always
		word='DSTYP%i' %i
		test=file[1].header[word]
		if(test==right):
			keynum=i
			i=num
		i+=1
	if(keynum==0):  #DSKEYS start numbering at 1, if this value hasn't been updated, KEYword doesn't exist
		print 'Error: No position keyword found in fits header (assuming position is RA and DEC.  Exiting...'
		exit()
	keyword='DSVAL%i' %keynum
	try:
		ra,dec,rad=header[keyword].strip('CIRCLE()').split(',') #gets rid of the circle and parenthesis part and splits around the comma
		float(ra)
	except:
		ra,dec,rad=header[keyword].strip('circle()').split(',')
	file.close()
	return float(ra),float(dec),float(rad)
	
#calculates the angular separation between two points on the sky
def angsep(ra1,dec1,ra2,dec2):
	ra1*=d2r
	dec1*=d2r
	ra2*=d2r
	dec2*=d2r
	diffCosine=cos(dec1)*cos(dec2)*cos(ra1-ra2)+sin(dec1)*sin(dec2)
	dC='%.10f'%diffCosine#when the source is right at the center of the roi python sometimes adds extraneous digits at the end of the value i.e. instead of 1.0
	#it returns 1.0000000000000024, which throws an error with the acos function
	return acos(float(dC))/d2r #returns values between 0 and pi radians

#Check if a given file exists or not
def fileCheck(file):
	if (not os.access(file,os.F_OK)):
		return 0
	return 1

#copied from Damien's macro
def parameter_element(free, name, maximum, minimum, scale, value):
    """Create an XML document parameter description element"""
    impl = minidom.getDOMImplementation()
    xmldoc_out = impl.createDocument(None,None,None)
    parameter = xmldoc_out.createElement('parameter')
    parameter.setAttribute('free', str(free))
    parameter.setAttribute('name', str(name))
    parameter.setAttribute('max', str(maximum))
    parameter.setAttribute('min', str(minimum))
    parameter.setAttribute('scale', str(scale))
    parameter.setAttribute('value', str(value))
    return parameter

def mybool(Input):
	return {'True':True,'False':False,'T':True,'F':False,'t':True,'f':False,'TRUE':True,'FALSE':False,"true":True,"false":False,"1":True,"0":False}.get(Input)

def cli():
	import argparse
	
	helpString="Creates an xml model from the 3FGL catalog (FITS or xml version) for a specific ROI,\
		    coordinates of the ROI center are taken from an input event file,\
		    the radius for including sources is 10 degrees beyond the extraction radius used in the event file,\
		    sources with free parameters within the original extraction radius are chosen based on nearness to center, significance, and variability."
	parser=argparse.ArgumentParser(description=helpString)
	parser.add_argument("catalog",type=str,help="Catalog file to use, can be FITS or xml.")
	parser.add_argument("ev",type=str,help="Event file with ROI information in header.")
	parser.add_argument("-o","--outputxml",type=str,default='mymodel.xml',help="Name of output xml file, is set to overwrite files of same name.")
	parser.add_argument("-G","--galfile",type=str,default='$(FERMI_DIR)/refdata/fermi/galdiffuse/gll_iem_v06.fits',help="Name and location of Galactic diffuse model to use, will default to 3FGL model.")
	parser.add_argument("-g","--galname",type=str,default='gll_iem_v06',help="Name of Galactic diffuse component in output model, will default to gll_iem_v06.")
	parser.add_argument("-I","--isofile",type=str,default='$(FERMI_DIR)/refdata/fermi/galdiffuse/iso_P8R2_SOURCE_V6_v06.txt',help="Name of isotropic diffuse template for output model, will default to 3FGL model.")
	parser.add_argument("-i","--isoname",type=str,default='iso_P8R2_SOURCE_V6_v06',help="Name of isotropic diffuse component in output model, will default to iso_source_v05.")
	parser.add_argument("-N","--normsonly",type=mybool,default=False,help="Flag to only let the normalizations of parameters be free, default is False.",nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-e","--extDir",type=str,default='',help="Path to directory with LAT extended source templates, will default to STs default.")#need to figure out what that is
	parser.add_argument("-r","--radLim",type=float,default=-1.,help="Radius, in degrees, from ROI center beyond which all source parameters should be fixed, will default to selection radius.")
	parser.add_argument("-R","--maxRad",type=float,default=None,help="Absolute maximum radius, in degrees, from ROI center beyond which all source parameters should be fixed, even variable sources will not be freed beyond this radius, defaults to radLim value.")
	parser.add_argument("-ER","--ExtraRad",type=float,default=10.,help="Radius beyond event file ROI out to which sources will be included in the model with all parameters fixed, default is 10, good for analyses starting around a few hundred MeV, can be decreased for high energy only fits.")
	parser.add_argument("-s","--sigFree",type=float,default=5.,help="Average significance below which all source parameters are fixed, defaults to 5.  Note, if using the 3FGL catalog xml file as input, this is actually a cut on TS, so adjust accordingly.")
	parser.add_argument("-v","--varFree",type=mybool,default=True,help="Flag to set normalization of significantly variable sources, even if source is beyond radius limit or below TS limit, default is True.",choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-p","--psForce",type=mybool,default=False,help="Flag to cast extended sources as point sources, default is False.",nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-m","--makeRegion",type=mybool,default=True,help="Flag to create ds9 region file as well as the xml model, default is True.",choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-GIF","--GIndexFree",type=mybool,default=False,help="Flag to use a power-law modification to the Galactic diffuse model spectrum and have the index be free, default is False.",nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-ED","--edisp",type=mybool,default=False,help="Flag to turn on energy dispersion for free point and extended sources, never for diffuse backgrounds, default is False.",nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-wd","--writeDir",type=str,default='',help="Directory to write the output ds9 region file in if not the current working directory or if you are specifying the full path to the newly made XML file.")
	parser.add_argument("-ON","--oldNames",type=mybool,default=False,help="Flag to use the make2FLGxml style naming convention, underscore before name and no spaces, default is False.",nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-P7","--pass7",type=mybool,default=False,help="Flag to say you're making a model for analysis of P7 data, default is False.  The only reason to use this is to switch the defaults for the diffuse components.",nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	
	args=parser.parse_args()
	
	if args.pass7:
		args.galfile='$(FERMI_DIR)/refdata/fermi/galdiffuse/gll_iem_v05_rev1.fits'
		args.galname='gll_iem_v05_rev1'
		args.isofile='$(FERMI_DIR)/refdata/fermi/galdiffuse/iso_source_v05.txt'
		args.isoname='iso_source_v05'
	
	print(args.catalog,args.ev,args.outputxml)
	print(args.galfile,args.galname,args.isofile,args.isoname,args.normsonly,args.extDir,args.radLim,args.maxRad,args.ExtraRad,args.sigFree,args.varFree,args.psForce,args.makeRegion,args.GIndexFree,args.edisp,args.writeDir,args.oldNames)

	sL=srcList(args.catalog,args.ev,args.outputxml)
	sL.makeModel(args.galfile,args.galname,args.isofile,args.isoname,args.normsonly,args.extDir,args.radLim,args.maxRad,args.ExtraRad,args.sigFree,args.varFree,args.psForce,args.makeRegion,args.GIndexFree,args.edisp,args.writeDir,args.oldNames)
	
	
if __name__=='__main__': cli()
