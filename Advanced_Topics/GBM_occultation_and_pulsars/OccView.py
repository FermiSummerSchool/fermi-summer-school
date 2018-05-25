#Sample python code to read and plot GBM Earth Occultation FITS files
#Authors J.M. Burgess and C.A. Wilson-Hodge 2014 Jul 18
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from numpy import array, sqrt, sum, mean, logical_and


class OccView(object):
    
    def __init__(self,data):
        
        data = fits.open(data) #read in FITS file
        
        self.fluxExt = data[1].data
        
        #if we have a single step file, convert the times from MET to MJD 
        if self.fluxExt['TIME'][0]>200000000.:  
           newtime = 51910+0.000742870370370+self.fluxExt['TIME']/86400.
           self.fluxExt['TIME'] = newtime
           
        #Create energy edges
        numEdges = len(data[2].data['E_MIN'])
        eneEdge = data[2].data
        #for i in range(numEdges):
        #    
        #    eneEdge.append(data[2].data['E_MIN'][i])
	#    if i==numEdges:
	#       eneEdge.append(data[2].data['E_MAX'][i])
        #       	   
        #
        self.eneEdge = array(eneEdge)
        self.figCount = 1
        
        
    def _Quality(self,q):
        
        if q=="good":
            q=0
        elif q=="bad":
            q=1
        else:
            print "Quality flag must be  [good/bad]"
            print "Defaulting to [good]"
            q=0
        
        return self.fluxExt["QUALITY"]==q
    
    
    def _FluxPlot(self,ax,time,flux,err,color):
        
        ax.errorbar(time,flux,yerr=err,fmt="-",color=color) #plot fluxes and errors
        ax.get_figure().autofmt_xdate()
        #ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=45)
        ax.set_xlabel("Time (MJD)")
        ax.set_ylabel("Flux (photons/cm**2/s)")
        print
        print "....Plotting flux"
        print
        
        
        
    def _BinData(self,flux,err,tstart,tstop,dt,quality):
        
        #Bin data in bins of user defined duration
        
        #First get the truth tables
        #tt1=self.fluxExt["TSTART"]>=tstart
        #tt2=self.fluxExt["TSTOP"]<=tstop
        #tt = logical_and(tt1,tt2)
        #time = self.fluxExt["TSTART"][tt]
        
        currentTime = tstart
        meanTime = []
        binnedFlux = []
        binnedError = []
        
        while(currentTime<tstop):
            
            
            tt1=self.fluxExt["TIME"][self._Quality(quality)]>=currentTime
            tt2=self.fluxExt["TIME"][self._Quality(quality)]<=currentTime+dt
            
            thisBin = logical_and(tt1,tt2)
            #added next line - want to use mid point times for each measurement (TSTART+TSTOP)/2
            meanTimes = array(map(lambda x,y:mean([x,y]),self.fluxExt["TIME"][self._Quality(quality)][thisBin],self.fluxExt["TIME"][self._Quality(quality)][thisBin]))
            fluxes = flux[thisBin] #deleted /dt
            errors = err[thisBin]  #deleted /dt
	    mask = errors[:,0] !=0 
            sumFlux = sum(fluxes[mask,:]/errors[mask,:]**2)/sum(1.0/errors]mask,:]**2)  #weighted average
              sumErr  = 1.0/sqrt(sum(1.0/errors[mask,:]**2))              #error on average
              binnedFlux.append(sumFlux)
              binnedError.append(sumErr)
            #meanTime.append(mean([currentTime,currentTime+dt]))
              meanTime.append(mean(meanTimes))  #Take the mean of the times included in the average
              currentTime+=dt      
        
        
        meanTime = array(meanTime)
        binnedFlux = array(binnedFlux)
        return [meanTime,binnedFlux,binnedError]
    
    
            
    def _PlotFluxBinned(self,ax,tstart,tstop,dt,channel,color,quality):
        #plot binned data
        flux = self.fluxExt["RATE"][:,channel][self._Quality(quality)]
        
        err = self.fluxExt["ERROR"][:,channel][self._Quality(quality)]
        meanTime, flux, err = self._BinData(flux,err,tstart,tstop,dt,quality)
        
        
        self._FluxPlot(ax,meanTime,flux,err,color)
        
    def _PlotFlux(self,ax,tstart,tstop,channel,color,quality):
        #plot raw data (either daily or single step depending on FITS file)
        flux = self.fluxExt["RATE"][:,channel][self._Quality(quality)]
        
        err = self.fluxExt["ERROR"][:,channel][self._Quality(quality)]
        
        tt1=self.fluxExt["TIME"][self._Quality(quality)]>=tstart
        tt2=self.fluxExt["TIME"][self._Quality(quality)]<=tstop
        
        tt = logical_and(tt1,tt2)
        flux = flux[tt]
        err=err[tt]
        
        meanTime = array(map(lambda x,y:mean([x,y]),self.fluxExt["TIME"][self._Quality(quality)][tt],self.fluxExt["TIME"][self._Quality(quality)][tt]))
        
        self._FluxPlot(ax,meanTime,flux,err,color)
    
    
    def GetFluxes(self,tstart,tstop,channel,quality="good"):
        #Return time, flux, and error for specified channel and time range, with quality good
        flux = self.fluxExt["RATE"][:,channel][self._Quality(quality)]
        
        err = self.fluxExt["ERROR"][:,channel][self._Quality(quality)]
        
        tt1=self.fluxExt["TIME"][self._Quality(quality)]>=tstart
        tt2=self.fluxExt["TIME"][self._Quality(quality)]<=tstop
        
        tt = logical_and(tt1,tt2)
        flux = flux[tt]
        err=err[tt]
        
        meanTime = array(map(lambda x,y:mean([x,y]),self.fluxExt["TIME"][self._Quality(quality)][tt],self.fluxExt["TIME"][self._Quality(quality)][tt]))
       
        data = {"time":meanTime,"flux":flux,"error":err}
        return data
    
    def GetBinnedFluxes(self,tstart,tstop,dt,channel,quality="good"):
        #Return binned times, fluxes, errors, within user specified time range, bin duration, and energy channel
        flux = self.fluxExt["RATE"][:,channel][self._Quality(quality)]
        err = self.fluxExt["ERROR"][:,channel][self._Quality(quality)]
        meanTime, flux, err = self._BinData(flux,err,tstart,tstop,dt,quality)
        
        data = {"time":meanTime,"flux":flux,"error":err}
        return data
        
    def PlotBinnedFluxes(self,tstart,tstop,dt,channels,quality="good",save=None):
        #plot binned fluxes for multiple channels if requested for specified time range and time bin duration
        fig = plt.figure(100+self.figCount)
        ax = fig.add_subplot(111)
        
        
        colors = ["red","green","blue","orange","cyan","k"]
        energies=[]
        if type(channels) == list:
            i=0
            
            for ch in channels:
                self._PlotFluxBinned(ax,tstart,tstop,dt,ch,colors[i],quality)
                edge = self.eneEdge[ch]
                s="%.1f-%.1f keV"%(edge[0],edge[1])
                energies.append(s)
                i+=1
            
        else:
            self._PlotFluxBinned(ax,tstart,tstop,dt,channels,colors[0],quality)
            edge = self.eneEdge[channels]
            s="%.1f-%.1f keV"%(edge[0],edge[1])
            energies.append(s)
        
        ax.legend(energies)
        ax.text(.05,.9,"Days binned: %.3f"%dt,transform=ax.transAxes)
        self.figCount+=1
        if save:
            fig.savefig(save)

        plt.show() 
    
    def PlotFluxes(self,tstart,tstop,channels,quality="good",save=None):
        #plot raw fluxes (daily or single step) for specificed time range and multiple energy channels
        fig = plt.figure(1000+self.figCount)
        ax = fig.add_subplot(111)
        
        
        colors = ["red","green","blue","orange","cyan","k"]
        energies=[]
        if type(channels) == list:
            i=0
            
            for ch in channels:
                self._PlotFlux(ax,tstart,tstop,ch,colors[i],quality)
                edge = self.eneEdge[ch]
                s="%.1f-%.1f keV"%(edge[0],edge[1])
                energies.append(s)
                i+=1
            
        else:
            self._PlotFlux(ax,tstart,tstop,channels,colors[0],quality)
            edge = self.eneEdge[channels]
            s="%.1f-%.1f keV"%(edge[0],edge[1])
            energies.append(s)
        
        ax.legend(energies)
        self.figCount+=1
        if save:
            fig.savefig(save)
        plt.show() 
