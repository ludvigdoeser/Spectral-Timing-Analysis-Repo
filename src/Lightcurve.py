#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from ImportPackages import *
from AuxiliaryFunctions import * 
from AuxiliaryFunctions import extract_fits

class lightcurve:
    """
    Make a light curve object from fits-file.
    
    **Parameters**:

    `filename`: string 
        The path to the .fits-file. If an empty string is given, a lightcurve object will 
        be initiated without values for its attributes.
        
    `keywords`: list of strings    
        Keywords (and corresponding values) apart from the ones mentioned below to return.    
        All data-keys are always automatically returned (for a lightcurve, these are: {"TIME", "RATE", "ERROR","FRACEXP"}).   
        Note also that the header-keys {"CONTENT","OBJECT"} are returned if they exist.    
    
    `p`: boolean    
        If True, print the headers of the fits-file.
        
    **Attributes**:

    `t`: np.ndarray    
        Time vector (in seconds).
        
    `rate`: np.ndarray    
        Rate vector (in counts/s).
        
    `err`: np.ndarray    
        Error vector in rate (in counts/s).
        
    `fracexp`: np.ndarray    
        Fractional exposure for each bin (in percent).
    
    `object`: string    
        Object for which we have constructed the lightcurve object.
        
    `dt`: float    
        Time resolution (in seconds).
        
    `R`: np.float    
        Mean count rate.
        
    `N`: int    
        Number of data points in observation.
    
    `Fvar`: float    
        Fractional variance amplitude (=rms) computed from the light curve.
        
    `Emin, Emax`: floats    
        Min and max energies of used energy band.
        
    `deltaE`: float    
        Energy range between min and max channels used.
    
    `E_mean`: float    
        Energy mean over the range between min and max channels used.
    """
      
    def __init__(self, filename, keywords=[], p=0):
        if isinstance(filename,str):
            
            if len(filename) != 0:
                # Initiate light curve object from .fits-file
                self.t, self.rate, self.err, self.fracexp, self.object = self.extract_lc(filename, keywords, p) #self.channels, self.minchan, self.maxchan
                self.dt = self.t[1]-self.t[0] # should not (but could) be a gap here... 
                self.R = np.mean(self.rate)
                self.N = len(self.t)  
                self.Emin, self.Emax = self.find_energyband(filename)
                self.deltaE = (self.Emax-self.Emin)/2
                self.Emean = (self.Emax+self.Emin)/2
                self.Fvar = self.Fvar_from_lc(m=self.N//100) #use 100 segments as default; can be changed.
                print('Light curve object for {} in Eband = {}-{} keV created.'.format(self.object,self.Emin,self.Emax))
                print('With parameters: N = {}, dt = {:.2}, R = {:.4}, and Fvar = {:.4}.'.format(self.N,self.dt,self.R,self.Fvar))
                print('-----------------------------------------------------------------------------------------------------------\n')
           
            else:
                # Initiate light curve object without any attribute values'
                self.t, self.rate, self.err, self.object = [], [], [], ''
                self.R, self.dt, self.N, self.Emin, self.Emax, self.deltaE, self.Emean, self.Fvar = 0,0,0,0,0,0,0,0
            
    def extract_lc(self, filename, keywords, p):
        """
        Extract light curve from .fits-file. For more info, see class doc.
        """
        
        print('-----------------------------------------------------------------------------------------------------------')
        print('               Importing lightcurve from f = {}'.format(filename))
        print('-----------------------------------------------------------------------------------------------------------')
        data = extract_fits(filename, keywords, p=p)

        t = np.array(data['TIME'])
        t -= t[0] #so that it starts from 0 
        rate = np.array(data['RATE'])
        err = np.array(data['ERROR'])
        fracexp = np.array(data['FRACEXP'])
        
        # Check for NaN in rate and remove such elements
        if len(np.argwhere(np.isnan(rate))) > 0:
            t = t[np.isfinite(rate)]
            err = err[np.isfinite(rate)]
            fracexp = fracexp[np.isfinite(rate)]
            rate = rate[np.isfinite(rate)]
        
        obj = data['OBJECT'].replace('_','')
        
        if len(data) > 5:
            keys = [key for key,val in data.items()]   
            for i in range(5,len(data)):
                setattr(self, keys[i].replace('-',''), data[keys[i]])
        
        return t, rate, err, fracexp, obj #, channels, minchan, maxchan
    
    def plot(self,m=0,n=0):
        """
        Display the full lightcurve.

        **Parameters**:

        `m`: int, optional, default = 0    
            Number of time bins per segment.
 
        `n`: int, optional, default: 0    
            The n:th segment will be displayed in the full light curve.
        """
        
        # Parameters 
        if m!=0:
            K = int(np.floor(self.N/m))
            print(r"Number of line segments of \approx {:.0f}s will be: ".format((self.t[1]-self.t[0])*m),K)
        else: 
            K = 1

        # Plot 
        ax = standard_plot()
        #plt.errorbar(t_seg,rate_seg,err_seg,fmt='o',markersize=2)
        plt.plot(self.t,self.rate,'-')
        ax.axhline(self.R,color='r',label='R = {:.0f}'.format(self.R))
        if m > 0:
            ax.axvline(self.t[n*m],color='k')
            ax.axvline(self.t[(n+1)*m],color='k')
            ax.axvspan(self.t[n*m], self.t[(n+1)*m], alpha=0.4, color='k',label='Line segment')

        plt.legend(loc='upper right')
        plt.xlabel('Time [$s$]')
        plt.ylabel('Rate [$c/s$]')
        plt.title('{} in {}-{} keV'.format(self.object,self.Emin,self.Emax))
        plt.show()
        
    def rebin(self,dt=1,force=False):
        """
        Rebin light curve into new time resolution. Note: inefficient... 
        Better to extract new light curves using XSELECT.

        **Parameters**:
        
        `dt`: float, optional, default: 1 (sec)    
            New time resolution.
            
        `force`: boolean, optional, default: False    
            If True, don't ask for permission to make time resolution changes
        """

        # New time vector
        eps = 10e-15 #to include the final point
        t = np.arange(self.t[0], self.t[-1]+eps,dt)
        self.N = len(t)
        s = self.dt/dt #scale factor

        if not force:
            still_update = input('You want to change the time resolution from dt = {:.3} (with N = {}) to dt = {:.3}. This will lead to a loss of information. \n Continue [y/n]? '.format(self.dt,self.N,float(dt)))
        else:
            still_update = 'y'
        
        if still_update == 'y':

            # Mid-points
            t_mid = [np.mean([t[i],t[i+1]]) for i in range(0,len(t)-1)]

            # Determine what old t-values goes into what new t-bin
            digitized = np.digitize(self.t, t)
            # New t, rate, err and fracexp vectors 

            rate, err, fracexp, t = [], [], [], []
            for i in range(1,self.N):
                if len(self.t[digitized == i])!=0:
                    rate.append(s*self.rate[digitized == i].sum())
                    err.append(s*np.sqrt(np.sum((self.err[digitized == i])**2)))
                    t.append(t_mid[i-1])

            # Update 
            self.t = np.array(t)
            self.rate = np.array(rate)
            self.err = np.array(err)
            self.dt = dt
            self.N = len(t)
            self.R = np.mean(self.rate)
            self.Fvar = self.Fvar_from_lc(m=self.N//100) #use 100 segments as default; can be changed.
            print('Time resolution has been changed.')
        else:
            print('Did not change the time resolution.')
        
        return self
        
    def extract_seg(self,m,n=0,bins_back=0,to_print=False,to_plot=False):
        """
        Extract (and potentially display) a light curve segment.

        **Parameters**:
        
        `m`: int    
            Number of time bins per segment. 

        `n`: int, optional, default: 0    
            The start (end) of the segment will be bin number m*n (m*(n+1)) 
            
        `bins_back`: int, optional, default: 0 
            Number of bins the segment will be moved backwards in the light curve. 
            
        `to_print`: boolean, optional, default: False 
            If True, prints out average count rate, totalt counts, and total time of segment. 
            
        `to_plot`: boolean, optional, default: False 
            If True, plot the segment. 

        **Returns**: 
        
        `t_seg, rate_seg, err_seg`: np.ndarrays    
            Time, rate, and error vectors for the segment.

        `N_gamma, R, T`: floats    
            Number of photons (in counts), mean count rate (in counts/s) and length (in seconds) respectively. 
        """

        # Pick out segment
        start = m*n - bins_back
        stop = m*(n+1) - bins_back
        
        t_seg = self.t[start:stop]
        rate_seg = self.rate[start:stop]
        err_seg = self.err[start:stop]

        # Relevant parameters
        T = t_seg[-1]-t_seg[0]
        R = np.mean(rate_seg)
        N_gamma = T*R

        # Print
        if to_print:
            print('R = average count rate = ',R)
            print('N_\gamma = total counts = ',N_gamma)
            print('T = total time = ',T) 

        # Plot
        if to_plot:
            ax = standard_plot()
            #plt.errorbar(t_seg,rate_seg,err_seg,fmt='o',markersize=2)
            plt.plot(t_seg,rate_seg,'-')
            ax.axhline(np.mean(rate_seg),color='r',label='R = {}'.format(R))
            plt.legend(loc='upper right')
            plt.xlabel('Time [$s$]')
            plt.ylabel('Rate [$c/s$]')
            plt.title('{} - light curve segment {}'.format(self.object,n+1))
            plt.show()
        else:
            return t_seg, rate_seg, err_seg, N_gamma, R, T
        
    def channel_to_energy_XTE(self, channel_to_kev):
        """
        Channel to corresponding energy for XTE-data.
        
        **Parameters**:
        
        `channel_to_kev`: np.ndarray    
            Conversion from channel (index) to energy (keV).
        """
        try: 
            minchan = self.MINCHAN
            maxchan = self.MAXCHAN
        except AttributeError: 
            print('Used the Emin/Emax attributes as channel limits.')
            minchan = int(self.Emin)
            maxchan = int(self.Emax)

        if minchan != 0:
            minchan -= 1
        minene = channel_to_kev[minchan]
        if minchan == 0:
            minene = 0
        maxene = channel_to_kev[maxchan]
        
        print('Energy channels (keV): ',minchan+1,'-',maxchan,'(',minene,'-',maxene,')')

        self.Emin = minene
        self.Emax = maxene
        setattr(self, 'deltaE', maxene-minene)
        setattr(self, 'E_mean', (maxene+minene)/2)
        # setattr(self, 'E_err', (maxene-minene)/2) = deltaE/2
    
    def find_energyband(self,f):
        """
        Find Emin and Emax from filename f. For more info, see class doc.
        """

        # -------------------- 1) Check filename f ------------------------- 
        # Please write filename in format: '[low_energy]to[high_energy]kev.lc'
        # e.g. f = '1.5to2kev.lc'
        try: 
            # Start by splitting at "to" or "-"
            if 'to' in f:
                split_at_to = f.split('/')[-1].split('to')
            elif '-' in f:
                split_at_to = f.split('/')[-1].split('-')

            # Find all Emin candidates
            Emin_candidates = re.split('\D+[.]*\D+',split_at_to[0][-6:]) #\D = non-digits, so this splits at letters 

            for minchan in Emin_candidates:
                if '_' in minchan:
                    minchan = minchan.split('_')[-1]
                if bool(re.search(r'\d', minchan)): #need to be a digit
                    Emin = float(minchan)

            # Try to split at 'ke' or other non-digits
            if 'ke' in f:
                Emax_candidates = split_at_to[1].split('k')[0]
            else:
                Emax_candidates = re.split('\D+[.]?\D+',split_at_to[1])

            # Go over all Emax candidates
            if isinstance(Emax_candidates,list):
                for maxcan in Emax_candidates:
                    if bool(re.search(r'\d', maxcan)):
                        Emax = float(maxcan)
            else: 
                Emax = float(Emax_candidates)

        # -------------------- 2) Insert manually ------------------------- 
        except ValueError:
            print('Had a hard time finding Emin/Emax from filename: ',f)
            print('Insert manually:')
            Emin = input('Emin [keV] = ')
            Emax = input('Emax [keV] = ')

        return Emin, Emax
        
    def Fvar_from_lc(self,m,percent_limit=0):
        """
        Calculate the fractional root mean square (rms) variability amplitude from 
        the variance (S^2) in the light curve, by averaging over K=int(np.floor(lc.N/m)) segments. 
        """
        
        Fvar = Fvar_from_lc(self,m,percent_limit)
        if math.isnan(Fvar): #might help to reduce the number of segments 
            Fvar = Fvar_from_lc(self,self.N//10,percent_limit)
        if math.isnan(Fvar): #if still nan, take full light curve 
            Fvar = Fvar_from_lc(self,self.N,percent_limit)
        if math.isnan(Fvar): #if still nan, print it
            print('NOTE: The rms cannot be computed since signal-to-noise ratio is too low; the expectation value of the Poisson variance term is larger than the \measured average variance term, producing negative average excess variances.\n')
        
        return Fvar 

