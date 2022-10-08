#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from ImportPackages import *
from AuxiliaryFunctions import *
from Lightcurve import *

class PowerSpectrum():
    """
    Find the power spectrum by splitting the light curve into K segments. A power spectrum for 
    each segment is computed and the final power spectrum is the average over all these and the error
    is the standard error (sigma/K) over all these.

    **Parameters**:
    
    `lc`: class: 'Lightcurve'-object   
        The light curve data to be Fourier-transformed.

    `m`: int   
        Number of time bins per segment. 

    `normalization`: {'rms' (Miyamoto), 'abs', 'Leahy', or 'none'}, optional, default: 'rms'.     
        What normalization to use, see Vaughan(2003, MNRAS 345).

    `noise`: {'Poisson','Gaussian'}, optional, default: 'Poisson'.     
        For a light curve with 'Poisson'/'Gaussian' errors. 
        
    `B_noise`: float, optional, default: 0     
        Background noise level to be used in the formula for the noise, see Eq. (A2) of Vaughan (2003, MNRAS 345).
        B_noise < 1 means B=B_noise*R.
        B_noise = 2 means that B=np.sqrt(R).
        Else, B=B_noise.

    `percent_limit`: float, optional, default: 90     
        Lower limit for percent of time bins having to be filled, i.e. without gaps, for that segment to be used.

    `timer_on`: boolean, optional, default: True     
        If True, print out progress during computation.
        
    `return_noise_wo_sub`: boolean, optional, default: False     
        If True, noise is returned but not subtracted from powspec.
        If False, returned and subtracted.
        
    `save_all`: boolean, optional, default: False     
        If True, save rate_seg, err_seg, R_seg. These are needed for covariance-computation.
        
    `use_max_data`: boolean, optional, default: True
        If True: if a segment is disregarded due to a larger gap, the next segment will start directly after the gap. 
        If False: if a segment is disregarded, the placement of the next segment is not affected, which means that the data between the gap in the previous segment and the next segment is lost. 
    
    **Attributes**:
    
    `m`: int     
        Number of time bins per segment. 
    
    `channels`: string     
        Channels used during the observation.
    
    `Emin, Emax`: floats     
        Min and max energies of used energy band.
    
    `xf`: np.ndarray     
        The Fourier frequencies. 
    
    `df`: float     
        Frequency resolution.
    
    `fft_rate`: np.ndarray     
        Power spectra. Fast fourier transformation of, and average over all, rate_seg.
    
    `fft_rate_v`: list of np.ndarrays     
        Power spectra for each segment.
        
    `averagePnoise`: float     
        Average noise power for whole light curve.
    
    `Pnoise_v`: list of floats     
        Noise power for each segment.
    
    `Fvar`: float     
        Fractional variance amplitude (=rms) computed from the power spectra by integration.

    `middle_of_log_bins`: np.ndarray     
        The logarithmic midpoint of each log-bin, computed as: 10**(1/2*(np.log10(log_bins[i])+np.log10(log_bins[i+1])))
    
    `fPf`: np.ndarray     
        Power spectra multiplied with f and re-binned logarithmically.
    
    `fPferror`: np.ndarray     
        The standard deviation of all points within a log-bin for the fPf-vector.
    
    `Pf`: np.ndarray     
        Power spectra re-binned logarithmically.
    
    `Pferror`: np.ndarray     
        The standard deviation of all points within a log-bin for the Pf-vector.      
        
    ***In addition, if save_all == True***:
    
    `rate_seg_v`: np.ndarray     
        All segments' rate-vectors.

    `err_seg_v`: np.ndarray     
        All segments' rate-error-vectors.

    `R_seg_v`: np.ndarray     
        All segments' mean count rates.
    """
    
    def __init__(self, lc, m=2**13, normalization='rms', noise='Poisson', B_noise=0,                 percent_limit=90, timer_on=True, return_noise_wo_sub=False, save_all=False,                 use_max_data=True):
        self.m = m
        self.Emin, self.Emax = lc.Emin, lc.Emax
        self.xf, self.fft_rate, self.fft_rate_err, self.fft_rate_v, self.Pnoise_v = self.powspec(lc, normalization, noise, B_noise,percent_limit, timer_on, return_noise_wo_sub, save_all, use_max_data)
        self.df = self.xf[1]-self.xf[0]
        self.averagePnoise = np.mean(self.Pnoise_v)
        self.Fvar = self.Fvar_from_ps()
        self.middle_of_log_bins, self.fPf, self.fPferror, self.Pf, self.Pferror = self.rebin(init=True)
    
    def powspec(self,lc,normalization,noise,B_noise,percent_limit,timer_on,return_noise_wo_sub,save_all,use_max_data):
        """
        Compute the power spectrum. For more info, see class doc.
        """
        
        print('Computing the power spectra using {} bins per segment, normalization "{}", and noise dist "{}"...'.format(self.m,normalization,noise))
        
        # Useful parameters
        dt = lc.dt
        K = int(np.floor(lc.N/self.m)) #num_of_line_seg

        # Average over time segments 
        fft_rate_v = []
        P_noise_v = []
        rate_seg_v = []
        err_seg_v = []
        R_seg_v = []
        num_discarded_segments = 0
        start = timeit.default_timer()
        percent_wo_gaps_v = []
        
        i = 0
        num_of_iterations = K
        bins_back = 0 
        
        while i < num_of_iterations:

            if timer_on:
                timer(i,num_of_iterations-1,start,clear=True)

            # Extract segment 
            t_seg, rate_seg, err_seg, N_gamma, R_seg, T_seg = lc.extract_seg(self.m,n=i,bins_back=bins_back,to_print=False,to_plot=False)

            if save_all:
                rate_seg_v.append(rate_seg)
                err_seg_v.append(err_seg)
                R_seg_v.append(R_seg)

            # Check gaps:
            ## Old version: if abs(dt*m - T_seg) < 1e1: 
            percent_wo_gaps, hist = percent_of_filled_time_bins(t_seg,dt,to_return=True)
            
            if timer_on:
                print('percent of filled time bins (segment {}): {:.2f}'.format(i,percent_wo_gaps))

            if percent_wo_gaps > percent_limit: # if True, then there is no large gap in the segment
                
                # Perform FFT to find Power spectra for one seg
                xf, fft_rate_normalized, P_noise = self.fft_seg(dt,t_seg,rate_seg,err_seg,N_gamma,R_seg,normalization=normalization,noise=noise,B_noise=B_noise,return_noise_wo_sub=return_noise_wo_sub)
                fft_rate_v.append(fft_rate_normalized)
                P_noise_v.append(P_noise)
            
            else:
                percent_wo_gaps_v.append([i,percent_wo_gaps])
                num_discarded_segments += 1

                if use_max_data: # then start new segment directly after the end of a gap
                    
                    # Find last gap in segment
                    for j in range(1,len(hist)):
                        if hist[j-1] == 1 and hist[j] == 0:
                            bin_nr_after_which_no_gaps_exist = j
                    
                    # Take into account all previous gaps so that indexation becomes right 
                    # only relevant when we have more than 1 gap
                    zeros_until_last_gap = np.count_nonzero(hist[:bin_nr_after_which_no_gaps_exist]==0)
                    bin_nr_after_which_no_gaps_exist -= zeros_until_last_gap 
                    
                    # Compute how many steps back we should take in the next segment
                    bins_back += len(t_seg) - bin_nr_after_which_no_gaps_exist
                
                    # Might need to update the total number of segments we have
                    num_of_iterations = int(K + np.floor(bins_back/self.m))
            
            i += 1 

        print('{} of {} segments were disregarded due to lower percent limit set to {:.2f}%:'.format(num_discarded_segments,num_of_iterations,percent_limit))
        if use_max_data:
            print('{} additional segments were used thanks to use_max_data == True'.format(int(np.floor(bins_back/self.m))))
        print('\n'.join('Seg nr = {}, percent of filled time bins = {:.2f}'.format(k[1][0],k[1][1]) for k in enumerate(percent_wo_gaps_v)))
        print('Power spectra done! \n')

        if save_all:
            setattr(self, 'rate_seg_v', rate_seg_v)
            setattr(self, 'err_seg_v', err_seg_v)
            setattr(self, 'R_seg_v', R_seg_v)
        
        fft_rate_mean = np.mean(fft_rate_v,axis=0)
        fft_rate_err = np.std(fft_rate_v,axis=0)/np.sqrt(num_of_iterations-num_discarded_segments) # Standard Error

        return xf, fft_rate_mean, fft_rate_err, fft_rate_v, P_noise_v

    def fft_seg(self,dt,t_seg,rate_seg,err_seg,N_gamma,R_seg,normalization='rms',noise='Poisson',B_noise=0,t_d=0,return_noise_wo_sub=False,to_plot=True):
        """
        Perform discrete FFT (fast-Fourier transform) on a light curve (lc) segment.

        **Parameters**: 
        
        `dt`: float     
            Time resolution.
        
        `t_seg`: np.ndarray     
            The segment's time-vector in seconds.

        `rate_seg`: np.ndarray     
            The segment's count rate(=flux)-vector in counts/seconds.

        `N_gamma`: int     
            Number of counted photons in the segment.

        `R_seg`: float     
            Mean count rate in the segment.

        `normalization`: {'rms' (Miyamoto), 'abs', 'Leahy', or 'none'}, optional, default: 'rms'.     
            What normalization to use, see Vaughan(2003, MNRAS 345).

        `noise`: {'Poisson','Gaussian'}, optional, default: 'Poisson'.     
            For a light curve with "noise" errors.
            
        `B_noise`: float, optional, default: 0     
            Background noise level to be used in the formula for the noise, see Eq. (A2) of Vaughan (2003, MNRAS 345).
            B_noise < 1 means B=B_noise*R.
            B_noise = 2 means that B=np.sqrt(R).
            else, B=B_noise.

        `t_d`: float (not implemented yet)     
            Dead-time of the instrument. 

        **Returns**:
        
        `xf`: np.array     
            The Fourier frequencies. 

        `fft_rate_noise_subtracted`: np.array     
            Fast fourier transformation of rate_seg. Noise subtraction has been made.
        """
        
        # Perform FFT
        xf = np.array(fftfreq(self.m, dt)[1:self.m//2]) #only want positive freq
        fft_rate = fft(rate_seg)[1:self.m//2]

        # Normalize
        fft_rate = self.normalize_power(fft_rate,R_seg,dt,normalization)

        # Remove Noise
        P_noise = self.noise_power(dt,err_seg,R_seg,noise,normalization,B_noise)
        if not return_noise_wo_sub:
            fft_rate = np.array([x-P_noise for x in fft_rate])

        return xf, fft_rate, P_noise

    def noise_power(self,dt,err_seg,R,noise='Poisson',normalization='rms',B_noise=0,t_d=0):
        """
        Calculate and return the Poisson noise power. See Appendix A of Vaughan (2003, MNRAS 345). Do not take dead time into account.

        **Parameters**:
        
        `R`: float     
            Mean count rate in the segment.

        `noise`: {'Poisson','Gaussian'}, optional, default: 'Poisson'.     
            For a light curve with "noise" errors.

        `normalization`: {'rms' (Miyamoto), 'abs', 'Leahy', or 'none'}, optional, default: 'rms'.     
            What normalization to use, see Vaughan(2003, MNRAS 345).

        `B_noise`: float, optional, default: 0     
            Background noise level to be used in the formula for the noise, see Eq. (A2) of Vaughan (2003, MNRAS 345).
            B_noise < 1 means B=B_noise*R.
            B_noise = 2 means that B=np.sqrt(R).
            else, B=B_noise.

        `t_d`: float (not implemented)     
            Dead time of instrument in seconds. 

        **Returns**:
        
        `P_noise`: float     
            Poisson noise power.
        """

        dT_samp_over_dT_bin = 1
        
        if noise == 'Poisson': 
            if isinstance(B_noise, int) or isinstance(B_noise, float):
                if B_noise < 1:
                    B=B_noise*R #B_noise=0.1 works for "first" CYG X-1 data
                elif B_noise == 2:
                    B=np.sqrt(R)
                else:
                    B=B_noise

            if normalization == 'rms':
                P_noise = 2*(R+B)/R**2*dT_samp_over_dT_bin
            elif normalization == 'Leahy': 
                P_noise = 2*(R+B)/R*dT_samp_over_dT_bin #should simply be =2 in most cases
            elif normalization == 'abs': 
                P_noise = 2*(R+B)/R**2*dT_samp_over_dT_bin
            else: #Poisson noise but w/o normalization
                P_noise = (R+B)*self.m/dt* dT_samp_over_dT_bin
        
        elif noise == 'Gaussian':
            if normalization == 'rms':
                P_noise = dT_samp_over_dT_bin*np.mean(err_seg**2)/(R**2*1/(2*dt))
            elif normalization == 'Leahy' or normalization == 'abs':
                print('Sorry, cannot perform that normalization! Not implemented...!')
            else: #Gaussian noise but w/o normalization
                P_noise = np.mean(np.abs(err_seg)**2)*self.m
        
        else: #noise == 'None':
            P_noise = 0
            
        return P_noise

    def normalize_power(self,fft_rate,R,dt,normalization='rms'):
        """
        Normalize the power spectra. 

        **Parameters**:
        
        `fft_rate`: np.array     
            Fast fourier transformation of light curve. To be normalized.

        `R`: float     
            Mean count rate in the segment.

        `dt`: float     
            Time resoltuion in observation.

        `normalization`: {'rms' (Miyamoto), 'abs', 'Leahy', or 'none'}, optional, default: 'rms'.     
            What normalization to use, see Vaughan(2003, MNRAS 345).

        **Returns**:
        
        `normalized_power`: np.ndarray     
            The normalized power. 
        """
        
        if normalization == 'rms':
            return np.array([2*dt/(R**2*self.m)*np.abs(x)**2 for x in fft_rate])
        elif normalization == 'Leahy':
            return np.array([2*dt/(R*self.m)*np.abs(x)**2 for x in fft_rate])    
        elif normalization == 'abs':      
            return np.array([2*dt/self.m*np.abs(x)**2 for x in fft_rate])    
        else: # no-normalization
            return np.array(fft_rate)
    
    def Fvar_from_ps(self):
        """
        Calculate the fractional root mean square (rms) variability amplitude by integrating 
        the power spectra over the full frequency interval. **If** wants the Fvar in a **smaller freq-band**, 
        first use remove_freq() to set power to zero outside given freq band and then use Fvar_from_ps.

        **Returns**:
        
        `F_var`: float     
            The fractional root mean square (rms) variability amplitude.
        """

        return Fvar_from_ps(self.xf,self.fft_rate)

    def rebin(self, num=50, init=False):
        """
        Perform logarithmic rebinning.
        
        **Parameters**:
        
        `num`: int     
            The number of logarithmic frequency bins to create.
        
        `init`: boolean, optional, default: False     
            If False, creates the "middle_of_log_bins, fPf, fPferror, Pf, Pferror" for the first time with num=50 bins.
            If True, updates these attributes, using num bins.
            
        **Returns**:
        
        `middle_of_log_bins, fPf, fPferror, Pf, Pferror`: np.ndarrays     
            See class documentation. 
        """
        
        middle_of_log_bins, fPf, fPferror = log_rebin(self.xf,self.fft_rate*self.xf,self.fft_rate_err*self.xf,num=num)
        middle_of_log_bins, Pf, Pferror = log_rebin(self.xf,self.fft_rate,self.fft_rate_err,num=num)
        
        if init:
            return middle_of_log_bins, fPf, fPferror, Pf, Pferror
        else:
            setattr(self, 'middle_of_log_bins', middle_of_log_bins)
            setattr(self, 'fPf', fPf)
            setattr(self, 'fPferror', fPferror)
            setattr(self, 'Pf', Pf)
            setattr(self, 'Pferror', Pferror)
    
    def plot(self,to_plot='fPf',w=10,with_noise=False,show=True,first=True,label='',color=None):
        """
        Plot power spectrum times freq vs freq. (fPf) or power spectrum vs freq. (Pf).
        
        **Parameters**:
        
        `to_plot`: {'fPf','Pf'}     
            What to plot.
            
        `w`: int, optional, default: 10.
            Width of figure.
            
        `with_noise`: boolean, optional, default: False
            Also display the noise power.
        
        `show`: True   
            If True, show the plot. Needs to be False to be able to add more spectra to the same figure.
        
        `first`: True
            If True, create the figure. Needs to be False when adding more spectra to the same figure.
            
        `label`: str, optional, default: ''
            Label to plot. If default, the energy range will be displayed. 
            
        `color`: optional, default: None
            Color of graph.
            
        """
        
        if first:
            ax = standard_plot(w=w)
        else:
            ax = plt.gca()
        
        if color == None:
            color = next(ax._get_lines.prop_cycler)['color']
        
        if label == '':
            label = '{}-{} keV'.format(self.Emin,self.Emax)
        else:
            print('The energy band is: {}-{} keV'.format(self.Emin,self.Emax))
        
        if to_plot == 'fPf':
            plt.step(self.middle_of_log_bins,self.fPf,where='mid',color=color,label=label)
            ax.errorbar(self.middle_of_log_bins,self.fPf,self.fPferror,fmt=',',color=color,elinewidth=1)
            plt.ylabel('Freq x Power') # P_f in [(RMS/Average)$^2$/Hz]
            if with_noise:
                plt.loglog(self.xf,self.xf*self.averagePnoise,label='Average Noise Power',color='k')
        elif to_plot == 'Pf':
            plt.step(self.middle_of_log_bins,self.Pf,where='mid',color=color,label=label)
            ax.errorbar(self.middle_of_log_bins,self.Pf,self.Pferror,fmt=',',color=color,elinewidth=1)
            plt.ylabel('Power') 
            if with_noise:
                ax.axhline(self.averagePnoise,label='Average Noise Power',color='k')        
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.xlabel('Frequency (Hz)')
        plt.legend()
        if show:
            plt.show()

