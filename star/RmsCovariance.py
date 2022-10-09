#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Rms and Covariance 

from .ImportPackages import *
from .AuxiliaryFunctions import *
from .Lightcurve import *
from .Powerspectrum import *
from .Coherence import *

def rms_vs_energy(lcs,ps_v,m=None,freq_low=None,freq_high=None,units='rms'):
    """
    Compute the fractional variance amplitude (rms) for multiple light curves in a given frequency band.
    
    **Parameters**:

    `lcs`: class: list of 'Lightcurve'-objects
        The light curves to be used.
    
    `ps_v`: class: list of 'PowerSpectrum'-objects
        The corresponding (important!) power spectras to be used.
    
    `m`: int
        Number of bins per segment.
        
    `freq_low/freq_high`: floats
        Lower and upper frequency limits.
        
    `units`: {'abs','frac'}, optional, default: 'abs'
        What unit to return the rms in.        

           
    **Returns**:

    `energy_mid_v`: np.ndarray 
        The average energy from each light curve's energy band. 
    
    `energy_err_v`: np.ndarray 
        Half the size of each light curve's energy band. 
        
    `rms_v`: np.ndarray 
        Absolute/fractional rms.
        
    `rms_err_v`: np.ndarray 
        The error in (absolute/fractional) rms.
    """
    
    print('---------------------------------------------------------------------------------------------------')
    print('                           Computing the rms...')
    print('---------------------------------------------------------------------------------------------------\n')
    start = timeit.default_timer()
    
    rms_v, rms_err_v = [], []
    
    for lc,ps in zip(lcs,ps_v):     
        # Find rms                
        rms, rms_err = rms_freqband(ps,lc,m,freq_low,freq_high,units=units)
        
        rms_v.append(rms)
        rms_err_v.append(rms_err)
        
    time_taken = timeit.default_timer()-start
    print('---------------------------------------------------------------------------------------------------')
    print('                           Rms found (in {:.2f} sec).'.format(time_taken))
    print('---------------------------------------------------------------------------------------------------')
     
    return np.array(rms_v),np.array(rms_err_v)

def rms_freqband(ps,lc=None,m=None,freq_low=None,freq_high=None,units='abs'):
    """
    Compute the rms in a given frequency band for a given energy (as given by 
    the lightcurve and corresponding power spectrum).
    
    **Parameters**:
    
    `ps`: class: 'PowerSpectrum'-object
        The corresponding (important!) power spectra to be used.
        
    `lc`: class 'Lightcurve'-object
        The light curve, whose rms is to be found.
    
    `m`: int
        Number of bins per segment.
        
    `freq_low/freq_high`: floats
        Lower and upper frequency limits.
        
    `units`: {'abs','frac'}, optional, default: 'abs'
        What unit to return the rms in.
        
    **Returns**:
    
    `sigma`: np.float
        Absolute/fractional rms.
    
    `sigma_err`: np.float
        The error in (absolute/fractional) rms. See Eq. (14) of Uttley (2014): "X-ray reverberation around accreting black holes"
    """
    
    # Compute rms via light curve first if no freq-boundaries 
    add_to_string = '\n'
    if freq_low == None and freq_high == None:
        if lc != None:
            rms_lc = Fvar_from_lc(lc,m,percent_limit=80) # --> yield the same value
            add_to_string = ', rms_lc = {:.3f}\n'.format(rms_lc)
        
    # Compute rms via power spectrum
    
    # Frequency range
    freq_low = ps.xf[0] if freq_low == None else freq_low
    freq_high = ps.xf[-1] if freq_high == None else freq_high
    dnu = freq_high-freq_low
    
    # Noise Power is constant in freq!
    P_noise = ps.averagePnoise
    
    # Error variance
    if units == 'abs':
        sigma_noise2 = P_noise*dnu*lc.R**2
    elif units == 'rms':
        sigma_noise2 = P_noise*dnu #now it is rather rms, not sigma...
    else:
        print("You need to pick 'abs' or 'rms' as unit.")
    
    # Find rms within the given interval
    xf, [fft_rate,fft_rate_err] = remove_freq(ps.xf,[ps.fft_rate,ps.fft_rate_err],limit=freq_low,geq=True,disregard=True)
    xf, [fft_rate,fft_rate_err] = remove_freq(xf,[fft_rate,fft_rate_err],limit=freq_high,leq=True,disregard=True)
    
    try: 
        rms = np.sqrt(dnu * np.mean(fft_rate)) #\approx same as Fvar_from_ps(xf,fft_rate)
        rms_err = ps.df*np.sqrt(np.sum(fft_rate_err**2))/(2*rms) #error propagation
        
        """
        # Version 1: Monte Carlo Estimation 
        # Given one power spectrum (for a given energy) we can use fft_rate (the averages) 
        # and fft_rate_err (the standard errors) to sample new power spectra by sampling new values for the 
        # power spectrum at each frequency point i according to the normal distribution N(fft_rate[i],fft_rate_err[i]).
        num_new_ps = 1000
        # Each frequency bin will yield a num_new_ps-long array that will be appended
        ps_new_v = []
        for i in range(0,len(xf)):
            ps_new_v.append(np.random.normal(fft_rate[i], fft_rate_err[i], num_new_ps))
        # Transpose to make each array become a power spectrum on its own
        ps_new_v = np.transpose(np.array(ps_new_v))
        rms_v = [Fvar_from_ps(xf,f) for f in ps_new_v]
        print('rms_err_ver1 = ',np.std(rms_v))
        rms_err = np.std(rms_v)

        # Version 2: 
        sigma_err_ver1 = np.sqrt((np.sqrt(1/(2*lc.N))*np.mean(lc.err**2)/(lc.R**2*rms))**2+\
                            (np.sqrt(np.mean(lc.err**2)/lc.N)*1/lc.R)**2) # Eq. B2 Vaughan2003 "On characterizing..." 

        print('rms_err_ver2 = ',sigma_err_ver1)
        
        # Version 3:
        # In case of unity coherence, the following formula should also work:
        sigma = rms
        sigma_err_ver2 = np.sqrt((2*sigma**2*sigma_noise2+sigma_noise2**2)/(len(lc.t)*sigma**2))
        print('rms_err_ver3 = ',sigma_err_ver2)
        """
 
        if units == 'abs': 
            rms = lc.R*rms
            rms_err = lc.R*rms_err
            """
            sigma_err_ver1 = np.sqrt((np.sqrt(1/(2*lc.N))*np.mean(lc.err**2)/(lc.R**2*rms))**2+\
                                     (np.sqrt(np.mean(lc.err**2)/lc.N)*1/lc.R)**2)    
            sigma_err_ver2 = np.sqrt((2*sigma**2*sigma_noise2+sigma_noise2**2)/(len(lc.t)*sigma**2))
            
            rms = sigma
            rms_err = sigma_err
            """
            print('Eband = {:.2f}-{:.2f} keV, Freq range = {:.3f}-{:.3f} Hz, rms_ps = {:.3f} pm {:.3f} ({:.1f}% error), R = {:.3f}'.format(lc.Emin,lc.Emax,freq_low,freq_high,rms,rms_err,rms_err/rms*100,lc.R)+add_to_string)
        else:
            print('Eband = {:.2f}-{:.2f} keV, Freq range = {:.3f}-{:.3f} Hz, rms_ps = {:.3f} pm {:.3f} ({:.1f}% error)'.format(lc.Emin,lc.Emax,freq_low,freq_high,rms,rms_err,rms_err/rms*100)+add_to_string)
    
    except FloatingPointError:
        rms, rms_err = float('NaN'), float('NaN')
    
    return rms, rms_err


def multiply_w_spectra(lc_v,rms_v,rms_err_v,channel_to_kev,spectral_data):
    """
    Multiply rms (or covariance) with spectra to obtain rms/covariance spectra.
    Might only work well for RXTE (where channels are given by their energy max): https://heasarc.gsfc.nasa.gov/docs/xte/e-c_table.html 
    
    **Parameters**:

    `lc_v`: class: list of 'Lightcurve'-objects    
        The light curves to be used.
            
    `rms_v`: np.ndarray    
        Absolute/fractional rms.
        
    `rms_err_v`: np.ndarray    
        The error in (absolute/fractional) rms.
        
    `channel_to_kev`: np.ndarray    
        Conversion from channel (index) to energy (keV).
    
    `spectral_data`: dict or None, optional, default: None    
        Spectral data of the observation. If dict, rms is scaled with the spectral_data to yield the rms spectra.
        
    **Returns**:

    `rms_v`: np.ndarray    
        Absolute/fractional rms mulitplied with the spectra.
        
    `rms_err_v`: np.ndarray    
        The error in (absolute/fractional) rms mulitplied with the spectra.
    """
    
    print('Multiplying with spectra. Note: only works for XTE-data atm.')
    
    for i in range(0,len(rms_v)):
        lc = lc_v[i]

        # Covert channel min/max to energy min/max
        minchan = lc.MINCHAN
        if minchan != 0:
            minchan -= 1 # For RXTE: since each channel only corresponds to its energy max, we need to sub to get its min (the former channel's max)
        minene = channel_to_kev[minchan]
        if minchan == 0:
            minene = 0
        maxchan = lc.MAXCHAN
        maxene = channel_to_kev[maxchan]

        scale_rms = np.mean(spectral_data['COUNTS/SEC'][minchan:maxchan])/lc.deltaE
        scale_err = np.mean(spectral_data['STATERR/SEC'][minchan:maxchan])/lc.deltaE
        #print('dkeV, counts, scale_rms, scale_err = ',lc.deltaE,', ',np.mean(data['COUNTS/SEC'][minchan:maxchan]),', ',scale_rms,', ',scale_err,'\n')

        rms_v[i] *= scale_rms
        rms_err_v[i] *= scale_err
        
    return rms_v, rms_err_v

def covariance(lc_interest_v,lc_ref,m,alt=1,freq_low=None,freq_high=None,noise='Poisson',percent_limit=90,units='abs',to_plot=False):
    """
    Compute the covariance, which is more robust than the rms. 
    
    **Parameters**:

    `lc_interest_v`: list of 'Lightcurve'-objects    
        The light curves of interest to be used in the covariance computation. 
         
    `lc_ref`: 'Lightcurve'-object       
        Reference light curve.
         
    `m`: int    
        Number of time bins per segment. 
        
    `alt`: int {1, 2}, optional, default: 1    
        What method to compute the covariance with.
        If full freq range: 1 = for whole light curve directly, 2 = segment-wise. 
        If smaller freq range: 1 = using FFT and inverseFFT, 2 = using coherence.
    
    `freq_low/freq_high`: floats, optional, default: None    
        Lower and upper frequency limits.
        
    `noise`: {'Poisson','Gaussian'}, optional, default: 'Poisson'.    
        For a light curve with 'Poisson'/'Gaussian' errors. 
     
    `percent_limit`: float, optional, default: 90    
        Lower limit for percent of time bins having to be filled, i.e. without gaps, for that segment to be used.
    
    `units`: {'abs','frac'}, optional, default: 'abs'    
        What unit to return the rms in.
    
    `to_plot`: boolean (default: False)    
        If True, a figure for different ways to compute the covariance is displayed. 
        
    **Returns**:

    `cov_v`: array    
        The normalised coviarances as calculated by Eq(2) and Eq(3) from Wilkinson(2009).
        
    `cov_err_v`: array    
        The statistical errors of the covariance.  
    """
    
    cov_v, cov_err_v = [], []

    for i in range(0,len(lc_interest_v)):
        lc_v = [lc_interest_v[i],lc_ref]
    
        print('---------------------------------------------------------------------------------------------------')
        print('                           Covariance computation about to begin...')
        print('---------------------------------------------------------------------------------------------------\n')

        start = timeit.default_timer()

        if np.all(lc_v[0].rate == lc_v[1].rate):
            print('The light curves are the same! Try again.\n')
            continue
            
        print('Light curve 1: {}-{} keV'.format(lc_v[0].Emin,lc_v[0].Emax))
        print('Light curve 2: {}-{} keV (this is the reference band)\n'.format(lc_v[1].Emin,lc_v[1].Emax))

        # If overlapping energy bands, this needs to be handled
        lc_v = subtract_overlapping_energybands(lc_v)

        # Extract data
        lc_X = lc_v[0]
        t_X, rate_X, err_X, R_X = lc_X.t, lc_X.rate, lc_X.err, lc_X.R
        # and from Reference Band:
        lc_Y = lc_v[1]
        t_Y, rate_Y, err_Y, R_Y = lc_Y.t, lc_Y.rate, lc_Y.err, lc_Y.R

        errX = np.mean(err_X**2)
        errY = np.mean(err_Y**2)

        comparison = t_X == t_Y
        assert comparison.all(), 'Time arrays are not identical. The two light curves must come from the same observation...'

        # Useful parameters (take from reference band, but is the same as the corresponding X-band values)
        dt = t_Y[1]-t_Y[0]
        N = len(t_Y)

        # Full frequency range
        if freq_low == None and freq_high == None: 

            if alt == 1:
                ## ------------------------------------ Cov for whole light curve ---------------------------------------------
                print('---------------------------------------------------------------------------------------------------')
                print('                           Computing covariance for whole light curve directly...')
                print('---------------------------------------------------------------------------------------------------\n')

                # Prep Work to normalize covariance and to compute covariance error
                errX = np.mean(err_X**2)
                errY = np.mean(err_Y**2)

                var_ex_X_lc = 1/(len(rate_X)-1)*np.sum((rate_X-np.mean(rate_X))**2)-errX
                var_ex_Y_lc = 1/(len(rate_Y)-1)*np.sum((rate_Y-np.mean(rate_Y))**2)-errY
                print('var_ex_X_lc = {}, var_ex_Y_lc (ref) = {}'.format(var_ex_X_lc,var_ex_Y_lc))
                if var_ex_X_lc > var_ex_Y_lc:
                    print("Note that the reference band should be the band with highest absolute variability; this is not the case here.")
                    print('var_ex_X_lc = {}, var_ex_Y_lc (ref) = {}'.format(var_ex_X_lc,var_ex_Y_lc))
                Fvar_lc = np.sqrt(var_ex_Y_lc)/np.mean(rate_Y) 

                # Compute Covariance
                cov = 1/(len(rate_Y)-1)*np.sum((rate_X-np.mean(rate_X))*(rate_Y-np.mean(rate_Y)))

                # Normalize
                cov = cov/np.sqrt(var_ex_Y_lc)
                # Error
                cov_err = np.sqrt((var_ex_X_lc*errY+var_ex_Y_lc*errX+errX*errY)/(len(t_X)*var_ex_Y_lc))

                if units != 'abs': 
                    cov /= R_X 
                    cov_err /= R_X

                print('Cov = ',cov,' pm ',cov_err,'\n')

            elif alt == 2:
                ## ------------------------------------- Cov for segments ---------------------------------------------------------
                print('---------------------------------------------------------------------------------------------------')
                print('                           Computing covariance segment-wise...')
                print('---------------------------------------------------------------------------------------------------\n')

                # Prep Work to normalize covariance and to compute covariance error
                ps_X = PowerSpectrum(lc_v[0],m=m,percent_limit=percent_limit,noise=noise,timer_on=False,save_all=True)
                ps_Y = PowerSpectrum(lc_v[1],m=m,percent_limit=percent_limit,noise=noise,timer_on=False,save_all=True)

                fft_rate_meanX = ps_X.fft_rate
                fft_rate_meanY = ps_Y.fft_rate

                var_ex_X = [1/(m-1)*np.sum((rate_seg-R_seg)**2)-np.mean(err_seg**2) for rate_seg, R_seg, err_seg in zip(ps_X.rate_seg, ps_X.R_seg, ps_X.err_seg)]
                var_ex_Y = [1/(m-1)*np.sum((rate_seg-R_seg)**2)-np.mean(err_seg**2) for rate_seg, R_seg, err_seg in zip(ps_Y.rate_seg, ps_Y.R_seg, ps_Y.err_seg)]

                var_err_X = [np.mean(e**2) for e in ps_X.err_seg]
                var_err_Y = [np.mean(e**2) for e in ps_Y.err_seg]

                P_X_mean = np.mean(ps_X.fft_rate)
                P_Y_mean = np.mean(ps_Y.fft_rate)

                # Should we take the average over all segments to get the final variance error? 
                var_err_X_mean = np.mean(var_err_X,axis=0)
                var_err_Y_mean = np.mean(var_err_Y,axis=0)

                # Excess Variance 
                var_ex_X_lc = np.mean(var_ex_X)
                var_ex_Y_lc = np.mean(var_ex_Y)
                if var_ex_X_lc > var_ex_Y_lc:
                    print("Note that the reference band should be the band with highest absolute variability; this is not the case here.")
                    print('var_ex_X_lc = {}, var_ex_Y_lc (ref) = {}'.format(var_ex_X_lc,var_ex_Y_lc))

                # Quantities neeeded:
                cov = [1/(m-1)*np.sum((rate_seg_X-R_X_seg)*(rate_seg_Y-R_Y_seg)) for rate_seg_X, rate_seg_Y, R_X_seg, R_Y_seg in zip(ps_X.rate_seg, ps_Y.rate_seg, ps_X.R_seg, ps_Y.R_seg)]
                if units != 'abs': 
                    cov = [c/R for c,R in zip(cov,ps_X.R_seg)]

                # 1) Small difference between taking average over all covariances and using full excess variance to normalize
                cov_mean = np.mean(cov,axis=0)
                cov_norm_alt1 = cov_mean/np.sqrt(var_ex_Y_lc)
                print('cov_norm = ',cov_norm_alt1)
                # 2) vs using excess variance from each segment to normalize and then taking the average 
                cov_seg_norm = np.array(cov)/np.sqrt(var_ex_Y)
                cov_norm_alt2 = np.mean(cov_seg_norm,axis=0)
                print('cov_norm = ',cov_norm_alt2)

                cov = cov_norm_alt2
                cov_err = np.sqrt((var_ex_X_lc*var_err_Y_mean+var_ex_Y_lc*var_ex_X_lc+var_err_X_mean*var_err_Y_mean)/(N*var_ex_Y_lc))
                if units != 'abs': 
                    cov_err /= R_X
                print('Cov = ',cov,' pm ',cov_err,'\n')

        # If not full freq. range
        else:
            if alt == 1:
                ## ------------------------------------ Cov using FFT and inverseFFT ---------------------------------------------       
                print('---------------------------------------------------------------------------------------------------')
                print('                           Computing covariance using FFT and inverseFFT...')
                print('---------------------------------------------------------------------------------------------------\n')

                print('Perform FFT and IFFT to extract light curves with variability only in the given freq band.\n')

                xf = np.array(fftfreq(N, dt))[1:N//2]
                # Pick out frequency range 
                if freq_low == None:
                    freq_low = xf[0]
                if freq_high == None:
                    freq_high = xf[-1]

                rate_X = pick_out_freq_from_lc(lc_v[0], freq_low, freq_high)
                rate_Y = pick_out_freq_from_lc(lc_v[1], freq_low, freq_high)

                print('Find excess variance through the rms (i.e. from a normalized power spectra) in the given freq range.\n')
                # rms needs to be taken from normalized power spectra in the relevant freq range
                ps_X, ps_Y = PowerSpectrum(lc_v[0],m=N,timer_on=False,noise=noise,save_all=True,percent_limit=0), PowerSpectrum(lc_v[1],m=N,timer_on=False,noise=noise,save_all=True,percent_limit=0)
                fft_rate_meanX, fft_rate_meanY = ps_X.fft_rate, ps_Y.fft_rate

                xf, fft_rates = remove_freq(xf,[fft_rate_meanX,fft_rate_meanY],freq_low,geq=True,disregard=False)
                xf, fft_rates = remove_freq(xf,[fft_rates[0],fft_rates[1]],freq_high,leq=True,disregard=False)

                FvarX, FvarY = Fvar_from_ps(xf, fft_rates[0]), Fvar_from_ps(xf, fft_rates[1])
                #print('FvarY = ',FvarY)
                var_ex_X_lc, var_ex_Y_lc = (FvarX * R_X)**2, (FvarY * R_Y)**2
                if var_ex_X_lc > var_ex_Y_lc:
                    print("Note that the reference band should be the band with highest absolute variability; this is not the case here.")
                    print('var_ex_X_lc = {}, var_ex_Y_lc (ref) = {}'.format(var_ex_X_lc,var_ex_Y_lc))

                # Compute cov
                cov = 1/(len(rate_Y)-1)*np.sum((rate_X-np.mean(rate_X))*(rate_Y-np.mean(rate_Y)))
                # Normalize                                     
                cov = cov/np.sqrt(var_ex_Y_lc)

                # Error
                P_noise_X_full = np.mean(err_X**2)/(R_X**2*1/(2*dt))
                P_noise_Y_full = np.mean(err_Y**2)/(R_Y**2*1/(2*dt))
                errX = P_noise_X_full*R_X**2*(freq_high-freq_low)
                errY = P_noise_Y_full*R_Y**2*(freq_high-freq_low)
                cov_err = np.sqrt((var_ex_X_lc*errY+var_ex_Y_lc*errX+errX*errY)/(len(t_X)*var_ex_Y_lc))

                if units != 'abs': 
                    cov /= R_X 
                    cov_err /= R_X
                print('Cov = ',cov,' pm ',cov_err,'\n')

            elif alt == 2:
                ## ----------------------------------------------- Cov using coherence ------------------------------------------------------------
                print('---------------------------------------------------------------------------------------------------')
                print('                          Covariance using coherence...')
                print('---------------------------------------------------------------------------------------------------\n')

                upper_freq_lim = freq_high
                assert upper_freq_lim != None, "You need to set an upper freq. limit; will be the same for all lc."

                # Compute intrinsic coherence 
                xf_coh, gamma2, delta_gamma2_int = coherence_intrinsic(lc_v,m_init=2**16)
                xf_coh, gamma2 = remove_freq(xf_coh,gamma2,limit=upper_freq_lim,leq=True)
                dnu = upper_freq_lim-xf_coh[0]

                # Should use the fft_rate_meanX/fft_rate_meanY from the coherence computation (that uses different m), 
                # but this is not implemented...
                # ---- This part will do until fixed.  Small difference... -----------
                ps_X = PowerSpectrum(lc_v[0],m=m,percent_limit=percent_limit,noise=noise,timer_on=False,save_all=True)
                ps_Y = PowerSpectrum(lc_v[1],m=m,percent_limit=percent_limit,noise=noise,timer_on=False,save_all=True)
                xf = ps_X.xf

                fft_rate_meanX = ps_X.fft_rate
                fft_rate_meanY = ps_Y.fft_rate

                # ----------------------------------------------------------------------------------------------------

                _, fft_rate_meanX_test = remove_freq(xf,fft_rate_meanX,limit=upper_freq_lim,leq=True)
                P_X_mean = np.mean(fft_rate_meanX_test)

                _, fft_rate_meanY_test = remove_freq(xf,fft_rate_meanY,limit=upper_freq_lim,leq=True)
                P_Y_mean = np.mean(fft_rate_meanY_test)

                # Compute covariance
                cov = R_X*np.sqrt(np.mean(gamma2)*P_X_mean*dnu)
                print('Cov = ',cov)

                # Error
                P_noise_X_full = np.mean(err_X**2)/(R_X**2*1/(2*dt))
                P_noise_Y_full = np.mean(err_Y**2)/(R_Y**2*1/(2*dt))

                errX = R_X*P_noise_X_full*dnu
                errY = R_Y*P_noise_Y_full*dnu
                var_ex_X_lc = R_X*P_X_mean*dnu
                var_ex_Y_lc = R_Y*P_Y_mean*dnu
                if var_ex_X_lc > var_ex_Y_lc:
                    print("Note that the reference band should be the band with highest absolute variability; this is not the case here.")
                    print('var_ex_X_lc = {}, var_ex_Y_lc (ref) = {}'.format(var_ex_X_lc,var_ex_Y_lc))

                cov_err = np.sqrt((cov**2*errY+var_ex_Y_lc*errX+errX*errY)/(len(t_X)*var_ex_Y_lc))

                if units != 'abs': 
                    cov /= R_X 
                    cov_err /= R_X

                print('Cov = ',cov,' pm ',cov_err,'\n')

        time_taken = timeit.default_timer()-start
        print('---------------------------------------------------------------------------------------------------')
        print('                           Covariance found (in {:.2f} sec).'.format(time_taken))
        print('---------------------------------------------------------------------------------------------------') 
    
        cov_v.append(cov)
        cov_err_v.append(cov_err)
    
    return cov_v, cov_err_v

