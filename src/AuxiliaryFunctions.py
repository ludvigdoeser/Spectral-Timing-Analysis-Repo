#!/usr/bin/env python
# coding: utf-8

# # Put log_rebin, standard_plot etc here

# In[ ]:


import sys
sys.path.append('..')
from ImportPackages import *
#from Lightcurve import *
# Auxiliary Functions

def standard_plot(h=4,w=10,fontsize=16):
    """
    Standard plot to enable use of the same figsize, fontsize and font.family in all figures.
    
    **Parameters**:   
    `h,w`: (float, float), optional, default: h=4, w=10     
        Height and width of the figure, i.e., figsize=(w,h).
    
    `fontsize`: int, optional, default: 16
        Fontsize of labels in figure.
        
    **Returns**: 
    
    `ax`: Axes,     
        Axes of the figure. 
    """
    
    fig = plt.figure(figsize=(w,h))
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'font.size': fontsize})
    plt.rcParams['font.family'] = 'Times'
    plt.rc('text', usetex=True) 
    
    return plt.gca()

def timer(i,nr,start,clear=True):
    """
    Print out the loop progress, and try to estimate the time remaining.
    
    **Parameters**:
    
    `i, nr`: ints      
        The current loop number, and the total number of loops.
    
    `start`: float      
        Start of the loop, as given by: timeit.default_timer(). 
    
    `clear`: boolean      
        Clear the output between loops. 
    """
    
    stop=timeit.default_timer()
    if (i/nr*100) < 10:
        expected_time="Calculating..."
    else:
        time_perc=timeit.default_timer()
        expected_time=np.round(((time_perc-start)/(i/nr))/60,2)
    if clear == True:
        clear_output(wait=True)
    
    print("Calculating the power spectra...")
    print("Current progress: {}%".format(np.round(i/nr*100,2)))
    print("Current run time: ",np.round((stop-start)/60,2)," minutes")
    print("Expected run time: ",expected_time," minutes")
    

def extract_fits(filename,keywords=[],p=0):
    """ 
    Extract (and print info about) a fits-file. 
    
    **Parameters**:
    
    `filename`: string     
        The path to the .fits-file.
    
    `keywords`: list of strings     
        Keywords (and corresponding values) apart from the ones mentioned below to return. 
        *All data-keys are always automatically returned (for a lightcurve, these are: {"TIME", "RATE", "ERROR","FRACEXP"}).
        *Note also that the header-keys {"CONTENT","OBJECT"} are returned if they exist. 
    
    `p`: boolean, optional, default: 0     
        If True (=1), print the headers of the fits-file.
        
    **Returns**:
    
    `data`: dictionary     
        Keys and the corresponding data values.
    """
    
    print('Loading fits from filename: ',filename)
    with fits.open(filename) as hdulist: 
    
        # HEADER
        header0 = hdulist[0].header
        header1 = hdulist[1].header
        
        # DATA KEYS
        binaryext = hdulist[1].data
        binarytable = Table(binaryext)
        keys = binarytable.keys()
        
        data = {}
        for key in keys:
            data[key] = hdulist[1].data[key]
        
        # HEADER KEYS
        try:
            data['CONTENT'] = header0['CONTENT']
            content_exist = True
        except KeyError:
            print('There is no content-information for this fits.')
            content_exist = False
        try:
            data['OBJECT'] = header0['OBJECT']
        except KeyError:   
            print('There is no object-information for this fits.')
          
        # Extract extra keyswords
        for extra_key in keywords:  
            try: 
                data[extra_key] = header0[extra_key]
                print("Found key {} in header0".format(extra_key))
                cannot_find = False
            except KeyError:
                print("There is no key {} in header0, let's have a look at header1.".format(extra_key))
            try: 
                data[extra_key] = header1[extra_key]
                print("Found key {} in header1.".format(extra_key))
                cannot_find = False
            except KeyError:
                print("There is no key {} in header1 either, sorry...".format(extra_key))
                cannot_find = True
            
            if cannot_find: 
                matchingK0 = matchingKeys(header0, extra_key)
                matchingK1 = matchingKeys(header1, extra_key)
                print('Matching keys in header0 = ',matchingK0)
                print('Matching keys in header1 = ',matchingK1)
                for K in matchingK0:
                    use = input('Do you want to extract the header0 key {}? [y/n] '.format(K))
                    if use == 'y' or use == 'yes' or use == 'Y':
                        data[K] = header0[K]
                for K in matchingK1:
                    use = input('Do you want to extract the header1 key {}? [y/n] '.format(K))
                    if use == 'y' or use == 'yes' or use == 'Y':
                        data[K] = header1[K]
         
        # Print out all info if true:
        if p:
            print('hdu.info()')
            print(hdulist.info(),'\n')
            print('Header0:')
            print(repr(header0),'\n') #repr() prints the info into neat columns
            print('Header1:')
            print(repr(header1),'\n') #repr() prints the info into neat columns
            print(binarytable[0:10],'\n')
        else:
            if content_exist:
                print('The keys to the {} data are: {}'.format(data['CONTENT'],data.keys()))
            else:
                print('The keys to the data are: {}'.format(data.keys()))
        print('Loading fits done. \n')
        
        return data

def matchingKeys(dictionary, searchString):
    """
    Find if a string is contained within any of the dictionary's keys.
    
    **Parameters**:
    
    `dictionary`: dict
        The dictionary to be searched.
        
    `searchString`: str
        The string to be looked after.
    
    ** Returns**:
    
    `matches`: list
        A list with all key-matches.
    """
    
    matches = [key for key,val in dictionary.items() if searchString in key]    
    
    return matches

def log_rebin(xf,take_average_of,err=[],low_lim=None,high_lim=None,num=50):
    """
    Distribute the data into logarithmic frequency bins and compute the bin-average of the data.
    
    **Parameters**:
    
    `xf`: np.ndarray     
        The frequency-vector.
    
    `take_average_of`: np.ndarray     
        The quantity to take average of. 
        
    `err`: np.ndarray, optional, default: empty list     
        Errors corresponding to the "take_average_of"-quantity.     
        If empty, compute the standard deviation of all values that end up in the same bin.
    
    `low_lim, high_lim`: floats, optional, default: None     
        The interval to bin into log-bins: 10^(low_lim) to 10^(high_lim).
        If None, limits are given by the respective end point of the frequency vector.
    
    `num`: int     
        The number of logarithmic frequency bins to create.
    
    **Returns**:
    
    `middle_of_log_bins`: np.ndarray     
        The logarithmic midpoint of each log-bin, computed as:     
    $$10^{\\frac{1}{2}\\left(\\texttt{np.log10}(\\texttt{log_bins}[i])+\\texttt{np.log10}(\\texttt{log_bins}[i+1])\\right)}$$
        
    `average`: np.ndarray     
        The average of the data within a log-bin.
         
    `error`: np.ndarray     
        The propagated error, or (in case empty error-input) the standard deviation of all points within a log-bin. 
    """
    
    if isinstance(xf,np.ndarray) and isinstance(take_average_of,np.ndarray):
        
        # Find limits and make log bins
        if low_lim == None:
            low_lim = np.log10(np.amin(xf))
        if high_lim == None:
            high_lim = np.log10(np.amax(xf))
        eps = 10e-15 #to include the final point
        log_bins = np.logspace(low_lim, high_lim+eps, num=num)

        # Find the logarithmic mid-points
        middle_of_log_bins = []
        for i in range(0,len(log_bins)-1):
            middle_of_log_bins.append(10**(1/2*(np.log10(log_bins[i])+np.log10(log_bins[i+1]))))

        # Determine what freq-values correspond to what bin
        digitized = np.digitize(xf, log_bins)

        # Sort into bins and make sure no bin is empty 
        # Note that num-1 is equal to len(middle_of_log_bins)
        # i-1 for middle_of_log_bins and i for average due to indexing.. 
        middle_of_log_bins = [middle_of_log_bins[i-1] for i in range(1, num) if len(take_average_of[digitized == i])!=0] 
        average = [take_average_of[digitized == i].mean() for i in range(1, num) if len(take_average_of[digitized == i])!=0]
        
        # Standard deviation for points in this bin if no error as input:
        if len(err) == 0:
            error = [take_average_of[digitized == i].std() for i in range(1, num) if len(take_average_of[digitized == i])!=0]
        # If error as input, compute error propagation (len(err[digitized==i]) has been moved out from root):
        else:
            error = [np.sqrt(np.sum((err[digitized == i])**2))/len(err[digitized == i]) for i in range(1, num) if len(take_average_of[digitized == i])!=0]

        return np.array(middle_of_log_bins), np.array(average), np.array(error)

    else:
        print('Input needs to be np.ndarrays, try again.')


def percent_of_filled_time_bins(t_seg,dt,to_return=True):
    """
    Compute percent of time bins being filled, i.e. without gaps.
    
    **Parameters**:
     
    `t_seg`: np.ndarray     
        Segment's time vector. Should start from zero, i.e. t_seg[0] = 0. 
        
    `dt`: np.float     
        Time resolution of observation.
        
    `to_return`: boolean (default: False)     
        If false, print gap percent, otherwise return it. 
        
    **Returns**:
    
    `percent_wo_gaps`: np.float     
        percent of time bins being filled, i.e. without gaps.
    
    `hist`: np.ndarray    
        A list containing 1 for filled bin, 0 for unfilled bin.
    """
    
    t_seg_temp = np.copy(t_seg)
    
    # Linear transformation of segment's first element to t_seg=0    
    if int(t_seg_temp[0]) != 0:
        t_seg_temp -= t_seg_temp[0]
    
    # Count number of filled (= non-empty = no gap) bins
    num_bins = math.ceil(t_seg_temp[-1]/dt)
    hist, edges = np.histogram(t_seg_temp,bins=num_bins,range=(0, dt*num_bins))
    
    # percent of filled bins (i.e. without gap)
    percent_wo_gaps = np.sum(hist)/num_bins*100

    if to_return:
        return percent_wo_gaps, hist 
    else: 
        print('perc_wo_gaps = {:.4f}'.format(percent_wo_gaps))

def Fvar_from_lc(lc, m, percent_limit):
    """
    Calculate the fractional root mean square (rms) variability amplitude (over all frequencies) from 
    the variance (S^2) in the light curve, by averaging over K=int(np.floor(lc.N/m)) segments. 
    
    **Parameters**:
    
    `lc`: class: 'Lightcurve'-object     
        The light curve data to be Fourier-transformed.
        
    `m`: int     
        Number of time bins per segment. 

    `percent_limit`: float, optional, default: 90     
        Lower limit for percent of time bins having to be filled, i.e. without gaps, for that segment to be used.

    **Returns**:
    
    `F_var`: float     
        The fractional root mean square (rms) variability amplitude.   
        See Eq. (10) from: "On characterizing the variability properties of X-ray light curves from active galaxies" (Vaughan, 2003)
    """
    
    # Find rms and return NaN if cannot be computed
    np.seterr(all='raise')
    try: 
        # Useful parameters
        dt = lc.dt
        K = int(np.floor(lc.N/m)) #num_of_line_seg

        # Average over time segments 
        rms_lc_seg_v = []
        for i in range(0,K):

            # Pick out one line segment (ls)
            t_seg, rate_seg, err_seg, N_gamma, R_seg, T_seg = lc.extract_seg(m,n=i,to_print=False,to_plot=False)

            if percent_limit > 0:
                percent_wo_gaps, _ = percent_of_filled_time_bins(t_seg,dt,to_return=True) 
            else: 
                percent_wo_gaps = 100
            
            if percent_wo_gaps > percent_limit: # if True, then there is no large gap in the segment
                # Perform FFT to find Power spectra for one seg
                S2 = 1/(len(rate_seg)-1)*np.sum((rate_seg-R_seg)**2)
                MSE = np.mean(err_seg**2)
                F_var = np.sqrt((S2-MSE)/R_seg**2)
                rms_lc_seg_v.append(F_var)
            else:
                pass

        return np.mean(rms_lc_seg_v)

    except FloatingPointError:
        return float("NaN")

def Fvar_from_ps(xf,fft_rate):
    """
    Calculate the fractional root mean square (rms) variability amplitude (over all frequencies) by 
    discretely integrating the power spectra over the full frequency interval.     
     
    Note: **for a smaller freq band**, call remove_freq() prior to calling this function.    
    
    **Parameters**:
    
    `xf, fft_rate`: np.ndarrays     
        Frequency vector and power spectrum (the Fourier transformed rate).
        
    **Returns**:
    
    `F_var`: float     
        The fractional root mean square (rms) variability amplitude.
    """
    np.seterr(all='raise')
    try: 
        df = xf[1]-xf[0]
        F_var = np.sqrt(df*np.sum(fft_rate))
        return F_var
    
    except FloatingPointError:
        return float("NaN")
    
def load_lightcurve(data):
    """ 
    Split the data from a light curve into time, flux, and error vectors. 
    
    **Parameters**:
    
    `data`: dictionary     
        Should contain the keys ['TIME', 'RATE', 'ERROR', 'FRACEXP', 'CONTENT'] with corresponding values.    
        The data['CONTENT'] ought to be 'LIGHT CURVE'.   
        
    **Returns**:
    
    `t, rate, error`: np.ndarrays     
        The data values for a light curve: time, rate, rate_error, respectively.
    """

    assert data['CONTENT']=='LIGHT CURVE', 'Data does not come from a light curve object.'

    t = np.array(data['TIME'])
    rate = np.array(data['RATE'])
    err = np.array(data['ERROR'])
    
    return t, rate, err

def remove_freq(xf,other_quantites,limit,leq=False,l=False,geq=False,g=False,disregard=True):
    """
    Remove frequencies (f) and other quantites' (OQ) corresponding values for those f.    
    Example: limit = 10 and leq = True: means that only f <= 10 are saved.
    
    **Parameters**:
    
    `xf`: np.ndarray     
        Frequency vector.
        
    `other_quantites`: list of np.ndarrays
        The quantities having a value for each frequnecy in xf.
    
    `limit`: float     
        The numerical value of the limit. The type of limit is determined next:   
    
    `leq, l, geq, g`: Boolean (default: False)     
        less or equal than (leq), less than (l), greater or equal to (geq), greater than (g) the limit = will be SAVED.
        
    `disregard`: Boolean (default: False)     
        If True, then the frequency and OQ will be cropped and returned as shorter vectors.
        If False, returned in the same length with OQs' elements outside given freq range being set to 0.

    **Returns**:
    
    `freq`: np.ndarray     
        The frequency vector. 
        
    `other_quantites`: list of np.ndarrays or np.ndarray (if just one quantity)     
        After disregarding values (or setting them to zero) corresponding to frequencies outside the desired interval. 
    """
    
    assert type(limit) != list, 'Unfortunately, you cannot fix more than one limit at once.'
    assert int(leq)+int(l)+int(geq)+int(g) <= 1, 'Unfortunately, you cannot fix more than one limit at once.'
    
    try: 
        return_as_nparray = False
        if type(other_quantites) != list: # i.e. we only have one quantity 
            other_quantites = [other_quantites]
            return_as_nparray = True

        if leq:
            for i in range(0,len(other_quantites)):
                if disregard:
                    other_quantites[i] = np.array(other_quantites[i][xf<=limit])
                else:
                    other_quantites[i][xf>limit] = 0
            if disregard:
                xf = xf[xf<=limit]
            else:
                pass

        if l:
            for i in range(0,len(other_quantites)):
                if disregard:
                    other_quantites[i] = np.array(other_quantites[i][xf<limit])
                else:
                    other_quantites[i][xf>=limit] = 0
            if disregard:
                xf = xf[xf<limit]
            else:
                pass

        if geq:
            for i in range(0,len(other_quantites)):
                if disregard: 
                    other_quantites[i] = np.array(other_quantites[i][xf>=limit])
                else:
                    other_quantites[i][xf<limit] = 0
            if disregard:
                xf = xf[xf>=limit]
            else:
                pass

        if g:
            for i in range(0,len(other_quantites)):
                if disregard: 
                    other_quantites[i] = np.array(other_quantites[i][xf>limit])
                else:
                    other_quantites[i][xf<=limit] = 0
            if disregard:
                xf = xf[xf>limit]
            else:
                pass

        if return_as_nparray:
            return np.array(xf), other_quantites[0]
        else:
            return np.array(xf), other_quantites
    
    except IndexError:
        print('Could not change anything, try again with new limits.')
        
        return xf, other_quantites
    
def error_change(err,err_lim=1):
    """
    Set error to zero if exceeds err_lim.
    
    **Parameters**:
    
    `err`: np.ndarray     
        Error vector.
        
    `err_lim`: float, optional, default: 1    
        Error limit. If an error element is larger than err_lim, it is set to zero. 
    
    **Returns**:
    
    `err`: np.ndarray     
        Error vector after setting elements that exceeds err_lim to zero.
    """
    
    e_temp = []
    for e in err:
        if e < err_lim:
            e_temp.append(e)
        else:
            e_temp.append(0)
    err = e_temp
    return err

    
def pick_out_freq_from_lc(lc, freq_low, freq_high, to_plot=False):
    """
    Fourier transforms a light curve, sets the power = 0 for all freq outside freq_range: freq_low-freq_high,
    and then performs an inverse Fourier transform back to time-domain.
    
    **Parameters**:
    
    `lc`: class 'Lightcurve'-object     
            The light curve, whose rms is to be found.
            
    `freq_low/freq_high`: floats, optional, default: None     
        Lower and upper frequency limits.
    
    **Returns**:
    
    `rate_ifft`: np.ndarray     
        The inverse Fourier transformed rate-vector, with mean = 0.
    """
    
    # Pick out relevant quantties from the lightcurves
    t, dt, rate, err, R, N = lc.t, lc.dt, lc.rate, lc.err, lc.R, lc.N
    
    # FFT on full light curve
    xf = np.array(fftfreq(N, dt))[1:N//2]
    fft_rate_unnormalized = fft(rate-R)[1:N//2]
                    
    if freq_low == xf[0] and freq_high == xf[-1]: 
        print("You're using the full freq. range.")
    else:
        print('Inverse FFT using only the freq interval: [{},{}]'.format(freq_low, freq_high)) 
    
    # Set fft_rate = 0 for frequencies outside given range
    xf, fft_rate = remove_freq(xf,fft_rate_unnormalized,freq_low,geq=True,disregard=False)
    xf, fft_rate = remove_freq(xf,fft_rate,freq_high,leq=True,disregard=False)
    
    # Transform back 
    rate_ifft = ifft_smallfreqband(fft_rate)
    print('Inverse FFT performed.\n')

    if to_plot:
        standard_plot()
        plt.plot(t,rate-R,label='Lc with all freq')
        plt.plot(t,rate_ifft,label='Lc using only f = [{},{}]'.format(freq_low,freq_high))
        plt.legend()
        ax = plt.gca()
        ax.set_xlim([t[0],t[500]])
        plt.show()
        
    return rate_ifft
    
def ifft_smallfreqband(fft_rate):
    """
    Perform inverse transformation from freq domain to time domain. Is called upon in pick_out_freq_from_lc().
    
    **Parameters**:
     
    `fft_rate`: np.ndarray     
        The power spectra after having set values outside freq range to 0.
        
    **Returns**:
    
    `rate`: np.ndarray     
        New rate vector, now only containing the frequencies in fft_rate. 
    """
    
    # Add the negative freq again (that were removed due to [1:m//2])
    y_together = np.append(np.zeros(1),fft_rate)
    y_together = np.append(y_together,np.zeros(1))
    y_together = np.append(y_together,np.conjugate(np.flip(fft_rate)))
    # Perform ifft
    yinv = ifft(y_together)
    # Make to np.ndarray and extract only the real values
    yinv = np.array(yinv)
    rate = yinv.real  
    return rate 

def rms_vs_channels(lc_v, rms_v, rms_err_v, channel_to_kev, overlap=0):
    """
    Create a list with rms for each channel. If the energy bands are chosen to 50-100, 100-150 etc, then
    channels 50-99 will be given the rms for the first band, and 100-149 will be given the rms of the second band.
    
    **Parameters**:
    
    `lc_v`: class: list of 'Lightcurve'-objects        
        The light curves to be used.
            
    `rms_v`: np.ndarray     
        Absolute/fractional rms mulitplied with the spectra.
        
    `rms_err_v`: np.ndarray     
        The error in (absolute/fractional) rms mulitplied with the spectra.
        
    `channel_to_kev`: np.ndarray     
        Conversion from channel (index) to energy (keV).
    
    **Returns**:
    
    `rms_list_channels`: np.ndarray <br>
        With rms for each channel.
    
    `rms_list_channels_err`: np.ndarray      
        With rms error for each channel.
        
    `overlap`: int, optional, default: 0 
        The number of channels that are in two adjacent energy bands.
    """
    
    # Convert to dictionary, where keys = Emax of corresponding channel and values = channels
    if isinstance(channel_to_kev,np.ndarray):
        channel_to_kev_dict = {k: v for v,k in enumerate(channel_to_kev)}
    elif isinstance(channel_to_kev_dict,dict):
        pass
    else:
        print('Channel_to_kev is neither a np.ndarray nor a dictionary.')
        return
    
    # Fill all channels with corresponding rms
    Emin = lc_v[0].Emin 
    Emax_v = [lc.Emax for lc in lc_v]
    
    assert np.all(Emax_v==sorted(Emax_v)), 'The energies are not in increasing order.'
    
    rms_list_channels, rms_list_channels_err = np.zeros(len(channel_to_kev_dict)), np.zeros(len(channel_to_kev_dict))
    grouping = -np.ones(len(channel_to_kev_dict))

    # Where does first Eband start?
    if Emin == 0:
        start_index = 0
    else:
        start_index = channel_to_kev_dict[Emin]
    
    # Go over all energy bands 
    for i,Emax in enumerate(Emax_v): 
        grouping[start_index] = 1
        end_index = channel_to_kev_dict[Emax]+1
        rms_list_channels[start_index:end_index] = rms_v[i] if not math.isnan(rms_v[i]) else 0 # if rms = NaN then put rms-values for these channels to 0
        rms_list_channels_err[start_index:end_index] = rms_err_v[i] if not math.isnan(rms_v[i]) else 0
        # Update start index
        start_index = end_index-overlap
        
    return rms_list_channels, rms_list_channels_err, grouping

def subtract_overlapping_energybands(lc_v):
    """
    Check if two light curve objects overlap (in terms of energy bands) and subtract one 
    from the other if one lies in the other.   
    
    Reason: When the cross spectral properties is being calculated for an energy channel inside another, 
    the lc that resides within another is subtracted from the other. The reasoning behind this is 
    that if one lc is duplicated in the other lc, the error contribution for that channel will 
    not cancel and will contaminate the result.
    
    **Parameters**:

    `lc_v`: list of two 'Lightcurve'-objects    
        The light curves to be used in the covariance computation. 
        lc_v[0] = 1st lightcurve-object, lc_v[1] = 2nd lightcurve-object.
            
    **Returns**:
    
    `lc_v`: list of deep copies of the two 'Lightcurve'-objects (as to not affect the original lightcurves)     
        If one lc was inside another (energy wise) this has been corrected for.
    """
    
    assert len(lc_v) == 2, "You can only compare two light curves."
    
    lc_X = copy.deepcopy(lc_v[0])
    lc_Y = copy.deepcopy(lc_v[1])
    
    # X entirely within Y
    if lc_X.Emin >= lc_Y.Emin and lc_X.Emax <= lc_Y.Emax: 
        print('1st lightcurve-object with Emin-Emax = {}-{} keV lies within 2nd lightcurve-object with Emin-Emax = {}-{} keV'.format(lc_X.Emin,lc_X.Emax,lc_Y.Emin,lc_Y.Emax))
        lc_Y.rate -= lc_X.rate
        lc_Y.err = np.sqrt(lc_Y.err**2-lc_X.err**2)
        lc_Y.R = np.mean(lc_Y.rate)
        print('Removed rate and err of 1st lightcurve-object from the 2nd lightcurve-object.\n')
        
    # Y entirely within X
    if lc_X.Emin <= lc_Y.Emin and lc_X.Emax >= lc_Y.Emax: 
        print('2nd lightcurve-object with Emin-Emax = {}-{} keV lies within 1st lightcurve-object with Emin-Emax = {}-{} keV'.format(lc_Y.Emin,lc_Y.Emax,lc_X.Emin,lc_X.Emax))
        lc_X.rate -= lc_Y.rate
        lc_X.err = np.sqrt(lc_X.err**2-lc_Y.err**2)
        lc_X.R = np.mean(lc_X.rate)
        print('Removed rate and err of the 2nd lightcurve-object from the 1st lightcurve-object.\n')
        
    return [lc_X,lc_Y]

def frs2pha(spectral_data,FRS,FRS_err,grouping,save_path):
    """
    Create a .pha-file to be read by XSPEC.      
    - Important to have the correct HEADERs: https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node6.html  
    - Important to have the correct names for the data: https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node7.html  
    
    **Parameters**:   
    
    `spectral_data`: dict    
        Spectral data extracted from energy spectra with extract_fits(). E.g.: spectral_data = extract_fits(filename_s,keywords=['DATE','EXPOSURE','TELESCOP','INSTRUME'],p=0) 
    
    `FRS`: np.ndarray       
        Frequency resolved spectrum, i.e. fractional variance in freq band multiplied with energy spectra. Units: counts
     
    `FRS_err`: np.ndarray       
        The corresponding error to the Frequency resolved spectrum.
     
    `grouping`: np.ndarray       
        The channel that is the start of a new energy band has "-1", o/w "1".
        
    `save_path`: str       
        Filename including the path to saving directory. E.g.: save_path = "../../Data/FrequencyResolvedSpectra/MAXIJ1535_571/{}_{}to{}Hz.pha".format(obs_id,freq_low,freq_high)    
    """
    
    hdu1 = fits.PrimaryHDU()

    hdu1.header['DATE'] = (datetime.today().strftime('%Y-%m-%d, %H:%M:%S'), 'Creation date (YYYY-MM-DD, hh:mm:ss CET)')

    col1 = fits.Column(name='CHANNEL', format='I', array=spectral_data['CHANNEL'])
    col2 = fits.Column(name='COUNTS', format='J', unit='count', array=FRS)
    col3 = fits.Column(name='STAT_ERR', format='E', unit='count', array=FRS_err)
    col4 = fits.Column(name='QUALITY', format='I', array=np.zeros(len(grouping)))
    col5 = fits.Column(name='GROUPING', format='I', array=grouping)

    coldefs = fits.ColDefs([col1, col2, col3, col4, col5])

    pha_head = fits.Header()

    pha_head['EXTNAME'] = 'SPECTRUM'
    pha_head['DATE'] = (spectral_data['DATE'], 'energy spectrum from (YYYY-MM-DDThh:mm:ss UT)')
    pha_head['EXPOSURE'] = (spectral_data['EXPOSURE'], 'Exposure time (s)')
    pha_head['TELESCOP'] = (spectral_data['TELESCOP'], 'mission/satellite name')
    pha_head['INSTRUME'] = (spectral_data['INSTRUME'], 'instrument/detector name')
    pha_head['BACKFILE'] = ('NONE    ', 'associated background filename')
    pha_head['BACKSCAL'] = (1.000000e+00, 'background file scaling factor')
    pha_head['CORRFILE'] = ('NONE    ', 'associated correction filename')
    pha_head['CORRSCAL'] = (1.000000e+00, 'correction file scaling factor')
    pha_head['RESPFILE'] = ('NONE    ', 'associated redistrib matrix filename')
    pha_head['ANCRFILE'] = ('NONE    ', 'associated ancillary response filename')
    pha_head['AREASCAL'] = (1.000000e+00, 'area scaling factor ')
    pha_head['HDUCLASS'] = ('OGIP    ', 'format conforms to OGIP standard')
    pha_head['HDUCLAS1'] = ('SPECTRUM', 'PHA dataset (OGIP memo OGIP-92-007)')
    pha_head['HDUVERS'] = ('1.1.0   ', 'Obsolete - included for backwards compatibility')
    pha_head['POISSERR'] = (False, 'Poissonian errors not applicable')
    pha_head['CHANTYPE'] = (spectral_data['CHANTYPE'], 'channel type (PHA, PI etc)')
    pha_head['DETCHANS'] = (len(spectral_data['CHANNEL']), 'total number possible channels')

    pha_data = fits.BinTableHDU.from_columns(coldefs,header=pha_head)
    phafile = fits.HDUList([hdu1, pha_data])
    phafile.writeto(save_path,overwrite=True)
    
def merge(lc_v):
    """
    Merge lightcurves into one light curve object.
    
    **Parameters**:
    
    ``lc_v``: list of light curve objects
        At least 2 light curves to be merged.
    """
    
    
    assert len(lc_v) >= 2, 'You need at least two light curves to be merged.'
    
    # Create a new light curve object
    lc_broad = lightcurve('')
    
    # Initilaize 
    print('Initializing with lc in {}-{} keV'.format(lc_v[0].Emin,lc_v[0].Emax))
    lc_broad.rate = np.copy(lc_v[0].rate)
    lc_broad.err = np.copy(lc_v[0].err)
    lc_broad.N = np.copy(lc_v[0].N)
    lc_broad.dt = np.copy(lc_v[0].dt)
    lc_broad.t = np.copy(lc_v[0].t)

    # Merge with the rest
    for i in range(1,len(lc_v)):
        print('Merging with lc in {}-{} keV'.format(lc_v[i].Emin,lc_v[i].Emax))
        assert lc_v[i].Emin >= lc_v[i-1].Emax, "No overlapping energy bands allowed"
        lc_broad.rate += lc_v[i].rate
        lc_broad.err = np.sqrt(lc_broad.err**2+lc_v[i].err**2)
            
    # Update again 
    lc_broad.R = np.mean(lc_broad.rate)
    lc_broad.Emin = lc_v[0].Emin
    lc_broad.Emax = lc_v[-1].Emax
    
    #lc_broad.err = np.sqrt(lc_broad.rate/lc.dt)

    return lc_broad

def merge_energies_lc(lc_v):
    """
    Merge light curves of different energies with each other.
    
    **Parameter**:
    
    ``lc_v``: list of lightcurve objects   
        Light curves to be merged.
        
    **Return**:
    
    ``lc_broad``: lightcurve object   
        The merged light curve, now spanning a larger energy range.
    """
    
    
    assert len(lc_v) >= 2, 'You need at least two light curves to be merged.'
    
    # Create a new light curve object
    lc_broad = lightcurve('')
    
    # Initilaize 
    print('Initializing with lc in {}-{} keV'.format(lc_v[0].Emin,lc_v[0].Emax))
    lc_broad.rate = np.copy(lc_v[0].rate)
    lc_broad.err = np.copy(lc_v[0].err)
    lc_broad.N = np.copy(lc_v[0].N)
    lc_broad.dt = np.copy(lc_v[0].dt)
    lc_broad.t = np.copy(lc_v[0].t)

    # Merge with the rest
    for i in range(1,len(lc_v)):
        print('Merging with lc in {}-{} keV'.format(lc_v[i].Emin,lc_v[i].Emax))
        assert lc_v[i].Emin >= lc_v[i-1].Emax, "No overlapping energy bands allowed"
        lc_broad.rate += lc_v[i].rate
        lc_broad.err = np.sqrt(lc_broad.err**2+lc_v[i].err**2)
            
    # Update again 
    lc_broad.R = np.mean(lc_broad.rate)
    lc_broad.Emin = lc_v[0].Emin
    lc_broad.Emax = lc_v[-1].Emax
    
    #lc_broad.err = np.sqrt(lc_broad.rate/lc.dt)

    return lc_broad

def split_time_lc(lc,equal=True,m=None,step=None,i=None,start=None,stop=None):
    """
    Split light curve into several shorter parts.
    
    **Parameters**:
    
    ``equal``: boolean, default, True
        If true, split light curve into equally long parts.
        If false, split light curve according to start and stop indices. 
        
    ``m,step,i``: ints
        Bins per segment, segments per part and part number respectively.
    
    ``start,stop``: ints
        Index for start and stop of the part.
        
    **Return**:
    
    ``lc_part``: lightcurve object   
        A part of the original light curve. 
    """
    
    if equal:
        start, stop = i*m*step, (i+1)*m*step
    
    print('New part for time indices: ({})-({})'.format(start,stop))
    
    lc_part = copy.deepcopy(lc)
    lc_part.t = lc_part.t[start:stop]
    lc_part.rate = lc_part.rate[start:stop]
    lc_part.err = lc_part.err[start:stop]
    lc_part.R = np.mean(lc_part.rate)
    lc_part.N = len(lc_part.t)

    return lc_part

def find_where_to_split_lc(lc_v,stops,m=2**13):
    """
    If want to split up an observation in several parts.
    
    Needs to be updated: should just need to input one lightcurve, not a list.
    
    **Parameters**:
    
    ``lc_v``: list of light curve objects.
        The lightcurves to be split., 
        
    **Returns**:
    
    ``start_v,stop_v``: list of arrays
        Vectors containing the indices for start and stop for each part.   
    """
    
    ax = standard_plot()
    plt.plot(lc.t,lc.rate)
    N = lc_v[0].N

    stops.append(lc_v[0].t[-1]) #the last part
    start = 0
    
    # Determine where to split parts
    break_points = []
    len_of_parts = []
    for stop in stops:
        part = [t for t in lc_v[0].t if start < t < stop]
        ax.axvline(part[-1],color='r',label='end of a part')
        break_points.append(np.where(lc_v[0].t == part[-1])[0][0])
        
        start = stop
        
        len_of_parts_temp = len(part)/m
        len_of_parts.append(len_of_parts_temp)
        print('\nLength of part: ',len_of_parts_temp)
    
    # Determine start and stop for each part
    break_points.insert(0,0)
    
    start_v = []
    stop_v = []
    for i in range(0,len(break_points)-1):
        #if len_of_parts[i] >= 10:
        start_v.append(break_points[i])
        stop_v.append(break_points[i+1])

    return start_v, stop_v 

def print_datetime_UT(lc_v,obs_start,stops):
    """
    If want to split up an observation in several parts, each of length m*step bins, 
    where m=bins/seg and step = num of seg.
    
    **Parameters**:
    
    ``obs_start``: str, 
        Format: YYYY-MM-DDThh:mm:ss.sss 
    """
    
    obs_time = dati.fromisoformat(obs_start)
    print('obs start time = ',obs_time,'\n')
    
    print('The different parts are:')
    
    start = 0
    for stop in stops:
    
        time_change_to_start = datetime.timedelta(seconds=start)
        time_change_to_end = datetime.timedelta(seconds=stop)
        
        part_start = obs_time + time_change_to_start
        part_end = obs_time + time_change_to_end
        
        print(str(part_start.isoformat(sep='T', timespec='milliseconds'))+',',part_end.isoformat(sep='T', timespec='milliseconds'))

        start = stop

