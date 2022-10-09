#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from .ImportPackages import *
from .AuxiliaryFunctions import *
from .Lightcurve import *
from .Powerspectrum import *

def cross_spec(data,m,percent_limit=90,noise=None,return_noise_wo_sub=False):
    """
    Compute the cross spectrum.
    
    **Parameters**:

    `data`: Either lc (1) or ps (2).    
        (1) class: list of two 'Lightcurve'-objects; the light curves to be used in the coherence computation. <br>
        (2) class: list of two 'PowerSpectrum'-objects; the power spectras to be used in the coherence computation.
        
    `m`: int    
        Number of time bins per segment.
        
    `percent_limit`: float, optional, default: 90    
            Lower limit for percent of time bins having to be filled, i.e. without gaps, for that segment to be used.
    
    `noise`: {'Poisson','Gaussian'}, optional, default: None.    
        For a light curve with 'Poisson'/'Gaussian' errors.
        
    `return_noise_wo_sub`: boolean, optional, default: False    
        If True, noise is returned but not subtracted from powspec.
        If False, returned and subtracted.
    
    **Returns**:

    `ps_v`: class: list of two 'PowerSpectrum'-objects    
        The power spectra used in the coherence computation.
    
    `C_v`: np.ndarray    
        The cross-spectrum (see e.g. section 3 of Epitropakis, A. (2017)) after averaging over K segments.
    """
    
    print('Computing the cross spectrum...\n')

    if isinstance(data[0], lightcurve) and isinstance(data[1], lightcurve):
    
        # Segment-wise: 
        ps_v = [PowerSpectrum(lc, m=m, normalization=None, noise=noise, percent_limit=percent_limit,                               timer_on=False, return_noise_wo_sub=return_noise_wo_sub) for lc in data]

    if isinstance(data[0], PowerSpectrum) and isinstance(data[1], PowerSpectrum):
        ps_v = data
        
    S_v = np.array([ps_v[0].fft_rate_v,ps_v[1].fft_rate_v])
    C_v = np.array([np.conjugate(S1)*S2 for S1,S2 in zip(S_v[0],S_v[1])])

    print('Cross spectrum computed.')    
    
    return ps_v, C_v

def coherence_noiseless(data,m,percent_limit=90,C_v=None):
    """
    Compute the coherence function when noise is not present in the signals.
    
    **Parameters**:

    `data`: Either lc (1) or ps (2).    
        (1) class: list of two 'Lightcurve'-objects; the light curves to be used in the coherence computation. <br>
        (2) class: list of two 'PowerSpectrum'-objects; the power spectras to be used in the coherence computation.
        
    `m`: int    
        Number of time bins per segment.
    
    `percent_limit`: float, optional, default: 90    
            Lower limit for percent of time bins having to be filled, i.e. without gaps, for that segment to be used.
    
    `C_v`: np.ndarray    
        If the cross spectrum has already been found, it can be used. In this case, we require data to be a list of PowerSpectrum-object.
    
    **Returns**:

    `xf`: np.ndarray    
        The frequency-vector (with the current binning).
        
    `gamma2`: np.ndarray    
        The coherence function.
     
    `delta_gamma2`: np.ndarray    
        One sigma uncertainty in the coherence function. 
    """

    start = timeit.default_timer()
    print('Computing the noiseless coherence...')
    
    if isinstance(data[0], lightcurve) and isinstance(data[1], lightcurve):
        # Need to subtract if one lightcurve lies within the other
        lc_v = subtract_overlapping_energybands(data)
    
    # Find the Cross Spectrum
    if not isinstance(C_v, np.ndarray):
        ps_v, C_v = cross_spec(data, m=m, percent_limit=percent_limit)
        S_v = np.array([ps_v[0].fft_rate_v,ps_v[1].fft_rate_v])
        xf = ps_v[0].xf
    
    # Have already found the Cross Spectrum
    elif isinstance(C_v, np.ndarray):
        print('Cross spectra already found.')
        if isinstance(data[0], PowerSpectrum) and isinstance(data[1], PowerSpectrum):
            xf = data[0].xf
            S_v = np.array([data[0].fft_rate_v,data[1].fft_rate_v])
        else:
            print('If you provide me with the cross spectra C_v, the input data need to be a list of PowerSpectrum-objects.')
    
    # The Periodograms w/o normalization
    P1_v = [np.abs(x)**2 for x in S_v[0]]
    P2_v = [np.abs(x)**2 for x in S_v[1]]

    # The Coherence and its error
    K = np.size(C_v,axis=0) #number of segments
    gamma2 = abs(np.mean(C_v,axis=0))**2/(np.mean(P1_v,axis=0)*np.mean(P2_v,axis=0))
    delta_gamma2 = np.sqrt(2)*(1-gamma2)/(np.sqrt(abs(gamma2))*np.sqrt(K))
    
    time_taken = timeit.default_timer()-start
    print('Noiseless coherence found (in {:.2f} sec).'.format(time_taken))
    
    return xf, gamma2, delta_gamma2, C_v

# See Appendix A of (Epitropakis,2017)

def compute_coherence_intrinsic(lc_v,m,noise,percent_limit,output=False):
    """
    Compute the intrinsic coherence for a given number of bins per segment (m).
    
    **Parameters**:
 
    `lc_v`: class: list of two 'Lightcurve'-objects    
        The light curves to be used in the coherence computation.
        
    `m`: int    
        Number of bins per segment.
    
    `noise`: {'Poisson','Gaussian'}.    
        For a light curve with 'Poisson'/'Gaussian' errors. 
    
    `percent_limit`: float
            Lower limit for percent of time bins having to be filled, i.e. without gaps, for that segment to be used.
    
    `output`: boolean    
        If True, plot the whole light curve as well as the "conditions_for_useful_est_int_coh met"-figures. Also, helpful statements are printed.
    
    **Returns**:

    `xf`: np.ndarray    
        The frequency-vector (which depends on m).
        
    `gamma2`: np.ndarray    
        The coherence function.
    
    `delta_gamma2`: np.ndarray    
        One sigma uncertainty in the coherence function. 
        
    `C_mean`: np.ndarray    
        The cross spectra, can be used for time lag estimations.
    """
    
    # Compute power spectra, cross spectra, and signal powers
    ps_v, C_v = cross_spec(lc_v,m,percent_limit,noise=noise,return_noise_wo_sub=True)
    S_v = np.array([ps_v[0].fft_rate_v,ps_v[1].fft_rate_v])

    # Need to extract freq.vector (same for ps_v[0] and ps_v[1]), number of segments K and mean cross spectra
    xf = ps_v[0].xf 
    K = len(C_v)
    C2_mean = np.array(np.abs(np.mean(C_v,axis=0))**2 )
    
    # Compute average for all power terms over all segments and then calculate n2 once
    P1_mean = np.mean(np.abs(S_v[0])**2,axis=0) #power of signal: |S|^2
    P2_mean = np.mean(np.abs(S_v[1])**2,axis=0)
    N1_mean = np.mean(ps_v[0].Pnoise_v,axis=0) #power of noise
    N2_mean = np.mean(ps_v[1].Pnoise_v,axis=0)
    n2 = (P1_mean*N2_mean + P2_mean*N1_mean - N1_mean*N2_mean)/(K)
    
    # Compute coherence
    gamma2 = comp_gamma2(C2_mean,n2,P1_mean,N1_mean,P2_mean,N2_mean)
    
    # Compute coherence error
    delta_gamma2_int = compute_delta_gamma2_int(gamma2,C2_mean,n2,P1_mean,N1_mean,P2_mean,N2_mean,K)
    
    # Check for what frequencies the intrinsic coherence can be usefully estimaed.
    upper_freq_lim = conditions_for_useful_est_int_coh(xf,gamma2,C2_mean,P1_mean,P2_mean,N1_mean,N2_mean,n2,K,output)
    
    return xf, gamma2, delta_gamma2_int, upper_freq_lim
    
def comp_gamma2(C2,n2,P1,N1,P2,N2):    
    """
    Compute gamma2. See: Epitropakis, A. (2017), Eq. (A1). Compare with Eq. (3) (the noiseless case).
    """
    
    return (C2-n2)/((P1-N1)*(P2-N2))
    
    
def compute_delta_gamma2_int(gamma2,C2,n2,P1,N1,P2,N2,K):
    """
    Compute the error delta_gamma2_int. See: Epitropakis, A. (2017), Eq. (A3)
    """
    
    delta_gamma2 = np.sqrt(2/K)*(1-gamma2)/(np.sqrt(abs(gamma2)))
    term1 = 2*K*(n2/(C2-n2))**2 #2*n2**2*K/(C2-n2)**2
    term2 = (N1/(P1-N1))**2 #N1**2/(P1-N1)**2
    term3 = (N2/(P2-N2))**2 #N2**2/(P2-N2)**2
    term4 = K*(delta_gamma2/gamma2)**2 #K*delta_gamma2**2/gamma2**2
    
    return gamma2/np.sqrt(K)*(term1+term2+term3+term4)**(1/2)
    
    
def coherence_intrinsic(lc_v,m_init,noise='Poisson',return_jointly=True,percent_limit=90,output=False,k_lowlim=8):
    """
    Compute the coherence function (Cf) when noise is present in the signals. The Cf is computed at least
    twice, once for m (number of bins per segment) being very high to cover low frequencies and 
    once for m very small (256 as smallest) to cover high frequencies. 
    
    **Parameters**:

    `lc_v`: class: list of two 'Lightcurve'-objects    
        The light curves to be used in the coherence computation.
        
    `m_init`: int    
        Number of bins per segment to start with. m_init will be lowered during each interation. 
        The higher m, the less number of segments. The fewer segments, the lower frequencies will 
        be considered, as f_min = 1/(m_init*dt).
       
    `noise`: {'Poisson','Gaussian'}, optional, default: 'Poisson'.    
        For a light curve with 'Poisson'/'Gaussian' errors. 
       
    `return_jointly`: boolean    
        If True, the coherence computed for different m:s are merged. If False, the coherence 
        is returned as a list of np.ndarrays, so that the computation for different m:s can be compared.
        
    `percent_limit`: float, optional, default: 90    
        Lower limit for percent of time bins having to be filled, i.e. without gaps, for that segment to be used.
    
    `output`: boolean    
        If True, plot the whole light curve as well as the "conditions_for_useful_est_int_coh met"-figures. Also, helpful statements are printed.
    
    `k_lowlim`: int, optional, default: 8    
        The lowest permitted m is 2**k_lowlim. 
    
    **Returns**:

    `xf_v`: np.ndarray (or a list of np.ndarrays if return_jointly is False)    
        The frequency-vector.
        
    `gamma2_v`: np.ndarray (or a list of np.ndarrays if return_jointly is False)    
        The coherence function.
    
    `delta_gamma2_v`: np.ndarray (or a list of np.ndarrays if return_jointly is False)    
        One sigma uncertainty in the coherence function. 
    """
    
    assert isinstance(lc_v[0], lightcurve) and isinstance(lc_v[1], lightcurve), 'The data-input does not contain two light curve objects.'
    assert len(lc_v)==2, 'Coherence should be computed with 2 light curves, not {}'.format(len(lc_v))
    
    print('---------------------------------------------------------------------------------------------------')
    print('                           Computing the intrinsic coherence...')
    print('---------------------------------------------------------------------------------------------------\n')
    start = timeit.default_timer()
    
    # Need to subtract if one lightcurve lies within the other
    lc_v = subtract_overlapping_energybands(lc_v)
    
    # Needed parameters
    m = m_init
    k = int(np.log2(m))
    dt = lc_v[0].dt
    find_again = True
    xf_v, gamma2_v, delta_gamma2_int_v, upper_freq_lim_v = [], [], [], []
    
    # Find intrinsic coherence for current m
    i=0
    while find_again:
        i+=1
        print('Iteration {}) Computing using m = {} bins per segment, i.e. f in [{:.3f},{:.3f}]'.format(i,m,1/(m*dt),1/(2*dt)))
        print('---------------------------------------------------------------------------------------------------')
        xf, gamma2, delta_gamma2_int, upper_freq_lim = compute_coherence_intrinsic(lc_v,m,noise,percent_limit,output=output)
    
        if return_jointly: 
            # Don't use the full freq band found; only use freq up to multi_factor*upper_freq_lim, 
            # i.e. slightly under uppper freq lim.
            if k!=k_lowlim:
                multi_factor = 0.9
                if upper_freq_lim != -1:
                    lim = upper_freq_lim*multi_factor
                else:
                    lim = xf[-1]
            else:
                lim = xf[-1]
                
            gamma2 = gamma2[xf < lim]
            delta_gamma2_int = delta_gamma2_int[xf < lim]
            xf = xf[xf < lim]
        
            # Append to lists
            xf_v = np.append(xf_v,xf)
            gamma2_v = np.append(gamma2_v,gamma2) 
            delta_gamma2_int_v = np.append(delta_gamma2_int_v,delta_gamma2_int) 
            upper_freq_lim_v = np.append(upper_freq_lim_v,upper_freq_lim)
        
        # If not return_jointly; can see what freq-range the different m caught.
        else:
            # Rebin directly; only want to get a feeling for the different m
            num = 50
            lim = xf[-1]
            xf,gamma2,delta_gamma2_int = log_rebin(xf,gamma2,delta_gamma2_int,num=num)
            delta_gamma2_int = error_change(delta_gamma2_int)
            
            # Append to lists
            xf_v.append(xf)
            gamma2_v.append(gamma2) 
            delta_gamma2_int_v.append(delta_gamma2_int) 
            upper_freq_lim_v.append(upper_freq_lim)
    
        print('Intrinsic coherence found for f in [{:.3f},{:.3f}] \n'.format(1/(m*dt),lim))
        
        # Change number of bins / segment until next time
        if k == k_lowlim:
            find_again = False
        else:
            k_temp = np.copy(k)
            multi_factor = 10 #discuss what this is... 
            while 2**(k-1) > multi_factor/(dt*upper_freq_lim):
                k -= 1
                if k == k_lowlim:
                    break 
            if k_temp == k: #to make sure we always make k smaller
                k -= 1
            m = 2**k
    
    time_taken = timeit.default_timer()-start
    print('---------------------------------------------------------------------------------------------------')
    print('             Intrinsic coherence found (in {:.2f} sec). return_jointly = {}'.format(time_taken,return_jointly))
    print('---------------------------------------------------------------------------------------------------')
    
    while True:
        to_clear = input("Do you want to clear the standard out from prints [y/n]? ")
        if to_clear not in ("n", "y"):
            print("Not an appropriate choice.")
        else:
            break
    
    if to_clear == 'y':
        clear_output(True)
        
    while True:
        to_plot = input("Do you want to plot the intrinsic coherence [y/n]? ")
        if to_plot not in ("n", "y"):
            print("Not an appropriate choice.")
        else:
            break
        
    if to_plot == 'y': 
        Ebands = ('{}-{}'.format(lc_v[0].Emin,lc_v[0].Emax),'{}-{}'.format(lc_v[1].Emin,lc_v[1].Emax))
        plot_coherence(xf_v, gamma2_v, delta_gamma2_int_v,Ebands=Ebands)
    
    return xf_v, gamma2_v, delta_gamma2_int_v

def conditions_for_useful_est_int_coh(xf,gamma2,C2,P1,P2,N1,N2,n2,K,output):
    """
    The intrinsic coherence can be usefully estimated when the following conditions_for_useful_est_int_coh are met: 
    1) sqrt(C2) > sqrt(n2)
    2) |S_1|^2/|N_1|^2 > 1/\sqrt{m}
    3) |S_2|^2/|N_2|^2 > 1/\sqrt{m}
    """
    
    upper_freq_lim = []
    
    C = np.sqrt(np.abs(C2))
    n = np.sqrt(np.abs(n2))
    
    if output:
        standard_plot(h=8)
        plt.subplot(3,1,1)

        plt.loglog(xf,C,label=r"$|<C(f)>|^2$")
        plt.loglog(xf,n,label=r"$n^2$")
    
    try:
        # Condition 1
        idx = np.argwhere(np.diff(np.sign(C-n))).flatten()
        upper_freq_lim.append(xf[idx[0]])
        if output:
            plt.plot(upper_freq_lim[-1], n[idx[0]], 'ro',label=str(xf[idx[0]]))
    except:
        pass
    
    if output:
        plt.legend()
        plt.title('Conditions to be met for useful estimation of $\gamma_I^2$')

        plt.subplot(3,2,3)
        plt.title(r'High Power (Term 2 of $\delta \gamma^2_{int}$)')
        plt.loglog(xf,(P1-N1),label=r'$|S_1|^2$')
        plt.loglog(xf,N1*np.ones(np.size(xf))/np.sqrt(K),label=r'$N_1/\sqrt{m}$')
    
    try:
        # Condition 2
        temp = N1*np.ones(np.size(xf))/np.sqrt(K)
        idx = np.argwhere(np.diff(np.sign(temp-(P1-N1)))).flatten()
        upper_freq_lim.append(xf[idx[0]])
        
        if output:
            plt.plot(upper_freq_lim[-1], temp[idx[0]], 'ro',label=str(xf[idx[0]]))
    except:
        pass
    
    if output:
        plt.legend()

        plt.subplot(3,2,4)
        plt.title(r'High Power (Term 3 of $\delta \gamma^2_{int}$)')
        plt.loglog(xf,(P2-N2),label=r'$|S_2|^2$')
        plt.loglog(xf,N2*np.ones(np.size(xf))/np.sqrt(K),label=r'$N_2/\sqrt{m}$')
        
    try:
        # Condition 3
        temp = N2*np.ones(np.size(xf))/np.sqrt(K)
        idx = np.argwhere(np.diff(np.sign(temp-(P2-N2)))).flatten()
        upper_freq_lim.append(xf[idx[0]])
        if output:
            plt.plot(upper_freq_lim[-1], temp[idx[0]], 'ro',label=str(xf[idx[0]]))
    except:
        pass
    
    if output:
        plt.legend()

        plt.subplot(3,1,3)
        plt.title('High Coherence')
        plt.loglog(xf,gamma2,label=r'$\gamma^2_{int}$')
        plt.loglog(xf,n**2/(P1*P2),label=r'$n^2/(P_1P_2)$')
    try:
        idx = np.argwhere(np.diff(np.sign(gamma2-n**2/(P1*P2)))).flatten()
        upper_freq_lim.append(xf[idx[0]])
        if output:
            plt.plot(upper_freq_lim[-1], gamma2[idx[0]], 'ro',label=str(xf[idx[0]]))
    except:
        pass
    
    if output:
        plt.legend()
        plt.tight_layout()
        plt.show()
    
    if len(upper_freq_lim) == 0:
        return -1
    else:
        return np.amin(upper_freq_lim)
    
def plot_coherence(xf,gamma2,delta_gamma2_int,Ebands=None,err_lim=1,num=75,save_fig = False):
    """
    Plot the coherence.
    
    **Parameters**:

    `xf`: np.ndarray    
        The frequency-vector (which depends on m).
        
    `gamma2`: np.ndarray    
        The coherence function.
    
    `delta_gamma2`: np.ndarray    
        One sigma uncertainty in the coherence function. 
     
    `Ebands`: tuple of strings, (Eband1,Eband2), optional, default: None    
        The energy bands of the two light curves compared.
    
    `err_lim`: float, optional, default: 1    
        Error limit. If an error element is larger than err_lim, it is set to zero. 
    
    `num`: int, optional, default: 75    
        The number of logarithmic frequency bins to create before plotting.
        
    `save_fig`: boolean, optional, default: False    
         If True, asks for path to location to save plot.
    """
    
    standard_plot()
    ax = plt.gca()

    if Ebands == None:
        Eband1 = input('First energyband in keV [?-?]): ')
        Eband2 = input('Second energyband in keV [?-?]: ')
    else:
        Eband1, Eband2 = Ebands[0], Ebands[1]
    ax.text(0.3,0.1,'({} keV) vs. ({} keV)'.format(Eband1,Eband2),fontsize=14,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    
    if len(xf) <= 3:
        for i in range(0,len(xf)):
            ax.errorbar(xf[i],gamma2[i],yerr=delta_gamma2_int[i], fmt = '.',mfc='w',capsize=2, elinewidth=1, markeredgewidth=1,label=r'$N=2^6$')
    else:
        xf, gamma2, error = log_rebin(xf,gamma2,delta_gamma2_int,num=num)
        error = error_change(error,err_lim)
        ax.errorbar(xf,gamma2,yerr=error, fmt = '.k',mfc='w',capsize=2, elinewidth=1, markeredgewidth=1,label=r'$N=2^6$')
    
    ax.set_ylim([0,1.4])
    ax.axhline(1,color='k',linewidth=1,alpha=0.7)
    ax.set_xscale("log")
    plt.xlabel('Frequency [Hz]')

    plt.ylabel('$\gamma_{int}^2$')
    if save_fig:
        path = input('Path to location to save plot [end with .png]: ') 
        plt.savefig(path,bbox_inches='tight')
    plt.show()

