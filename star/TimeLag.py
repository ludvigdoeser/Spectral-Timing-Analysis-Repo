#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from .ImportPackages import *
from .AuxiliaryFunctions import *
from .Lightcurve import *
from .Powerspectrum import *
from .Coherence import *

def time_lag(lc_v,m_v,percent_limit=90,auto=False):
    """
    Find the time lag between two light curves. 
    
    **Parameters**:

    `lc_v`: class: list of two 'Lightcurve'-objects    
        The light curves to be used in the coherence computation.
    
    `m_v`: np.ndarray    
        Array of number of segments to use, in increasing order.
        format: m_v[m_1,m_2,...], where type(m_i) = int, preferably a power of 2, and m_1 < m_2 etc.  
        
    `percent_limit`: float, optional, default: 90    
        Lower limit for percent of time bins having to be filled, i.e. without gaps, for that segment to be used.
    
    `auto`: boolean, default: False,
        If True, does not ask for user input (if want to clear all output and display time lag plot).
    
    **Returns**:

    `xf_v`: np.ndarray    
        Frequencies.
    
    `tau_v`: np.ndarray    
        Time lag. IMPORTANT NOTE: positive sign here means that hard leads (soft lags). Note, however,
        that the plot_time_lag-function makes clear what is what. Sorry for the confusion, but didn't want 
        to mess up somewhere by fixing the sign here and forgetting to fix it elsewhere.
    
    `dtau_v`: np.ndarray    
        Error in time lag.
        
    `num`: int
        The number of logarithmic frequency bins to create before plotting.
    """
    
    print('---------------------------------------------------------------------------------------------------')
    print('                           Computing the time lag...')
    print('---------------------------------------------------------------------------------------------------\n')
    start = timeit.default_timer()
        
    assert isinstance(lc_v[0], lightcurve) and isinstance(lc_v[1], lightcurve), 'The lc_v-input does not contain two light curve objects.'
    if len(m_v) >= 2:
        for i in range(0,len(m_v)-1):
            assert m_v[i] < m_v[i+1], 'Make sure that m_v[i] < m_v[i+1] for all i'
    
    # Make sure first lc is the one with lowest energy (the soft one); if not, switch place
    if lc_v[0].Emin > lc_v[1].Emin:
        print('Light curve with lowest energy was placed second in the list. I will switch and put it first.\n')
        lc_temp = copy.deepcopy(lc_v[0])
        lc_v[0] = copy.deepcopy(lc_v[1])
        lc_v[1] = lc_temp
    
    # Need to subtract if one lightcurve lies within the other
    lc_v = subtract_overlapping_energybands(lc_v)
    
    xf_v = np.array([])
    tau_v = np.array([])
    dtau_v = np.array([])
    dt = lc_v[0].dt
    
    freq_uplim = 1/(2*dt)
    for i in range(0,len(m_v)):
        m = m_v[i]
        print('Iteration {}) Computing using m = {} bins per segment, i.e. f in [{:.3f},{:.3f}]'.format(i+1,m,1/(m*dt),freq_uplim))
        print('---------------------------------------------------------------------------------------------------')
        
        # Find cross spectra
        ps_v, C_v = cross_spec(lc_v,m=m,percent_limit=percent_limit)

        # Coherence
        xf, gamma2, delta_gamma2, _ = coherence_noiseless(ps_v,m=m,C_v=C_v)
        
        K = np.size(C_v,axis=0) #number of segments 
        C = np.mean(C_v,axis=0) #mean cross spectra
        
        tau = np.angle(C)/(2*np.pi*xf)
        dtau = np.sqrt((1-gamma2)/(2*gamma2*K))/(2*np.pi*xf) 
        
        if freq_uplim != -1:
            print('Only use freq < {:.3f} to avoid to much overlaping.'.format(freq_uplim))
            tau = tau[xf < freq_uplim]
            dtau = dtau[xf < freq_uplim]
            xf = xf[xf < freq_uplim]
        
        xf_v = np.append(xf_v,xf)
        tau_v = np.append(tau_v,tau)
        dtau_v = np.append(dtau_v,dtau)
        
        # If len(m) >= 2, then we don't want the two computations to overlap too much... 
        # Now we make them overlap a factor 10 larger than the lowest frequency from the small m
        # E.g. if m = [2**8,2**17], dt = 0.002 --> freq_uplim = 1.96*10, meaning that even though 
        # m=2**17 covers f\in[0.0038,249] we only look at f\in[0.0038,19.6]. The high freq interval was
        # taken care of by m = [2**8].
        freq_uplim = np.amin(xf_v)*10
        print('Time lag for m = {} computed. \n'.format(m))
     
    time_taken = timeit.default_timer()-start
    print('---------------------------------------------------------------------------------------------------')
    print('                           Time lag found (in {:.2f} sec).'.format(time_taken))
    print('---------------------------------------------------------------------------------------------------')
        
    while True:
        if not auto:
            to_clear = input("Do you want to clear the standard out from prints [y/n]? ")
        else:
            to_clear = 'y'
        if to_clear not in ("n", "y"):
            print("Not an appropriate choice.")
        else:
            break
    
    if to_clear == 'y':
        clear_output(True)
        
    while True:
        if not auto:
            to_plot = input("Do you want to plot the time lag [y/n]? ")
        else:
            to_plot = 'n'
            
        if to_plot not in ("n", "y"):
            print("Not an appropriate choice.")
        else:
            break
        
    if to_plot == 'y': 
        Ebands = ('{}-{}'.format(lc_v[0].Emin,lc_v[0].Emax),'{}-{}'.format(lc_v[1].Emin,lc_v[1].Emax))
        plot_timelag(xf_v, tau_v, dtau_v,Ebands=Ebands)
        
    return xf_v, tau_v, dtau_v

def plot_timelag(xf_v,tau_v,dtau_v,Ebands=None,num=50,save_fig=False):
    """
    Plot the time lag.
    
    **Parameters**:

    `xf_v`: np.ndarray    
        Frequencies.
    
    `tau_v`: np.ndarray    
        Time lag.
    
    `dtau_v`: np.ndarray    
        Error in time lag.
        
    `Ebands`: tuple of strings, (Eband1,Eband2), optional, default: None    
        The energy bands of the two light curves compared.
    
    `num`: int, optional, default: 50
        The number of logarithmic frequency bins to create before plotting.
       
    `save_fig`: boolean, optional, default: False    
         If True, asks for path to location to save plot.
    """
    
    # Rebin
    xf_rebin,tau,dtau = log_rebin(xf_v, tau_v, dtau_v,num=num)
    
    ax = standard_plot()
    
    # Separate into pos/neg values of tau
    tau_n, tau_p = [-t for t in tau if t<0 ],[t for t in tau if t>0]
    x_n, x_p = [x for x,t in zip(xf_rebin,tau) if t<0],[x for x,t in zip(xf_rebin,tau) if t>0]
    dtau_n, dtau_p = [x for x,t in zip(dtau,tau) if t<0],[x for x,t in zip(dtau,tau) if t>0]
    
    # Loop over neg vs pos lags
    for x,tau,dtau,c,s,l in zip([x_n, x_p],[tau_n, tau_p],[dtau_n, dtau_p],['w','k'],[5,4.5],['soft leads','hard leads']):
        
        # Give label to first
        first = True
        
        # Replace errors that go below 0 with a downarrow
        for i in range(0,len(x)):
            
            if abs(dtau[i]) > abs(tau[i]): #if error is larger than tau-value
                # Plot lower error as an arrow
                ax.errorbar(x[i],tau[i],yerr=[[0.7*tau[i]],[0]], fmt = 'ok', uplims = True, mfc=c, markersize=s, elinewidth=1)
                # Plot upper error as normal
                ax.errorbar(x[i],tau[i],yerr=[[0],[dtau[i]]], fmt = 'ok', mfc=c, markersize=s, capsize=2, elinewidth=1, markeredgewidth=1)
            else:
                pass
                # Plot as normal
                if first:
                    ax.errorbar(x[i],tau[i],yerr=[[dtau[i]],[dtau[i]]], fmt = 'ok', mfc=c, markersize=s, capsize=2, elinewidth=1, markeredgewidth=1,label=l)
                    first = False
                else:
                    ax.errorbar(x[i],tau[i],yerr=[[dtau[i]],[dtau[i]]], fmt = 'ok', mfc=c, markersize=s, capsize=2, elinewidth=1, markeredgewidth=1)

    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.legend()
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Lag (sec.)')
    
    if Ebands == None:
        Eband1 = input('First energyband in keV [?-?]): ')
        Eband2 = input('Second energyband in keV [?-?]: ')
    else:
        Eband1, Eband2 = Ebands[0], Ebands[1]
    ax.text(0.3,0.1,'({} keV) vs. ({} keV)'.format(Eband1,Eband2),fontsize=14,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    # Also good spot to place text:
    #ax.text(0.73,0.93,'({} keV) vs. ({} keV)'.format(first_band,second_band),fontsize=14,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    
    if save_fig:
        path = input('Path to location to save plot [end with .png]: ') 
        plt.savefig(path,bbox_inches='tight')
    plt.show()

