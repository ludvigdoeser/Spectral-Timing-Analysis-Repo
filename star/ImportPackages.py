#!/usr/bin/env python
# coding: utf-8

# In[2]:


## Import packages
from astropy.table import Table  
from lightkurve import search_lightcurve
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
from scipy.fft import fft, fftfreq, ifft
from scipy.signal import blackman
from scipy.signal import hann
from astropy.io import fits
from IPython.display import clear_output ; import timeit ; start=timeit.default_timer()
import time
import math
from functools import reduce
from pprint import pprint
import re
import operator
import copy
import pickle 
from itertools import combinations
import datetime
from datetime import datetime as dati

