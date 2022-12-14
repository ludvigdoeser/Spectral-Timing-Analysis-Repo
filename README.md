# Spectral Timing Analysis Repo (STAR)
Spectral Timing Analysis Software Package.

Spectral and timing properties of black hole X-ray binaries (BHXRB) evolve significantly  in response to changes in mass accretion rate and the configuration of the accretion flow. Tracking their evolution has therefore played a major role in understanding these systems' geometry and behaviour.

For basic introduction to timing analysis, check out e.g.: [Light Curves and What They Can Tell Us](https://imagine.gsfc.nasa.gov/science/toolbox/timing1.html) (NASA Goddard Space Flight Center, August 2013) 

## Installation and setup
[Install conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) if you don't have conda installed.

Clone repository and then install all dependencies using conda via:

```bash
cd star
conda env create -f environment.yml
```

Activate the environment by invoking:

```bash
conda activate star
```

Next, to initialize the correct paths for the scripts in this project run the setup tools using:

```bash
pip install -e .
```

The `-e` flag here enables you to alter the package in the src/ folder and the star-package will be updated automatically. **NOTE**, however, that if you're in a jupyter notebook and you do modify the source code, then you will have to restart the kernel for the change to take effect.

Or if you don't intend to modify the source code, just run pip without the `-e` flag.

In python, you should now be able to import all modules by:

```python
from star import *
```

## Documentation

Find the path to the repository on your local machine and click on `open with` and choose your browser of preference. If that does not work, you can try accessing my [documentation page](https://ludvigdoeser.github.io/stasp/index.html).

To understand what the code does, please have a look at the summary notes at:

```bash
doc/SpectralTimingAnalysis_Summary.pdf
```

It can also be opened in your browser using:

```bash
file:///*path-to-repo*/doc/SpectralTimingAnalysis_Summary.pdf
```

If it complains about it not being readable, you might have to manually go to finder (on Mac), locate the repository and the .pdf, then click on `open with` and choose your browser of preference.

## Installation of other softwares

Go [here](https://github.com/ludvigdoeser/Spectral-Timing-Analysis/tree/main/doc/HEASoft) for notes on how to work with HeaSoft, XSPEC, XSELECT, nicerl2, setting up CALDB etc.
