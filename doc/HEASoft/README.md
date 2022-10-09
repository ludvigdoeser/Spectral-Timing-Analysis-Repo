<h1>
  How to work with HEASOFT; XSELECT; XSPEC; nicerl2 etc
</h1>

What's the purpose of NICER? Have a look at: [NICER mission guide](https://sites.astro.caltech.edu/~srk/XC/Notes/NICER_Mission_Guide.pdf).

---

<h1> Download data: </h1>

1. Go to [NICER Archive](https://heasarc.gsfc.nasa.gov/docs/nicer/nicer_archive.html) and click on **[HEASARC Browse interface](https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dnicermastr&Action=More+Options)**. Alternatively, go directly to: [HEASARC Browse](https://heasarc.gsfc.nasa.gov/cgi-bin/W3Browse/w3browse.pl). Search for an object, e.g. GX339-4, and if successful, you'll get a data table:

| ![NICERarchive](NICERarchive.png)|
|-|

2. Click on D under services.

3. Then go to: XTI Cleaned Data (event_cl):

| ![NICERarchive_XTIcleanedEVT](NICERarchive_XTIcleanedEVT.png)|
|-|

4. Click on the *cl.evt.gz-file, which will re-direct you. Copy paste the link to that webpage and run it with `wget` in the terminal, e.g.

```bash
wget https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2021_10/4133010273/xti/event_cl/ni4133010273_0mpu7_cl.evt.gz
```

You can also do `!wget *link*` directly within a jupyter notebook if you want to download it from there.

5. To unzip correctly, run:

```bash
gzip -dk filename.gz
```

6. You can then use the HeaSoft-package to extract **light curves** or **energy spectra**. If you want you can also run the star.extract_fits()-method to get the data:

```python
extract_fits('ni4133010273_0mpu7_cl.evt')
```

<h2>Downloading HeaSoft</h2>

1: [Download the HEASOFT Software](https://heasarc.gsfc.nasa.gov/lheasoft/download.html). At that page, follow steps 1-3. Pick **SOURCE CODE DISTRIBUTION** and pick "All" packages.

2 (on Intel Mac): Make sure that you have installed  "XQuartz" samt "Homebrew"; if not, start with installing those. Then, [install HEASoft - Intel Mac](https://heasarc.gsfc.nasa.gov/lheasoft/macos.html). Follow the "Bourne shell variants (bash/sh/zsh)"-column.

**Some notes:**
* Before running the EXPORT-commands, make sure the files/executables are at the correct paths. E.g. I had to change `export FC=/usr/local/bin/gfortran-11` to `export FC=/usr/local/gfortran/bin/gfortran`.
* You can run `tail -f build.log` to see how the build is doing. Wait for it to finish before running the installation.

---

<h1> XSELECT </h1>

[Documentation for XSELECT](https://heasarc.gsfc.nasa.gov/docs/rosat/ros_xselect_guide/#tth_chAp2).

<h2> In terminal to extract lightcurves: </h2>

**To initialize software:**

Create environment HEADAS and source:

```bash
export HEADAS=/Applications/heasoft-6.29/x86_64-apple-darwin21.1.0  

. $HEADAS/headas-init.sh
```

**To start software, e.g. XSELECT, just run:**

```bash
XSELECT
```

**Then you need to input:**

```
[session name] (does only matter if want to run several parallell sessions)
read events (to tell software to read an event file)
[enter directory] ("./" if you started the software in the correct directory)
[enter filename] (should be on form "*.evt"
yes (to the q: "Reset the mission?")
@script (e.g @../../HeaSoft/lc0p5_5ms.xco)
```

A .xco-script can look like:

```
set binsize 0.005  
set phaname PI

filter pha_cutoff 200 250  
extract curve  
save curve 2to2.5kev.lc    
```

The last 3 lines can then be repeated for all light curves to be extracted

Exit program by:

```
exit
```

Then the program will ask if you want to save the session (if yes, then you can continue that session later).

<h2> In terminal to extract spectrum: </h2>

Like above (for light curves), with the difference that we know need to specify the channels we want to extract the spectrum for (e.g. 0-1500 = 0-15 keV)

```
set phaname PI

filter pha_cutoff 0 1500  
extract spectrum  
save spectrum fullspec.pha   
```

If we want to filter time: skip "set phaname PI" for some reason...

To get the new times, one can use the following function from the package:

```python
print_datetime_UT(lc_v,obs_start,stops)
```

which will print out the times to provide in XSELECT:

```
filter pha_cutoff 0 1500  

filter time ut
[enter start and stop] > yyyy-mm-ddThh:mm:ss.sss, yyyy-mm-ddThh:mm:ss.sss
[enter start and stop] > x #to save and exit "filter time"
extract spectrum  
save spectrum fullspec_partX.pha

clear time all #clears the filter, otherwise you will have two filters next time
# and then repeat from filter time ut to filter next part

# save spectrum fullspec_part1.pha
# save spectrum fullspec_part2.pha
# etc...
```

---

<h1> XSPEC </h1>

Useful links:

* http://polywww.in2p3.fr/activites/physique/glast/workbook/pages/sciTools_latGrbAnalysis/2_runXSPEC_part01.html
* https://www.rri.res.in/~bpaul/asw/html/spec_exe1.html
* http://polywww.in2p3.fr/activites/physique/glast/workbook/pages/sciTools_latGrbAnalysis/example01_xspecHelp.html

<h3> How to make .pha-files correctly?</h3>

<u>Alternative 1</u>, **In python**:

```python
a = np.array([spectral_data['CHANNEL'], FRS, FRS_err])
mat = np.transpose(np.matrix(a))
with open(path-to+filename.txt,'wb') as f:
for line in mat:
    np.savetxt(f, line, fmt='%.0f')
```

Then in bash:

```
ascii2pha filename.txt
```

Note: "Errors present?"" will always have [no] as default, while the others will be remembered until next time, as:

| ![ascii2pha](ascii2pha.png)|
|-|

Then run:

```
grppha filename.pha

group 0 1500 50
```

where filename.pha is the entered output filename from ascii2pha. The group-command then gives two extra columns (quality, grouping) to the .pha-file (a fits-file format).

<u>Alternative 2</u>:

Do everything in **python** using FITS-format there. Simply call:

```
frs2pha(spectral_data,FRS,FRS_err,grouping,save_path)
```

<h3> Model fitting and plotting of spectra</h3>

Start software, by running:

```
XSPEC
```

If it doesn't work, use the EXPORT command as above under XSELECT. Then open a figure:

```
cpd /xw
```

Import the energy spectra:

```
data filename.pha
```

It will tell you that no response has been loaded. So these files need to be downloaded and then we load them, e.g.:

```
resp nixtiref20170601v002.rmf

arf nixtiaveonaxis20170601v004.arf
```

or, in case that resp och arf has been combined into .rsp-file:

```
response _.rsp
```

Then to ignore certain channels:

```
ig 0-50, 1000-**
```

To plot the data:

```
plot ldata
```

To fit a model to the data:

```
model powerl

1:powerlaw:PhoIndex>1.6

2:powerlaw:norm>1
```

The main command to plot is then:

```
plot eeuf
```

If need to override something, e.g. what channels to ignore, just rerun such a command:

```
ig 0-150
```

If want to notice some channels again

```
not 100-150
```

If we want to import another file, then we can do:

```
data 2:2 test.pha
```

where 2:2 specifies the "channel" of the figure to use. We then need to use the response-files again, but now for this file, which we specify with a "2" for this "channel":

```
resp 2 nixtiref20170601v002.rmf  
```

and similar for the .arf-file. And then a third example:

```
data 3:3 testg.pha

load response files as above...

ig 3: 1-3, 20-**
```

To show all data and all model-parameters:

```
show data

show param
```

To change a param value:

```
new 7 1
```
will change param7 to value 1.

To modify plot you can do:

```
setplot rebin 100 10

setplot energy
```

For simple modeling of X-ray binaries:

```
model wabs*(diskbb+nthcomp)
```

1. diskbb is an extended black body to suit accretion discs. It only has two parameters (the inner temperature, approx 0.5 keV) and the normalization.
2. For nthcomp you need to specify more params:

* Gamma (the slope, same as for a power-law)
* kT_e (electron temperature. Can be "freezed" to 100 keV.)
* kT_bb (temperature of incoming photons. Can be linked to Tin in diskbb.)
* inp_type: let it be 0
* Redshift: let it be 0
* norm: normalization

To link two parameters, run `show par` to see which parameter numbers the relevant parameters have. If we want connect e.g. param 6 to 2, then do `new 6=2`. To freeze a parameter, do `fre N` (with N being the parameter number); this value will then not change during the fit.

---

<h1> CALDB </h1>

CALDB, e.g. for making arf and rmf files: [NICER Responses (ARFs and RMFs)](https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/arf-rmf/).

Start by sourcing CALDB:

```bash
export CALDB=/Applications/heasoft-6.29/caldb
source $CALDB/software/tools/caldbinit.csh
```

Make sure you have **data** for the correct mission and instrument; check at: [How to Install a Calibration Database](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/install.html).

If does not work, just run:

```bash
export CALDB=/Applications/heasoft-6.29/caldb
export CALDBCONFIG=$CALDB/software/tools/caldb.config
export CALDBALIAS=$CALDB/software/tools/alias_config.fits
```

Then you can do:

```
nicerarf fullspec.pha 275.0915000 7.1853889 ni1200120103.mkf ni1200120103.mkf nixtia1200120103.arf outwtfile=nixtia1200120103_wt.lis
```

where

* fullspec.pha is the spectrum generated by XSELECT
* 275.0915000 7.1853889 is RA and DEC in J2000 (can be found in the Hesarc archive, e.g. by clicking at **O** to get RA and DEC in J2000 (is displayed in hrs and deg at the main page))
* ni1200120103.mkf is the filter file, which we can download at the Heasarc archive by going into **D** and then into the auxil-folder
* nixtia1200120103.arf is the output arf-file
* outwtfile=nixtia1200120103_wt.lis is a file needed by the rmf-generation

**NOTE:** according to [NICER Responses (ARFs and RMFs)](https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/arf-rmf/) we should use the cleaned event file ni1200120103_0mpu7_cl.evt instead of the first ni1200120103.mkf. However, according to https://www.youtube.com/watch?v=mRlH17R0bmY and https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/nicerarf.html the .mkf-file should be put there twice..?

Then, for rmf:

```
nicerrmf fullspec.pha ni1200120103.mkf nixtir1200120103.rmf detlist=@nixtia1200120103_wt.lis
```

**(writing some parts of the docs in retrospect, so I don't remember if this worked...)**

---

<h1> nicerl2 </h1>

If no cleaned data - extract cleaned data using nicerl2:

1. Extract "Full Observation Dataset" as tar-file from HeaSarc, then follow [nicerl2 docs](https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/nicerl2.html).

2. Need to export environments to use the CALDB-package:

```bash
export HEADAS=/Applications/heasoft-6.29/x86_64-apple-darwin21.1.0  

. $HEADAS/headas-init.sh

export CALDB=/Applications/heasoft-6.29/caldb
export CALDBCONFIG=$CALDB/software/tools/caldb.config
export CALDBALIAS=$CALDB/software/tools/alias_config.fits
```

and then run:

```bash
nicerl2 indir=OBSERVATION_ID clobber=YES
# e.g. nicerl2 indir=4133010105 clobber=YES
```

where:
* indir = path to directory with all data! E.g. "Full Observation Dataset (1200120224)" (HeaSarc)
* clobber = YES = needed to overwrite files
