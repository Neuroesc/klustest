# klustest
Klustest, cluster analysis tool in Matlab.

Following automated or manual spike sorting using Klustakwik (see [kwiktint](https://github.com/Neuroesc/kwiktint)) and manual curation and refinement (see Tint), this Matlab function loads the spike sorted data, performs analyses (see below), creates a summary figure and saves an output in Matlab table format.

<img width="3627" height="2107" alt="Picture1" src="https://github.com/user-attachments/assets/e8478c3a-23c5-4caa-a2d5-72b273505dc3" />

> [!NOTE]
> This function is intended to run after spike sorting using [kwiktint](https://github.com/Neuroesc/kwiktint).

# Outputs
This function performs spatial analyses and creates spatial modulation visualisation plots:
- Plots positions and spikes
- Plots the firing rate map (see [ratemapper](https://github.com/Neuroesc/rate_mapper))
  - calculates [spatial information content](https://doi.org/10.1002%2F%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K), [sparsity](https://doi.org/10.1002%2F%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K), [coherence](https://dx.doi.org/10.1523%2FJNEUROSCI.1704-07.2007)
- Plots the spatial autocorrelation
  - calculates [grid score](https://doi.org/10.1126/science.1188210), wavelength, orientation

Waveform plots and quantification:
- Plots the mean and standard deviation waveforms and a random sample
  - calculates width of waveform
  - all waveforms are converted to volts

Cluster quality plots and quantification:
 - Mahalanobis distance (within- vs between-cluster distances)
   - calculates [isolation distance and l-ratio](https://doi.org/10.1016/j.neuroscience.2004.09.066)
 - Feature space (feature 1 vs 2 for highest amplitude channel)

Head direction plots and quantification:
- Linear head direction and dwell time plot
  - calculates [Rayleigh vector length](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics), preferred firing direction
- Polar head direction plot and dwell time
- Split, half-session, linear head direction plots 

Spike-time analyses and visualisation:
- Plots 50 ms inter-spike-interal
  - calculates ISI full width to half maximum
  - calculates [burst index](https://doi.org/10.1002/hipo.22002)
- Plots 50 ms higher order inter-spike-interval (i.e. spike autocorrelograms)
  - calculates proportional [refractory period violations](https://doi.org/10.1152/jn.00699.2015)
- Plots 500 ms higher order inter-spike-interval (i.e. spike autocorrelograms)
- Plots theta-phase spike relationship
  - calculates calculates [Rayleigh vector length](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics) (theta modulation)
  - calculates preferred theta phase
- Raster plot

Movement modulation plots and quantification:
- Plots running speed vs firing rate linear relationship
  - calculates [speed score, slope, y-intercept](https://doi.org/10.1038/nature14622)
- Plots angular head velocity vs firing rate

Local field potential:
- Plots the frequency of spikes vs theta phase
  - calculates Rayleigh vector length (theta modulation)
  - calculates preferred theta phase

# Usage
process files with default settings:
```
klustest()
```
process with Name-Value pairs used to control aspects of the cluster cutting:
```
kwiktint(Name,Value,...)
```
<img width="1201" height="671" alt="image" src="https://github.com/user-attachments/assets/fa4d86a3-f3c2-4cf3-8326-f6eeec5be179" />

# Parameters
'tetrodes'      -   Numeric vector that specified the tetrodes to run on (i.e. [1 2 3 4 5 6] would run on tetrodes 1 to 6). Default value is 1:16

'clusters'      -   Numeric vector that specifies the clusters to run on. Set to 0 to run on all clusters. Default value is 0

'rname'         -   String or charcter vector specifying the rat name This is used in the sdata structure so its important to ensure this is correct. Default value is the name of the parent directory

'outname'       -   String or charcter vector specifying the output filename used by klustest - i.e. figures will be saved in a folder named this. Default value is 'klustest'.

'cname'         -   String or charcter vector specifying the output filename used by kwiktint - or the name used for klustakwiked files. Default value is 'kwiktint'.

> [!NOTE]
> When merging and splitting sessions, additional information needs to be provided within the klustest code

# Examples
run function using default values
```
klustest()
```
run function using default values, but only on tetrodes 1 and 5
```
klustest('tetrodes',[1 5])
```
