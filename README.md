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
  - calculates spatial information content, sparsity, coherence
- Plots the spatial autocorrelation
  - calculates grid score, wavelength, orientation

Waveform plots and quantification:
- Plots the mean and standard deviation waveforms and a random sample
  - calculates width of waveform
  - all waveforms are converted to volts

Cluster quality plots and quantification:
 - Mahalanobis distance (within- vs between-cluster distances)
   - calculates isolation distance and l-ratio
 - Feature space (feature 1 vs 2 for highest amplitude channel)

Head direction plots and quantification:
- Linear head direction and dwell time plot
  - calculates Rayleigh vector length, preferred firing direction
- Polar head direction plot and dwell time
- Split, half-session, linear head direction plots 

Spike-time analyses and visualisation:
- Plots 50 ms inter-spike-interal
  - calculates ISI full width to half maximum
  - calculates burst index
- Plots 50 ms higher order inter-spike-interval (i.e. spike autocorrelograms)
  - calculates proportional refractory period violations
- Plots 500 ms higher order inter-spike-interval (i.e. spike autocorrelograms)
- Plots theta-phase spike relationship
  - calculates calculates Rayleigh vector length (theta modulation)
  - calculates preferred theta phase
- Raster plot

Movement modulation plots and quantification:
- Plots running speed vs firing rate linear relationship
  - calculates speed score, slope, y-intercept
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
'combine'          -   (default = true) Logical or scalar, set to 1 or true to combine multiple .set files into one output, set to 0 or false to analyse sessions individually. If sessions are to be combined, they should be named in numerically or alphabetically ascending order, matching the order they were recorded. TINT will always order them in this way when they are opened or kkwiked, so for continuity they should be named this way. I name recordings using this convention:
>(date in format yymmdd)(a-z order of recording)_(name of maze)

'tetrodes'         -   (default = 1:16) Vector of tetrodes to be analysed i.e. [1 2 3 4], the function will run on the included tetrodes if they are available, missing tetrodes are ignored

'outname'          -   (default = 'kwiktint') String, the file name to use for combined outputs

'assume_all_set'   -   (default = true) Logical or scalar, set this to 1 or true and the function will always just assume you want to analyse all available .set files, it will not ask you to select them I always separate sessions (i.e. all the recordings related to one data collection) into different directories, so I always want to combine all .set files in a directory. Some people have other conventions like saving all of the recordings for a day in a directory, in which case they would need to specify the files each time

'backup_cuts'      -   (default = true) Logical or scalar, set to 1 or true and kwiktint will backup .cut files if they already exist, the backups are named with the exact date/time they are backed up appended at the end of the extension and saved in a 'kwiktint' directory alongside the data

'max_tp'           -   (default = 3) Scalar specifying how many tetrodes can be analysed simultaneously or how many instances of TINT can be open simultaneously, includes ones opened by the user. This means if you open a copy of TINT manually to do some manual cluster cutting etc klustakwik will only run max_tp-1 copies of TINT.

> [!NOTE]
> Any klustakwik settings that can be personalised in Tint can also be specified in kwiktint.

# Examples
run function using default values
```
kwiktint()
```
run function using default values, but only on tetrodes 1 and 5
```
kwiktint('tetrodes',[1 5])
```
run function using default values, all specified
```
kwiktint('combine',1,'tetrodes',1:16,'outname','kwiktint','assume_all_set',1,'backup_cuts',1,'max_tp',3)
```
