<img width="303" height="102" alt="neuroesc_long" src="https://github.com/user-attachments/assets/e185a933-e27d-4436-ab1f-52a4fe389e38" />

# klustest
Klustest, cluster analysis tool in Matlab.

Following automated or manual spike sorting using Klustakwik (see [kwiktint](https://github.com/Neuroesc/kwiktint)) and manual curation and refinement (see Tint), this Matlab function loads the spike sorted data, performs analyses (see below), creates a summary figure and saves an output in Matlab table format.

> [!IMPORTANT]
> This function is intended to run after spike sorting using [kwiktint](https://github.com/Neuroesc/kwiktint).

# Contents
- [Specifying data 'parts'](https://github.com/Neuroesc/klustest/edit/main/README.md#specifying-data-parts)
- [Usage](https://github.com/Neuroesc/klustest/edit/main/README.md#usage)
- [Examples](https://github.com/Neuroesc/klustest/edit/main/README.md#examples)
- [Parameters](https://github.com/Neuroesc/klustest/edit/main/README.md#parameters)
- [Figure outputs](https://github.com/Neuroesc/klustest/edit/main/README.md#figure-outputs)
- [Data outputs](https://github.com/Neuroesc/klustest/edit/main/README.md#data-outputs)

# Specifying data 'parts'
This function is assumed to follow spike sorting using [kwiktint](https://github.com/Neuroesc/kwiktint). Kwiktint can merge multiple recording sessions into one output. Klustest can redivide the data back up into indiviual sessions or 'parts' for analysis and visualisation.

Each part is defined by a set of characteristics:
1. Each part needs a **name**, which is used when saving the output. Because these are used as structure field names they cannot start with a number. Common names might be 'session1' or 'part_3' or 'Arena_1'.
> [!NOTE]
> These part names are also used in the [figure](https://github.com/Neuroesc/klustest/edit/main/README.md#figure-outputs) and [data](https://github.com/Neuroesc/klustest/edit/main/README.md#data-outputs) outputs.

2. Each part needs a **method** for partitioning - do we want to split the merged file up according to recordings? Do we want to keep some recordings combined?:  
   - 1 = combine everything into this part  
   - 2 = take a recording session for this part  
   - 3 = use digital inputs for this part
> [!NOTE]
> Parts can overlap, for example, one part might include all data, while another part includes a single recording period. This allows different subsets of data to be analyses in one go.

3. If method is 2 or 3 we also need to specify which **intervals** to include. Intervals can be recording sessions (if method is set to 2) or could be a sequence of time periods that represent trials (if method is set to 3):  
   - if method = 1 value is ignored for this part  
   - if method = 2 specify which recording session(s) to include in this part, i.e. [2] would mean just the second recording, [2 3] would mean the second and third recordings.  
   - if method = 3 specify which digital input pairs to include (inf = use all) in this part, i.e. [1 3 4 5] would mean to use the first, third, fourth and fifth keypress pairs.
> [!NOTE]
> Recordings are ordered according to their order in the merged kwiktint file, which is dictated by the alphabetical order of the original files. i.e. If files named session1, session2 and session3 are combined by kwiktint, session1 will be recording 1, session2 will be recording 2 etc. Digital input pairs are ordered by the onset time of the first keypress.

5. Lastly, if the method is set to 3 (digital inputs), the intervals list tells us which digital input key pairs to use for the start/end of the trials we want to include, but we still need to specify what those digital input key pairs are, some people use 's' and 'e' to mark the start and end of each trial for example. If we are not using method 3, this value is ignored.
> [!NOTE]
> Keypresses used to delineate trial starts and ends {start,end} can be integers or characters i.e. {1,2} or {'s','e'}, {'1','2'} is the same as {1,2}. These are usually saved during recording using dacqUSB in .inp files.

For example, here we create a part named 'session1', it will include data from a single recording (method 2) which is the first recording made (interval 1) and we don't need any interval keys:  
```
part_name                               = 'session1';
part_config.(part_name).method          = 2;
part_config.(part_name).intervals       = [1];
part_config.(part_name).interval_keys   = {};
```

Here we create a part named 'session2', it will combine data from multiple recordings (method 2) which are the second and third recordings made (intervals 2 & 3) and we don't need any interval keys:  
```
part_name                               = 'session1';
part_config.(part_name).method          = 2;
part_config.(part_name).intervals       = [2 3];
part_config.(part_name).interval_keys   = {};
```

Here we create a part named 'session3', it will take data from any recording session that fell between keypresses (method 3). The interval keys it will use to specify the start and end of each trial are 's' for starts and 'e' for ends, it will combine the trials 1,2,3,4,5,6,10,11 and 15:  
```
part_name                               = 'session3';
part_config.(part_name).method          = 3;
part_config.(part_name).intervals       = [1 2 3 4 5 6 10 11 15];
part_config.(part_name).interval_keys   = {'s','e'};
```

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

# Examples
run function using default values
```
klustest()
```
run function using default values, but only on tetrodes 1 and 5
```
klustest('tetrodes',[1 5])
```

# Parameters
'tetrodes'      -   Numeric vector that specified the tetrodes to run on (i.e. [1 2 3 4 5 6] would run on tetrodes 1 to 6). Default value is 1:16

'clusters'      -   Numeric vector that specifies the clusters to run on. Set to 0 to run on all clusters. Default value is 0

'rname'         -   String or charcter vector specifying the rat name This is used in the sdata structure so its important to ensure this is correct. Default value is the name of the parent directory

'outname'       -   String or charcter vector specifying the output filename used by klustest - i.e. figures will be saved in a folder named this. Default value is 'klustest'.

'cname'         -   String or charcter vector specifying the output filename used by kwiktint - or the name used for klustakwiked files. Default value is 'kwiktint'.

> [!NOTE]
> When merging and splitting sessions, additional information needs to be provided within the klustest code, see [Specifying data 'parts'](https://github.com/Neuroesc/klustest/edit/main/README.md#specifying-data-parts).

# Figure outputs
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

Example of some output figures and a visual guide:
<img width="3627" height="2107" alt="Picture1" src="https://github.com/user-attachments/assets/e8478c3a-23c5-4caa-a2d5-72b273505dc3" />

# Data outputs
In addition to saving figure outputs, this function outputs spikes, position data and statistics in Matlab table format. This table contains different data types, some are integers while others are matrices or vectors.

Each row of this table contains:
- rat - String, the rat name, provided as an input ('rname') or assumed to be the last part of the current working directory  
- date - String, the date of recording, extracted from recording files if possible  
- partn - Scalar, the part number, given the parts specified in order at the beginning of the file, which part is this?  
- dir - String, directory of the data  
- uci - String, unique cell identifier: rat name, recording date, tetrode, cluster  
- tetrode - Scalar, tetrode this data was recorded on  
- electrode - Scalar, tetrode this data was recorded on   
- cluster - Scalar, cluster assignment in the .clu file  
- spt_pot_index - Nx1 numeric vector, where N = number of spikes, index into the position data, nearest position data point for each spike  
- spike_times - Nx1 numeric vector, where N = number of spikes, time (s) of each spike since the start of recording  
- nspikes - Scalar, total number of spikes  
- frate - Scalar, firing rate (Hz)  
- isod - 1x2, Cluster isolation distance and L-ratio (see [here](https://doi.org/10.1016/j.neuroscience.2004.09.066))  
- wave_width - Scalar, width of waveform, peak to trough of mean waveform (of the channel with the highest amplitude) (see [here](https://doi.org/10.7554/eLife.15986))  
- wave_amps - Scalar, waveform amplitude, peak of mean waveform (of the channel with the highest amplitude)  
- wave_mean - NxM, N = number of recording channels, M = spike samples, mean waveform for each channel (volts)  
- wave_std - NxM, N = number of recording channels, M = spike samples, standard deviation waveform for each channel (volts)  
- ratemap - NxM, firing rate map (see [here](https://github.com/Neuroesc/rate_mapper.git), [here](https://doi.org/10.1371/journal.pcbi.1011763))  
- spikemap - NxM, firing rate map (see [here](https://github.com/Neuroesc/rate_mapper.git), [here](https://doi.org/10.1371/journal.pcbi.1011763))    
- spatial_info - 1x4, spatial information (b/s), spatial information (b/sp), sparsity, coherence, calculated using firing rate map (see [here](https://doi.org/10.1002%2F%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K), [here](https://doi.org/10.1002%2F%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K), [here](https://dx.doi.org/10.1523%2FJNEUROSCI.1704-07.2007))  
- nfields - Scalar, number of place fields (regions greater than a set area and with a minimum peak firing rate)  
- field_data - Nx9 cell array, N = number of fields, for each field: Area, Centroid, WeightedCentroid, MajorAxisLength, MinorAxisLength, Orientation, ConvexHull, PixelIdxList, MaxIntensity  
- gridmap - NxM, spatial autocorrelation, calculated using firing rate map.  
- grid_info - 1x3, grid score, grid wavelength, grid orientation, calculated using gridmap (see [here](https://doi.org/10.1126/science.1188210))  
- grid_field_info - 1x6, field radius, major axis length, minimum axis length, height, width, orientation, calculated using the central peak of gridmap (see [here](https://doi.org/10.1073/pnas.1811867116))  
- hd_ratemap - Nx1, where N = number of HD bins, head direction firing rate map  
- hd_info - 1x4, Rayleigh vector length, preferred firing direction, mean angle, SD angle, calculated using hd_ratemap (see [here](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics))  
- hd_ratemap_half - Nx1 x two columns, where N = number of HD bins, head direction firing rate map for each half of the session (median split)  
- hd_info_half - 1x8, Rayleigh vector length, preferred firing direction, mean angle, SD angle repeated twice - once for each session half in hd_ratemap_half (see [here](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics))  
- spatial_info_z - 1x3, spatial information (b/s), grid score and Rayleigh vector length, z-scored relative to a shuffle using the Savelli method (see [here](https://doi.org/10.1002%2F%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K), [here](https://elifesciences.org/articles/21354))  
- spatial_info_p - 1x3, spatial information (b/s), grid score and Rayleigh vector length, p-values calculated using a shuffle using the Savelli method (see [here](https://doi.org/10.1002%2F%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K), [here](https://elifesciences.org/articles/21354))  
- ahv_curve - Nx1, where N = number of AHV bins, angular head velocity x firing rate map  
- theta_phase - 1xN, where N = number of theta phase bins, spike histogram relative to theta phase  
- theta_info - 1x3, Rayleigh vector length, preferred phase (bin with max count), mean angle (average phase) of spike theta phases (see [here](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics))  
- isi - 1xN, where N = number of interspike interval bins, histogram of interspike intervals  
- isi_info - not used  
- burst_index - proportion of ISIs less than 6ms (see [here](https://doi.org/10.1002/hipo.22002))  
- autocorr_theta - Nx1, where N = number of autocorrelogram bins, large window autocorrelation (usually 500ms, 10ms bins)  
- autocorr_refrac - Nx1, where N = number of autocorrelogram bins, small window autocorrelation (usually 50ms, 1ms bins)  
- autocorr_refrac_info - 1x2, number of refractory period violations (ISIs less than 2ms), and the proportion of spikes with a refractory period less than 2ms (see [here](https://doi.org/10.1152/jn.00699.2015))  
- speed_slope - 1xN, where N = number of speed bins, running speed x firing rate map (see [here](https://doi.org/10.1038/nature14622))  
- speed_info - 1x3, speed score, slope of linear fit, y-intercept of linear fit (see [here](https://doi.org/10.1038/nature14622))  
