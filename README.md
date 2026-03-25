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
- [Additional code](https://github.com/Neuroesc/klustest/edit/main/README.md#additional-code)

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
## sdata
In addition to saving figure outputs, this function outputs spikes and statistics in Matlab table format (sdata.mat). This table contains different data types, some are integers while others are matrices or vectors. This table can be loaded using:  
```
load('sdata.mat','-mat')
```  
Each row of this table contains:
- sdata.rat - String, the rat name, provided as an input ('rname') or assumed to be the last part of the current working directory  
- sdata.date - String, the date of recording, extracted from recording files if possible  
- sdata.partn - Scalar, the part number, given the parts specified in order at the beginning of the file, which part is this?  
- sdata.dir - String, directory of the data  
- sdata.uci - String, unique cell identifier: rat name, recording date, tetrode, cluster  
- sdata.tetrode - Scalar, tetrode this data was recorded on  
- sdata.electrode - Scalar, tetrode this data was recorded on   
- sdata.cluster - Scalar, cluster assignment in the .clu file  
- sdata.spt_pot_index - Nx1 numeric vector, where N = number of spikes, index into the position data, nearest position data point for each spike  
- sdata.spike_times - Nx1 numeric vector, where N = number of spikes, time (s) of each spike since the start of recording  
- sdata.nspikes - Scalar, total number of spikes  
- sdata.frate - Scalar, firing rate (Hz)  
- sdata.isod - 1x2, Cluster isolation distance and L-ratio (see [here](https://doi.org/10.1016/j.neuroscience.2004.09.066))  
- sdata.wave_width - Scalar, width of waveform, peak to trough of mean waveform (of the channel with the highest amplitude) (see [here](https://doi.org/10.7554/eLife.15986))  
- sdata.wave_amps - Scalar, waveform amplitude, peak of mean waveform (of the channel with the highest amplitude)  
- sdata.wave_mean - NxM, N = number of recording channels, M = spike samples, mean waveform for each channel (volts)  
- sdata.wave_std - NxM, N = number of recording channels, M = spike samples, standard deviation waveform for each channel (volts)  
- sdata.ratemap - NxM, firing rate map (see [here](https://github.com/Neuroesc/rate_mapper.git), [here](https://doi.org/10.1371/journal.pcbi.1011763))  
- sdata.spikemap - NxM, firing rate map (see [here](https://github.com/Neuroesc/rate_mapper.git), [here](https://doi.org/10.1371/journal.pcbi.1011763))    
- sdata.spatial_info - 1x4, spatial information (b/s), spatial information (b/sp), sparsity, coherence, calculated using firing rate map (see [here](https://doi.org/10.1002%2F%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K), [here](https://doi.org/10.1002%2F%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K), [here](https://dx.doi.org/10.1523%2FJNEUROSCI.1704-07.2007))  
- sdata.nfields - Scalar, number of place fields (regions greater than a set area and with a minimum peak firing rate)  
- sdata.field_data - Nx9 cell array, N = number of fields, for each field: Area, Centroid, WeightedCentroid, MajorAxisLength, MinorAxisLength, Orientation, ConvexHull, PixelIdxList, MaxIntensity  
- sdata.gridmap - NxM, spatial autocorrelation, calculated using firing rate map.  
- sdata.grid_info - 1x3, grid score, grid wavelength, grid orientation, calculated using gridmap (see [here](https://doi.org/10.1126/science.1188210))  
- sdata.grid_field_info - 1x6, field radius, major axis length, minimum axis length, height, width, orientation, calculated using the central peak of gridmap (see [here](https://doi.org/10.1073/pnas.1811867116))  
- sdata.hd_ratemap - Nx1, where N = number of HD bins, head direction firing rate map  
- sdata.hd_info - 1x4, Rayleigh vector length, preferred firing direction, mean angle, SD angle, calculated using hd_ratemap (see [here](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics))  
- sdata.hd_ratemap_half - Nx1 x two columns, where N = number of HD bins, head direction firing rate map for each half of the session (median split)  
- sdata.hd_info_half - 1x8, Rayleigh vector length, preferred firing direction, mean angle, SD angle repeated twice - once for each session half in hd_ratemap_half (see [here](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics))  
- sdata.spatial_info_z - 1x3, spatial information (b/s), grid score and Rayleigh vector length, z-scored relative to a shuffle using the Savelli method (see [here](https://doi.org/10.1002%2F%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K), [here](https://elifesciences.org/articles/21354))  
- sdata.spatial_info_p - 1x3, spatial information (b/s), grid score and Rayleigh vector length, p-values calculated using a shuffle using the Savelli method (see [here](https://doi.org/10.1002%2F%28SICI%291098-1063%281996%296%3A2%3C149%3A%3AAID-HIPO6%3E3.0.CO%3B2-K), [here](https://elifesciences.org/articles/21354))  
- sdata.ahv_curve - Nx1, where N = number of AHV bins, angular head velocity x firing rate map  
- sdata.theta_phase - 1xN, where N = number of theta phase bins, spike histogram relative to theta phase  
- sdata.theta_info - 1x3, Rayleigh vector length, preferred phase (bin with max count), mean angle (average phase) of spike theta phases (see [here](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics))  
- sdata.isi - 1xN, where N = number of interspike interval bins, histogram of interspike intervals  
- sdata.isi_info - not used  
- sdata.burst_index - proportion of ISIs less than 6ms (see [here](https://doi.org/10.1002/hipo.22002))  
- sdata.autocorr_theta - Nx1, where N = number of autocorrelogram bins, large window autocorrelation (usually 500ms, 10ms bins)  
- sdata.autocorr_refrac - Nx1, where N = number of autocorrelogram bins, small window autocorrelation (usually 50ms, 1ms bins)  
- sdata.autocorr_refrac_info - 1x2, number of refractory period violations (ISIs less than 2ms), and the proportion of spikes with a refractory period less than 2ms (see [here](https://doi.org/10.1152/jn.00699.2015))  
- sdata.speed_slope - 1xN, where N = number of speed bins, running speed x firing rate map (see [here](https://doi.org/10.1038/nature14622))  
- sdata.speed_info - 1x3, speed score, slope of linear fit, y-intercept of linear fit (see [here](https://doi.org/10.1038/nature14622))  

## bdata
Alongside spikes and firing rate statistics, this function also stores position data and other session level data in Matlab table format (bdata.mat). This table contains different data types, some are integers while others are matrices or vectors. For convenience, this table is stored in the custom property field of sdata:
```
bdata = sdata.Properties.CustomProperties.bdata;
```
Each row of this table contains:
- bdata.rat - String, the rat name, provided as an input ('rname') or assumed to be the last part of the current working directory, mainly included so data can be aligned with sdata  
- bdata.date - String, the date of recording, extracted from recording files if possible, mainly included so data can be aligned with sdata  
- bdata.partn - Scalar, the part number, given the parts specified in order at the beginning of the file, which part is this?, mainly included so data can be aligned with sdata  
- bdata.dir - String, directory of the data, mainly included so data can be aligned with sdata  
- bdata.pos - Nx7, where N = the number of position data points: time (s), x (cm), y (cm), head direction (deg), speed (cm/s), angular head velocity (deg/s), displacement (cm/s)
- bdata.dwellmap - NxM, firing rate map (see [here](https://github.com/Neuroesc/rate_mapper.git), [here](https://doi.org/10.1371/journal.pcbi.1011763))  
- bdata.hd_dwellmap - Nx1, where N = number of HD bins, head direction dwell time map (s)  
- bdata.hd_dwellmap_half - Nx1 x two columns, where N = number of HD bins, head direction dwell time map for each half of the session (median split, s)  
- bdata.ahv_dwellmap - Nx1, where N = number of AHV bins, angular head velocity x dwell time map (s)  
- bdata.speed_dwell_time - 1xN, where N = number of speed bins, running speed x dwell time map (see [here](https://doi.org/10.1038/nature14622))  

## pdata
Additionally, klustest also stores session specific settings and data in Matlab structure format (pdata.mat). This struct is stored in the custom property field of sdata:  
```
pdata = sdata.Properties.CustomProperties.pdata;
```
This struct contains fields such as:  
- pdata.rat - String, the rat name, provided as an input ('rname') or assumed to be the last part of the current working directory  
- pdata.date - String, the date of recording, extracted from recording files if possible  
- pdata.analysed - String, the date the data were last analysed using klustest  
- pdata.directory - String, directory of the data  
- pdata.tetrodes - Nx2, where N = number of electrodes, column 1 has the electrode number, column 2 contains the electrode type, 0 = tetrode, 1 = stereotrode  
- pdata.sessions - Scalar, the number of recording sessions contributing to the merged file  
- pdata.data_dirs - Nx1 cell array, directories of each of the recording sessions  
- pdata.snames - Nx1 cell array, the names of each of the recording sessions (filename format)  
- pdata.mapset - Structure, settings used to create firing rate maps and treat position data. Not all of these settings may be used. These parameters can be changed at the top of klustest and are stored in pdata for reference only.  
   - pdata.ppm - pixels per meter    
   - pdata.jumpcut - z-score to use when rejection jumps in position data  
   - pdata.jumpwindow - time window over which to compute jumps    
   - pdata.method - method for creating firing rate maps    
   - pdata.binsize - bin size for firing rate maps    
   - pdata.ssigma - smoothing strength for firing rate maps    
   - pdata.padding - padding to use for firing rate maps    
   - pdata.mindwell - minimum dwell time for a bin to be considered visited in firing rate maps    
   - pdata.mindist - minimum distance a bin must be from position data to be considered visited in firing rate maps      
   - pdata.smethod - smoothing method for firing rate maps    
   - pdata.zcut - z-score cutoff to use when thresholding firing rate maps to detect place fields  
   - pdata.frcut - minimum peak firing rate cutoff to use when excluding place fields  
   - pdata.arcut - minimum field area to use when excluding place fields  
   - pdata.fix_aspect - 1 means firing rate map etc are oriented with the long axis horizontal  
   - pdata.wave_window - time window over which to show waveforms    
   - pdata.hd_displace - if set to 1, will use displacement (movement direction) instead of head direction   
   - pdata.hd_type - 'density' corresponds to a circular kolmogorov smirnov density estimate plot, 'histogram' corresponds to a traditional histogram polar plot  
   - pdata.hd_bins - the number of bins to use when computing HD plot  
   - pdata.hd_sigma - the standard deviation of the gaussian kernel used in HD circular density estimate  
   - pdata.hd_boxcar - the number of bins over which to compute the HD histogram boxcar average  
- pdata.cname - Name used for the merged cluster cutting output, usually kwiktint  
- pdata.outname - Name to use for outputs and output directories, usually klustest  
- pdata.pos - Nx7 table, position data for entire session: time (s), x (cm), y (cm), head direction (deg), speed (cm/s), angular head velocity (deg/s), displacement (cm/s)  
- pdata.pos_srate - Sample rate of the position data (Hz)  
- pdata.tstart - Start time of the entire recording, Neuralynx files have a a start time equal to the actual time, we need to zero this, this value retains the original start time  
- pdata.clusters - Nx1 cell array, for each electrode, the cluster assignment of every spike  
- pdata.spike_times - Nx1 cell array, for each electrode, the time of every spike (s)  
- pdata.wavtime - 1xN, where N = spike samples, the time relative to the waveform peak, for plotting waveforms (i.e. the x-axis values for waveforms stored in sdata.wave_mean)  
- pdata.part_config - Structure, copy of the part_config structure which is also saved as a .json file  
- pdata.ahv_xvalues - 1xN, where N = number of AHV bins, bin positions for AHV (i.e. the x-axis values for AHV ratemaps in sdata.ahv_curve)  
- pdata.isi_xvalues - 1xN, where N = number of interspike interval bins, bin positions for ISI histogram (i.e. the x-axis values for ISI histograms in sdata.isi)  
- pdata.autocorr_theta_xvalues - 1xN, where N = number of autocorrelogram bins, bin positions for large window autocorrelogram (i.e. the x-axis values for autocorrelograms in sdata.autocorr_theta)
- pdata.autocorr_theta_evalues - 1xN, where N = number of autocorrelogram bins **+ 1**, bin **edges** for large window autocorrelogram  
- pdata.autocorr_refrac_xvalues - 1xN, where N = number of autocorrelogram bins, bin positions for small window autocorrelogram (i.e. the x-axis values for autocorrelograms in sdata.autocorr_refrac)  
- pdata.autocorr_refrac_evalues - 1xN, where N = number of autocorrelogram bins **+ 1**, bin **edges** for small window autocorrelogram

# Additional code
Some additional code is provided in the 'additional' directory which can be used to aid analysis.

## get_combined_sdata.m 
Can be used to combine multiple sdata tables (the output of Klustest) into one merged dataset. This merged dataset can be useful when analysing the overall results of a project. However it can also be used to run any general functions on every recording session in a dataset. This can be very useful if you would like to rerun kwiktint or klustest on your entire dataset.
  - This function expects data to be arranged in a specific file format:  
``project > experiment > rat > date``
  - 'project' is the overall project you want to combine the data from, 'experiments' are the sub experiments conducted as part of this project that you want to combine, 'rats' are the animals you want to be combined, and 'dates' are recording sessions which are each saved in their own directory. I would recommend using numeric date formats for the 'dates' directories, formatted YYYYMMDD.  
  - The input data_dir to get_combined_sdata (which defaults to the current working directory) should be the topmost project directory.  
> [!NOTE]
> get_combined_sdata will only work on data that has already been analysed using klustest, as it can only combine existing sdata tables

## add_manual_cell_type.m 
Can be used to automatically and manually add cell type identifications to all clusters in an sdata. Combined with get_combined_sdata, this code can be used to manually identify cells throughout an entire dataset, allowing easy filtering later by curated cell type.  
- This function is intended to be used in combination with get_combined_sdata, but you can also run it on any single recording directory - or direct it to a recording directory using the data_dir Name-Value input
> [!NOTE]
> add_manual_cell_type will only work on data that has already been analysed using klustest, as it can only add info to an existing sdata table

<img width="1808" height="767" alt="image" src="https://github.com/user-attachments/assets/fdf7da4b-9cf0-4cf6-863e-d0da4ac68b2f" />

When run, add_manual_cell_type moved through each cluster and presents the user with an interactive plot. This shows the spike plot, firing rate map, spatial autocorrelogram, head direction tuning curve, waveform, refractory period spike autocorrelogram and theta modulation spike autocorrelogram (see figure). There are also buttons providing the main cell type options that can be used as categories. The automated cell identification is shown by the green highlighted button.  
  - Both the automated cell type and manually curated cell type are saved into a new column of sdata  
  - Because this information is added to a new column, sdatas with no cell typing and sdatas with cell typing cannot be merged unless an empty cell type column is added to the former first  


