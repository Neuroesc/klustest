# klustest
Klustest, cluster analysis tool in Matlab.

Following automated or manual spike sorting using Klustakwik (see [kwiktint](https://github.com/Neuroesc/kwiktint)), manual curation and refinement (see Tint), this Matlab function loads the spike sorted data, performs analyses (see below), creates a summary figure and saves an output in Matlab table format.



prepares command line instructions and files necessary to automatically run Klustakwik (via Tint) on Axona recording files. The function can automatically move from one tetrode to the next and spike sort multiple tetrodes simultaneously, saving the user a lot of time. This function also allows merging separate sessions into a single, combined output. Because merging is performed by Tint there is no risk or errors during concatenation. Merging is preferable when neurons are recorded across multiple sessions as all of the sessions will be spike sorted together, ensuring consistency. Any klustakwik settings that can be personalised in Tint can be specified in kwiktint. Additionally, kwiktint will automatically exclude grounded/LFP recording channels from spike sorting.

<img width="2417" height="778" alt="image" src="https://github.com/user-attachments/assets/833dcdc8-e484-451c-8acd-d21d17621c4a" />

> [!NOTE]
> This function is intended to run after spike sorting using [kwiktint](https://github.com/Neuroesc/kwiktint).

# Outputs
This function performs spike-time analyses and creates spike-time visualisation plots:
- Plots 50 ms inter-spike-interal
  - calculates ISI full width to half maximum
  - calculates burst index
- Plots 50 ms higher order inter-spike-interval (i.e. spike autocorrelograms)
  - calculates proportional refractory period violations
- Plots 500 ms higher order inter-spike-interval (i.e. spike autocorrelograms)
  - calculates theta index
  - calculates theta frequency
- Plots theta-phase spike relationship
  - calculates calculates Rayleigh vector length (theta modulation)
  - calculates preferred theta phase


# Usage
process files with default settings:
```
kwiktint()
```

'Screening mode'  - If Matlab is in a directory containing only one session (1 .set file) and no inputs are given this function will automatically analyse that session and look for all tetrodes. Output files are named 'kwiktint'.  

'Experiment mode' - If multiple .set files are present the function will ask the user to identify which sessions to analyse. It will then assume these sessions should be kkwiked together and that all tetrodes should be analysed. Output files are named 'kwiktint'.  

'Batch mode'      - If you wish to kkwik multiple sessions, but keep the outputs seperate and in different .cut files then run the function with the 'combine' input set to 0. The function will continue to ask which sessions should be analysed, but these will be kkwiked seperately. Output files are named after each original session name, to avoid overwriting anything  

process with Name-Value pairs used to control aspects of the cluster cutting:
```
kwiktint(Name,Value,...)
```
<img width="724" height="596" alt="image" src="https://github.com/user-attachments/assets/415573c9-7ab5-43e5-bd4e-fc1ea3e6a5b4" />

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

<img width="1613" height="743" alt="image" src="https://github.com/user-attachments/assets/ee3dae1b-1c15-4bd6-9102-59676392968d" />

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
