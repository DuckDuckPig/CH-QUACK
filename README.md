# CH_QUACK <br> Quantifying Unipolarity via Active Contour Kenetics
This code develops segmentations of coronal holes (CHs) utilizing solar extreme ultraviolet (EUV) images and magnetogram data taken from the Atmospheric Imager Assembly (AIA) and Heliosismic and Magnetic Imager (HMI) (respectivly) aboard the Solar Dynamics Observatory (SDO). This code is an extension of the Active Contours Without Edges for Coronal Hole (CH-ACWE) segmentation method (the original public release for which can be found at [https://github.com/DuckDuckPig/CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE)) that adds additional forces to quantify the unipolarity of the underlying magnetic field in identified regions under the expectation that CHs, being regions of open magnetic field, will appear unipolar.

## Requirements: [environment.yml](environment.yml)
This updated environment file specifies the packages necessary to run this code. It is backwards compatible with [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE).

## General Notes About This repository
- Throughout all code/scrips the `Key Variables` cell (usually `In[2]`) will need to be updated to point to the correct directories.
- The folder `FinalPipeline` is a copy of this repository that only contains the tools for developing a segmentation according to the final EUV+HMI magnetogram segmentation process. Users interested in only the final segmentation process will find a self-contained version of the code, with updated functions to facilitate easier segmentations there.

## Downloading the dataset
The dataset used for this project is identical to the one from [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE). The code from that repository responsible for the download process is reproduced in the the `DatasetTools` folder. The instructions from the [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE) `README` file are reproduced below:

> - `DownloadLists`: This folder contains an organized lists of the dataset. These lists are organized into four `.csv` files, one for each Carrington Rotation (CR).
> - `Carrington Rotation Start Dates.csv`: This file is a list of the start dates for each Carrington Rotation from CR -10 through CR 2300. This file is used by `DownloadByRotation.py` for both downloading and organizing the dataset.
> - `DataManagmentTools.py`: Tools/functions for formatting dates to allow for request of data from [jsoc.stanford.edu](jsoc.stanford.edu).
> - `RebuildDataset.py`: Find and download any file that is missing from the dataset folder.
>   - This script will rebuild the dataset directly from the specified file present in the `DownloadLists` folder.
>   - User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
>   - User will need to register their email at [http://jsoc.stanford.edu/ajax/register_email.html](http://jsoc.stanford.edu/ajax/register_email.html) and add that email to the appropriate variable in the `Key Variables` cell.
> - `DownloadByRotation.py`: Download `aia.lev1_euv_12s` and `hmi.M_720s` images at a 1 hour cadence for the specified Carrington rotation(s). 
>   - This script will attempt to omit any time frame wherein at least one file is missing or does not and have a `QUALITY` key of `0`. 
>   - User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories.
>   - User will need to register their email at [http://jsoc.stanford.edu/ajax/register_email.html](http://jsoc.stanford.edu/ajax/register_email.html) and add that email to the appropriate variable in the `Key Variables` cell.
>   - This script can be used to speed up the process of rebuilding the dataset. This is achieved by
>     1. Creating a temporary subfolder within the `DatasetTools` folder
>     2. Ensuring that the variable `traceFolder` points to that temporary subfolder
>     3. Setting the remaining variables in the `Key Variables` cell to ensure the correct CR is downloaded and saved where the user wishes
>     4. Running `DownloadByRotation.py`
>     5. Deleting the temporary subfolder
>     6. Running `RebuildDataset.py` with `traceFolder = 'DownloadLists/'` to download any missing files
> - `GapCheck.py`: Inform the user as to the largest hour gap between entries in the specified CR within the dataset.

## EUV Single Channel Segmentations
The folder `EUV_SingleChannel` contains code and tools for comparing the new segmentations and the prior segmentation. Within this folder you will find:

### `EUV_SingleChannel/Analysis`
This folder contains tools for quantifying various aspects of the segmented regions, from initial seed through final, fully evolved, segmentation, for both EUV-only segmentations and EUV+HMI magnetogram segmentations.
- `Expanded Uni on Uni.py`: Calculate the unipolarity of each region from initial seed, through each iteration, until final segmentation for the EUV+HMI magnetogram segmentations.
  - Requires the output of `UniCompair Seed to Segmentation.ipynb`
  - Requires the segmentations from `HMI_Experiments/TestEvolutionMethods/Standard/runACWEdefaltSampleUnipolarity_history.py`
  - The `Key Variables` cell is labeled `In[1b]`
- `Expanded Uni`: Calculate the unipolarity of each region from initial seed, through each iteration, until final segmentation for the EUV-only segmentations.
  - Requires the output of `UniCompair Seed to Segmentation.ipynb`
  - Requires the segmentations from `EUV_SingleChannel/Standard/runACWEdefault_history.py`
  - The `Key Variables` cell is labeled `In[1b]`
- `Find Samples.ipynb`: Plot the unipolarity of individual regions from initial seed to final segmetation as a function of iteration
  - Requires the output of `Expanded Uni`
  - Requires the output of `Expanded Uni on Uni.py`
- `MapUni.py`: Calculate the unipolarity of each region in the final segmentation
- `SeedAndSegCounts.py`: Count the number of regions in the initial seed and final segmentation for the EUV-only segmentations
- `SeedMapUni.py`: Calculate the unipolarity of the each region in the initial seed
- `UniCompair Seed to Segmentation.ipynb`: Plot unipolarity of initial seed vs unipolarity of final segmentation for the specified CR
  - Note: Requires the output of both `MapUni.py` and `SeedMapUni.py` for both full-scale and one-eighth-scale segmentations

### `EUV_SingleChannel/ConfidenceMapping`
This folder contains a copy of the script `runACWEconfidenceLevelSet_Default.py` from [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE). The instructions from the [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE) `README` file are reproduced below:
> - User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories and desired EUV wavelength (193 angstroms is the assumed default).
> - The script will assume that the data are organized by CR, with a sub directory for each record time in the `.csv` file in the `DownloadLists` subfolder within the `DatasetTools` directory. Both `DownloadByRotation.py` and `RebuildDataset.py` will organize the dataset appropriately.

### `EUV_SingleChannel/Scaled`
This folder contains a copy of the script `runACWEscaledDefault.py` from [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE). The instructions from the [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE) `README` file are reproduced below:
> - [EUV-only] Segmentations generated at any spatial resolution other than 512x512 pixels should be performed using the script `runACWEscaledDefault.py`
>   - User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories and desired EUV wavelength (193 angstroms is the assumed default).
>   - The script will assume that the data are organized by CR, with a sub directory for each record time in the `.csv` file in the `DownloadLists` subfolder within the `DatasetTools` directory. Both `DownloadByRotation.py` and `RebuildDataset.py` will organize the dataset appropriately.
>   - Under the assumption that the input image is 4096x4096 pixels (the resolution of AIA), if the resize parameter variable `resize_param = 8`, this function will generate a standard [EUV-only] segmentation.

### `EUV_SingleChannel/Standard`
This folder contains two copies of the standard EUV-only segmentation script from [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE). The script `runACWEdefault.py` is the unaltered script. The script `runACWEdefault_history.py` is a modified version that saves the output of every single ACWE iteration to facilitate the analysis in `EUV_SingleChannel/Analysis/Expanded Uni.py`. The instructions from the [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE) `README` file are reproduced below:
> - User will need to adjust the variables in the `Key Variables` cell (`In[2]`) to point to the correct directories and desired EUV wavelength (193 angstroms is the assumed default).
> - The script will assume that the data are organized by CR, with a sub directory for each record time in the `.csv` file in the `DownloadLists` subfolder within the `DatasetTools` directory. Both `DownloadByRotation.py` and `RebuildDataset.py` will organize the dataset appropriately.
