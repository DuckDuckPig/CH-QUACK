# CH_QUACK
This code develops segmentations of coronal holes (CHs) utilizing solar extreme ultraviolet (EUV) images and magnetogram data taken from the Atmospheric Imager Assembly (AIA) and Heliosismic and Magnetic Imager (HMI) aboard the Solar Dynamics Observatory (SDO). This code is an extension of the Active Contours Without Edges for Coronal Hole (CH-ACWE) segmentation method (the original public release for which can be found at [https://github.com/DuckDuckPig/CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE)) that adds additional forces to gauge the unipolarity of the underlying magnetic field of the identified regions under the expectation that CHs, being regions of open magnetic field, will appear unipolar.

## Requirements: [environment.yml](environment.yml)
This updated environment file specifies the packages necessary to run this code. It is backwards compatible with [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE).

## General Notes About This repository
- Throughout all code/scrips the `Key Variables` cell (usually `In[2]`) will need to be updated to point to the correct directories.
- The folder `FinalPipeline` is a copy of this repository that only contains the tools for developing a segmentation according to the final EUV+HMI magnetogram segmentation process. Users interested in only the final segmentation process will find a self-contained version of the code, with updated functions to facilitate easier segmetnations there.

## Downloading the dataset
The dataset used for this project is identical to the one from [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE). The code from that repository respoinsible for the download process is reproduced in the the `DatasetTools` folder. The instructions from the [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE) `README` file are repoduced below:

> - `DownloadLists`: This folder contains an organized lists of the dataset. These lists are organized into four `.csv` files, one for each Carrington Rotation (CR).
> - `Carrington Rotation Start Dates.csv`: This file is a list of the start dates for each Carrington Rotation from CR -10 through CR 2300. This file is used by `DownloadByRotation.py` for both downloading and organizing the dataset.
> - `DataManagmentTools.py`: Tools/functions for formatting dates to allow for request of data from [jsoc.stanford.edu](jsoc.stanford.edu).
> - `RebuildDataset.py`: Find and download any file that is missing from the dataset folder.
>   - This script will rebuild the dataset direcly from the specified file present in the `DownloadLists` folder.
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
The folder `EUV_SingleChannel` contains code and tools for comparing the new segmenations and the prior segmentation
