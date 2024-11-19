# CH_QUACK <br> Quantifying Unipolarity via Active Contour Kenetics
This code develops segmentations of coronal holes (CHs) utilizing solar extreme ultraviolet (EUV) images and magnetogram data taken from the Atmospheric Imager Assembly (AIA) and Heliosismic and Magnetic Imager (HMI) (respectively) aboard the Solar Dynamics Observatory (SDO). This code is an extension of the Active Contours Without Edges for Coronal Hole (CH-ACWE) segmentation method (the original public release for which can be found at [https://github.com/DuckDuckPig/CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE)) that adds additional forces to quantify the unipolarity of the underlying magnetic field in identified regions under the expectation that CHs, being regions of open magnetic field, will appear unipolar.

## Requirements: [environment.yml](environment.yml)
This updated environment file specifies the packages necessary to run this code. It is backwards compatible with [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE).

## General Notes About This repository
- In every python (`.py`) or jupyter notebook (`ipynb`) file the `Key Variables` cell (usually `In[2]`) will need to be updated to point to the correct directories, unless noted below:
  - None of the tools in the `ACWE_python_fall_2023` folder need to be updated
  - None of the tools in `Metrics` folder need to be updated
  - The code `DatasetTools/DataManagmentTools.py` does not need to be updated
  - None of the tools in `SDO_tools` folder need to be updated
- The version of this code that you are seeing here is taken from the folder `FinalPipeline`. This version of the repository only contains the tools needed to develop a segmentation according to the final EUV+HMI magnetogram segmentation process. Users interested in all code related to CH_QUACK may visit [https://github.com/DuckDuckPig/CH-QUACK](https://github.com/DuckDuckPig/CH-QUACK) for the complete repository.

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

## General Tools
The folder `ACWE_python_fall_2023` contains updated versions of the functions used to generate segmentations, both with and without HMI magnetogram data, and for saving the resulting segmentations. The following scripts are identical between this code and [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE), as such the instructions from the [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE) `README` file are reproduced below:

> - `acweConfidenceMapTools_v3.py`: Tools/functions for combining a segmentation group (collection of segmentations from the same EUV observation) in order to generate a confidence map.
> - `acweFunctions_v6.py`: Tools/functions for preprocessing an EUV image, generating an initial mask, and running ACWE for both single output/segmentation and for a confidence map. 
>   - The function `run_acwe` performs all processing and returns the final segmentation and initial mask. 
>   - The function `run_acwe_confidenceMap` performs all processing and returns the final confidence map as a series of segmentations and initial mask.
>   - Additional functions are also provided to perform each step separately.
>   - These functions will work for both AIA and Solar Terrestrial RElations Observatory (STEREO) observations, however a resize parameter of 4 and seeding parameter `alpha` in the range \[0.8,0.9\] are recommended for STEREO data. 
> - `acweRestoreScale.py`: Tools/functions for resizing a segmentation to match the spatial resolution of the input image.
>   - Upscale a confidence map using `upscaleConMap`
>   - Upscale a single segmentation using `upscale` 
>   - Both functions take in the ACWE header and the segmentation or confidence map and return the same segmentation or confidence map, upscaled to match the resolution of the original EUV image.
> - `acweSaveSeg_v5.py`: Tools/functions for saving and opening segmentations. 
>   - The function `saveSeg` takes in the header of the original EUV image, the final segmentation(s), and the list of ACWE parameters. It generates an .npz file which saves the final segmentation with a header outlining the ACWE parameters and a copy of the header for the original EUV image. 
>   - The function `openSeg` opens and returns the header of the original EUV image, as a dictionary, the header outlining the options used to generate the ACWE segmentation, organized as a dictionary, and the final ACWE segmentation(s).
>   - Both functions work for both single segmentations and for confidence maps.

The folder `ACWE_python_fall_2023` also contains:
- `acweUniFunctions_v1.py`: A version of `acweFunctions_v6.py` that accepts vector-valued data and also allows the user to evolve the contour to maximize the unipolarity of the underlying region, as observed in the magnetic field.

In addition to this, the folder `ACWE_python_fall_2023/ACWE_python_v3`, which still contains the original ACWE functions, modified to work on python version 3.0 or higher, also includes the following new scripts/versions of ACWE:
- `acwe_exp1.py`: ACWE for vector-valued images
- `acwe_exp2.py`: ***The final CH-QUACK evolution function.*** This version has an additional pair of forces that evolve the contour to maximize the unipolarity of the underlying region, as observed in the magnetic field
- `acwe_exp3.py`: This version has an additional pair of forces that evolve the contour to maximize the flux imbalance of the underlying region, as observed in the magnetic field.
- `acwe_exp4*.py`: This version has an additional pair of forces that evolve the contour to maximize the absolute skew of the underlying region, as observed in the magnetic field.

The folder `SDO_tools` contains the file `read_fits.py`. This file provides tools for opening and alighting SDO-AIA and SDO-HMI data.

The folder `Metrics` contains tools for comparing segmentations to each other.

## Generating a QUACK Segmentation
The file `Standard/runACWEunipolarity.py` Generates a binary CH segmentation according to the final QUACK pipeline. Please note:
- The variables in the `Key Variables` cell (`In[2]`) will need to be adjusted to point to the correct directories.
- As written, this script assume that the user will be using the 193 angstroms EUV observation and the HMI magnetogram, however this can be changed to take in any combination of EUV and HMI magnetogram data by adjusting the `acweChoices` variable to list all image(s) chosen for the process.
  - The size of the foreground, background, and `alpha` weights must be adjusted to match the number of inputs and order chosen
  - Always list the magnetogram file last
  - Setting `alpha[i] = -1` will result in the seeding function ignoring channel `i` when seeding. The final seed will be the union of all seeds produced. The channel with the HMI magnetogram should therefore be set to -1.
- As written, this code will perform 50 EUV-only iterations of ACWE to enlarge regions of the initial seed in order to filter out filament regions. The number of EUV-only iterations can be changed by adjusting the variable `switch` located in the `Key Variable` cell.
- Like the scripts in [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE): "The script will assume that the data are organized by CR, with a sub directory for each record time in the `.csv` file in the `DownloadLists` subfolder within the `DatasetTools` directory. Both `DownloadByRotation.py` and `RebuildDataset.py` will organize the dataset appropriately."
