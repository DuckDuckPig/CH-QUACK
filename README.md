# CH_QUACK <br> Quantifying Unipolarity via Active Contour Kenetics
This code develops segmentations of coronal holes (CHs) utilizing solar extreme ultraviolet (EUV) images and magnetogram data taken from the Atmospheric Imager Assembly (AIA) and Heliosismic and Magnetic Imager (HMI) (respectively) aboard the Solar Dynamics Observatory (SDO). This code is an extension of the Active Contours Without Edges for Coronal Hole (CH-ACWE) segmentation method (the original public release for which can be found at [https://github.com/DuckDuckPig/CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE)) that adds additional forces to quantify the unipolarity of the underlying magnetic field in identified regions under the expectation that CHs, being regions of open magnetic field, will appear unipolar.

## Requirements: [environment.yml](environment.yml)
This updated environment file specifies the packages necessary to run this code. It is backwards compatible with [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE).

## General Notes About This repository
- In every python (`.py`) or jupyter notebook (`ipynb`) file the `Key Variables` cell (usually `In[2]`) will need to be updated to point to the correct directories, unless noted below:
  - None of the tools in the `ACWE_python_fall_2023` folder need to be updated
  - None of the tools in `Metrics` folder need to be updated
  - The code `DatasetTools/DataManagmentTools.py` does not need to be updated
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

In addition to this, the folder `ACWE_python_v3`, which still contains the original ACWE functions, modified to work on python version 3.0 or higher, also includes the following new scripts/versions of ACWE:
- `acwe_exp1.py`: ACWE for vector-valued images
- `acwe_exp2.py`: ***The final CH-QUACK evolution function.*** This version has an additional pair of forces that evolve the contour to maximize the unipolarity of the underlying region, as observed in the magnetic field
- `acwe_exp3.py`: This version has an additional pair of forces that evolve the contour to maximize the flux imbalance of the underlying region, as observed in the magnetic field.
- `acwe_exp4.py`: This version has an additional pair of forces that evolve the contour to maximize the absolute skew of the underlying region, as observed in the magnetic field.

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

## EUV + HMI Magnetogram Segmentations
The folder `HMI_Experiments` contains code for evaluating the effects of incorporating various metrics for constraining CH evolution based on the characteristics of the underlying magnetic field both during the evolution of the region, and in the initial seeding of the CH regions.

### `HMI_Experiments/TestEvolutionMethods/Scaled`
This folder contains tools for performing full-scale EUV+HMI magnetogram segmetnations, with additional tools for comparing segmetnations to the full-scale EUV-only segmentations of [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE), and to one-eighth-scale EUV+HMI magnetogram segmentations.

#### Test Evolution with Novel Force:
- `runACWEscaledSample*.py`: These four files will generate segmentations using the specified force to guide contour evolution based on the magnetic field data. We note that unipolarity is the final force employed in this work, however, the final implementation includes a prefiltering process that is explored in the `HMI_Experiments/TestSeedingMethods` folder. Across all four scripts, the following notes should be observed:
  - The variables in the `Key Variables` cell (`In[2]`) will need to be adjusted to point to the correct directories.
  - As written, these scripts assume that the user will be using the 193 angstroms EUV observation and the HMI magnetogram, however this can be changed to take in any combination of EUV and HMI magnetogram data by adjusting the `acweChoices` variable to list all image(s) chosen for the process.
    - The size of the foreground, background, and `alpha` weights must be adjusted to match the number of inputs and order chosen
    - Always list the magnetogram file last
    - Setting `alpha[i] = -1` will result in the seeding function ignoring channel `i` when seeding. The final seed will be the union of all seeds produced. The channel with the HMI magnetogram should therefore be set to -1.
  - Like the scripts in [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE): "The script will assume that the data are organized by CR, with a sub directory for each record time in the `.csv` file in the `DownloadLists` subfolder within the `DatasetTools` directory. Both `DownloadByRotation.py` and `RebuildDataset.py` will organize the dataset appropriately."
  - Also like [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE): "Under the assumption that the input image is 4096x4096 pixels (the resolution of AIA), if the resize parameter variable `resize_param = 8`, this function will generate a standard [EUV+HMI] segmentation."

#### Analysis Compared to EUV-Only Segmentations
- `Analysis_EUV/AreaChecks/Empty Segs.ipynb`: Report the names of segmentations that are empty. User will need to update the variables in the `Key Variables` cell (`In[2]`) to match both the correct directories and the image force employed for the HMI magnetogram data.
- `Analysis_EUV/AreaChecks/Large Area Look.ipynb`: Report and display cases that have a very large area identified as pertaining to a CH. User will need to update the variables in the `Key Variables` cell (`In[2]`) to match both the correct directories and the image force employed for the HMI magnetogram data.
- `Analysis_EUV/Single CR Lowest IOU Cases-SameSet.ipynb` and `Analysis_EUV/Single CR Lowest IOU Cases.ipynb`: Display the cases where the full-scale EUV-only and full-scale EUV+HMI magnetogram segmetnations differ the most.
  - The variables in the `Key Variables` cell (`In[2]`) will need to be adjusted to point to the correct directories and chosen CR.
  - The user will also need to update the `magMethod` variable (also in the `Key Variables` cell) to select the image force employed for the HMI magnetogram data.
  - In `Analysis_EUV/Single CR Lowest IOU Cases-SameSet.ipynb`, the variable `compMethod` allows the user to view the same cases automatically chosen for a different image force, in order to allow for an direct comparison.
  - Requires the output of `HMI_Experiments/TestEvolutionMethods/Scaled/Analysis_EUV/analizeACWEscaledMagToEUV.py` for the specified image force and CR.
- `Analysis_EUV/analizeACWEscaledMagToEUV.py`: Compare the full-scale EUV-only segmentations of [CH-ACWE](https://github.com/DuckDuckPig/CH-ACWE) and full-scale EUV+HMI magnetogram segmentations developed here.
  - The variables in the `Key Variables` cell (`In[2]`) will need to be adjusted to point to the correct directories and chosen CR.
  - The user will also need to update the `magMethod` variable (also in the `Key Variables` cell) to select the image force employed for the HMI magnetogram data.
- `Analysis_EUV/visualization-BW.ipynb`: Provides box-and-whisker plots for the the output of `HMI_Experiments/TestEvolutionMethods/Scaled/Analysis_EUV/analizeACWEscaledMagToEUV.py`.

#### Analysis of EUV+HMI Magnetogram Segmentations as a Function of Spatial Resolution
- `Analysis_HMI/AreaChecks/Empty Segs.ipynb`: Report the names of one-eighth scale ("standard") segmentations that are empty. User will need to update the variables in the `Key Variables` cell (`In[2]`) to match both the correct directories and the image force employed for the HMI magnetogram data.
- `Analysis_HMI/AreaChecks/Large Area Look.ipynb`: Report and display one-eighth scale ("standard") cases that have a very large area identified as pertaining to a CH. User will need to update the variables in the `Key Variables` cell (`In[2]`) to match both the correct directories and the image force employed for the HMI magnetogram data.
- `Analysis_HMI/Magnetic Unipolarity Scaling Samples.ipynb`: Display examples of EUV+HMI magnetogram segmetnations at different spatial resolutions. User will need to update the variables in the `Key Variables` cell (`In[2]`) to match both the correct directories and the image force employed for the HMI magnetogram data.
- `Analysis_HMI/Single CR Lowest IOU Cases with Seed.ipynb`: Display cases where the full-scale EUV+HMI magnetogram segmentations and the one-eighth-scale EUV+HMI magnetogram segmentations differ the most.
  - The variables in the `Key Variables` cell (`In[2]`) will need to be adjusted to point to the correct directories and chosen CR.
  - The user will also need to update the `magMethod` variable (also in the `Key Variables` cell) to select the image force employed for the HMI magnetogram data.
  - Requires the output of `HMI_Experiments/TestEvolutionMethods/Scaled/Analysis_HMI/analizeACWEscaledMag.py` for the specified image force and CR.
- `Analysis_HMI/analizeACWEscaledMag.py`: Compare the full-scale EUV+HMI magnetogram segmentations with the one-eighth-scale EUV+HMI magnetogram segmentations.
  - The variables in the `Key Variables` cell (`In[2]`) will need to be adjusted to point to the correct directories and chosen CR.
  - The user will also need to update the `magMethod` variable (also in the `Key Variables` cell) to select the image force employed for the HMI magnetogram data.
- `Analysis_HMI/visulization_*-BW*.ipynb`: Provides box-and-whisker plots for the the output of `HMI_Experiments/TestEvolutionMethods/Scaled/Analysis_HMI/analizeACWEscaledMag.py`.
