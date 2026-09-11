# Software_PolMFM
MATLAB software for processing and analysing PolMFM single-molecule orientation and localization microscopy data. 

## Requirements
### OS supported
Windows 10 and 11
### System tested
* Windows 10 with processor Intel(R) Xeon(R) CPU E5-1650 and Graphics Card	NVIDIA GeForce GTX 1630 (4 GB)
* Windows 11 with processeur 11th Gen Intel(R) Core(TM) i7-1165G7 (2.80 Ghz) and Graphics Card Intel(R) Iris(R) Xe Graphics (128 Mo)
* Windows 11 with processor 13th Gen Intel(R) Core(TM) i7-13700 (2.10 GHz) and Graphics Cards NVIDIA GeForce RTX 3050 OEM (8 GB) and Intel(R) UHD Graphics 770 (128 MB)

### Version
The code was tested in matlab R2021a (9.10), matlab R2022a (9.12) and matlab R2025b (25.2). 

### Required tollkits:

'Signal Processing Toolbox'	'8.6'

'Image Processing Toolbox'	'11.3'

'Statistics and Machine Learning Toolbox'	'12.1'

'Curve Fitting Toolbox'	'3.5.13'

'Parallel Computing Toolbox'	'7.4' (optional for parallel computation)

'MATLAB Parallel Server'	'7.4'(optional for parallel computation)

'Polyspace Bug Finder'	'3.4'

## Tutorial
A step-by-step tutorial describing the analysis workflow is provided in this repository: **Pipeline_PolMFM-tutorial.pdf**.

## Data
Image data should be in Tiff format.

To test the software, you must first download a testing dataset that can be found in the following Zenodo repository: _https://zenodo.org/records/22043565_. Download and unzip  the folder **0 - DATA_Tutorial** and copy-paste it into the main GitHub folder. The testing dataset includes:
- Raw PolMFM data of Nile Red-labelled Large Unilamellar Vesicles (LUVs): **DATA_SilicaBeads_SLB_NR_2024-01-05_002_substack**
- Corresponding calibration data: **Calibrations**

## Dependancies

This software includes parts of FISH-quant software, used for prelocalization and 3D Gaussian fitting(Mueller, F., Senecal, A., Tantale, K. et al. FISH-quant: automatic counting of transcripts in 3D FISH images. Nat Methods 10, 277–278 (2013). https://doi.org/10.1038/nmeth.2406). 

## License and Citation

