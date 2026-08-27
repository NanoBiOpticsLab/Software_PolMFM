# Software_PolMFM
MATLAB software for processing and analysing PolMFM single-molecule orientation and localization microscopy data. 

A step-by-step tutorial describing the analysis workflow is provided in this repository: Pipeline_PolMFM-tutorial.pdf.

To test the software, you must first download a testing dataset that can be found in the following Zenodo repository: _https://zenodo.org/records/22043565_. Download and unzip  the folder **0 - DATA_Tutorial** and copy-paste it into the main GitHub folder. The testing dataset includes:
- Raw PolMFM data of Nile Red-labelled Large Unilamellar Vesicles (LUVs): **DATA_SilicaBeads_SLB_NR_2024-01-05_002_substack**
- Corresponding calibration data: **Calibrations**

This software is based on FISH-quant, a previously published software in Matlab for quantifying dense mRNA aggregates at transcription sites in three dimensions (Mueller, F., Senecal, A., Tantale, K. et al. FISH-quant: automatic counting of transcripts in 3D FISH images. Nat Methods 10, 277–278 (2013). https://doi.org/10.1038/nmeth.2406). 
