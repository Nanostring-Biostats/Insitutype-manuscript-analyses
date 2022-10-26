This repo holds analyses used in the manuscript, "Insitutype: likelihood-based cell typing for single cell spatial transcriptomics".

Note: the datasets used here are too big to host on github. They can be downloaded from: https://nanostring.box.com/s/cx5j3dlrmext3guyy1bb2r4rm6p8d364.

How to use:
- First, download the data from the above url and place in the "data" folder. 
- Each folder holds a distinct analysis (a figure or sub-figure)
- All R scripts execute in the directory in which they lie
- These analyses use a development version of Insitutype, "MLEcell". This package can be installed using the "environment_setup" script.
- The benchmarking analysis (Fig 4) uses a more final version of the package, "InSituType". This is also installed in the "environment_setup" script.

