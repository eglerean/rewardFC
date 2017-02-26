# Functional connectivity for "Reward PET" study

Analysis is divided into steps. Order below.

## Step 0
### from SPM preprocessed data to NIFTI 4D files

Run BASH script step0_merge.sh. This script also takes care of copying motion parameters that will be used at step 1.
Furthermore, create the ROIs using the matlab script inside the ROIs subfolder. It requires marsbars, installed with
```git clone https://github.com/matthew-brett/marsbar```

## Step 1 - further preprocessing for connectivity




