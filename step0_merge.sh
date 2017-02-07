#!/bin/bash
basepath='/home/scratch/lauri/RewardPET/';
outpath='/home/scratch/eglerean/dataout/';

# set FSL output to NII without gz
export FSLOUTPUTTYPE=NIFTI

# First we find which folders contain data. There are two runs per subject. Each run has length 430
find $basepath -type f -name "sw*_430.nii"|sort -n|cut -d\/ -f1-8 > $outpath/folderlist.txt

subNum=1;
for f in $(cat $outpath/folderlist.txt); do
	# get subject ID from the study, needed for the regressors
	subjID=$(echo $f|cut -d\/ -f7);
	# get the run number for the subject
	runID=$(echo $f|cut -d\/ -f8|sed 's/[A-z]*//g')
	echo mkdir -p $outpath/$subNum/$runID/
	# identify string for subject data
	id=$(find $f|grep sw.*0001.nii);
	ide=$(echo $id|sed 's/0001.nii//g');
	# identify motion parameters
	mot=$(find $f|grep rp.*txt);
	echo cp $mot $outpath/$subNum/$runID/
	echo /scratch/braindata/eglerean/motfmri/git/merge_it.sh $ide 1 430 $outpath/$subNum/$runID/epi_prerocessed.nii  			
	if [ $? -ne 0 ]; then
		echo "error. stopping"
		exit 1	
	fi  
	let subNum=subNum+1;
done
