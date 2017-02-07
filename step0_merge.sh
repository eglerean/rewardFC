#!/bin/bash
basepath='/home/scratch/lauri/RewardPET/';
outpath='/home/scratch/eglerean/food/dataout/';

# set FSL output to NII without gz
export FSLOUTPUTTYPE=NIFTI

# First we find which folders contain data. There are two runs per subject. Each run has length 430
find $basepath -type f -name "sw*_430.nii"|sort -n|cut -d\/ -f1-8 > $outpath/folderlist.txt

subNum=0;
prevsubjID=-1;
for f in $(cat $outpath/folderlist.txt); do
	# get subject ID from the study, needed for the regressors
	subjID=$(echo $f|cut -d\/ -f7);
	if [ $prevsubjID -ne $subjID ]; then
		let subNum=subNum+1;
		prevsubjID=$subjID;
	fi
	# get the run number for the subject
	runID=$(echo $f|cut -d\/ -f8|sed 's/[A-z]*//g')
	echo mkdir -p $outpath/$subNum/$runID/
	#mkdir -p $outpath/$subNum/$runID/
	
	# identify string for subject data
	id=$(find $f|grep sw.*001.nii);
	ide=$(echo $id|sed 's/001.nii//g');
	
	# identify motion parameters
	mot=$(find $f|grep rp.*txt);
	echo cp $mot $outpath/$subNum/$runID/rp.txt

	# copy the regressor file
	regloc=$(echo $f|cut -d\/ -f1-6);
	regloc=$regloc"/RegressorFil*/";
	echo cp $regloc/$subjID-$runID-* $outpath/$subNum/$runID/	

	# do the actual merge
	echo merge_it.sh $ide 1 430 $outpath/$subNum/$runID/epi_prerocessed.nii 2
	if [ $? -ne 0 ]; then
		echo "error. stopping"
		exit 1	
	fi 
done
