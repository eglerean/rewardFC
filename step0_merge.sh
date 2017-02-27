#!/bin/bash
basepath='/scratch/lauri/RewardPET/';
outpath='/scratch/eglerean/food/dataout/';

# set FSL output to NII without gz
export FSLOUTPUTTYPE=NIFTI

# First we find which folders contain data. There are two runs per subject. Each run has length 430
find $basepath -type f -name "sw*_430.nii"|sort -n|cut -d\/ -f1-7 > $outpath/folderlist.txt

subNum=0;
prevsubjID=-1;
for f in $(cat $outpath/folderlist.txt); do
	 
	# get subject ID from the study, needed for the regressors
	subjID=$(echo $f|cut -d\/ -f6);
	if [ $prevsubjID -ne $subjID ]; then
		let subNum=subNum+1;
		prevsubjID=$subjID;
	fi
	

	#limit analysis to few subj for now
	if [ $subNum -gt 5 ]; then
		break
	fi
	# get the run number for the subject
	runID=$(echo $f|cut -d\/ -f7|sed 's/[A-z]*//g')
	echo mkdir -p $outpath/$subNum/$runID/
	mkdir -p $outpath/$subNum/$runID/
	
	# identify string for subject data
	id=$(find $f|grep sw.*001.nii);
	ide=$(echo $id|sed 's/001.nii//g');
	
	# identify motion parameters
	mot=$(find $f|grep rp.*txt);
	echo cp $mot $outpath/$subNum/$runID/rp.txt
	cp $mot $outpath/$subNum/$runID/rp.txt



	# copy the regressor file
	regloc=$(echo $f|cut -d\/ -f1-5);
	regloc=$regloc"/RegressorFil*/";
	for file in $(ls $regloc/$subjID-$runID-*.txt); do
		filename=$(echo $file|cut -d\/ -f7);
		echo $filename;
		newfilename=$(echo $filename|cut -d- -f2-);
		cat $file|tr '' '\n' >  $outpath/$subNum/$runID/$newfilename
		#cp $regloc/$subjID-$runID-* $outpath/$subNum/$runID/	
	done

	# convert regressor to unix format
	ls $outpath/$subNum/$runID/


	# do the actual merge
	echo merge_it.sh $ide 1 430 $outpath/$subNum/$runID/epi_preprocessed.nii 2
	./merge_it.sh $ide 1 430 $outpath/$subNum/$runID/epi_preprocessed.nii 2
	if [ $? -ne 0 ]; then
		echo "error. stopping"
		#exit 1	
	fi

	#convert the output to MNI FSL boundary box
	mni2mm=$FSLDIR"/data/standard/MNI152_T1_2mm_brain.nii.gz"
	ref=$mni2mm;

	echo flirt -in $outpath/$subNum/$runID/epi_preprocessed.nii -ref $ref -out $outpath/$subNum/$runID/epi_preprocessed"_FSLMNI.nii" -init spm2fsl2mm.mat -applyxfm
	flirt -in $outpath/$subNum/$runID/epi_preprocessed.nii -ref $ref -out $outpath/$subNum/$runID/epi_preprocessed"_FSLMNI.nii" -init spm2fsl2mm.mat -applyxfm
	rm -f $outpath/$subNum/$runID/epi_preprocessed.nii
done
