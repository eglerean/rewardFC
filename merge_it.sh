## bash function that merges data preprocessed with SPM into 4D nifti
## Uses fslmerge
#!/bin/bash

# $1 = swMitmot_Ami_01_20140320_001_002_EPI_64_mit_mot_544_
# $2 = 1
# $3 = 544
# $4 = mitmot1.nii
# $5 = TR 1.7

# comments: note that there is basically no input validation

if [ $# -lt 4 ]; then 
	echo "missing argument"
	exit 1
fi

args='';
for i in $(seq -f '%04g' $2 $3); do
	args="$args $1$i".nii;
done
echo fslmerge -t $4 $args
fslmerge -t $4 $args
