basepath='/scratch/braindata/eglerean/motfmri/niidata/';
for n in $(seq 4 26); do
	if [ -d aa$n ]; then
		echo "Skipping $n"
	else
		#extract in temp
		echo tar xvf $basepath$n.tar -C $basepath/temp/
		tar xvf $basepath$n.tar -C $basepath/temp/
		# process mot1 to 3
		for m in $(seq 1 3); do
			echo mkdir -p $basepath/$n/
			mkdir -p $basepath/$n/
			id=$(find $basepath/temp/$n/mitmot$m|grep sw.*0001.nii)
			ide=$(echo $id|sed 's/0001.nii//g');
			mot=$(find $basepath/temp/$n/mitmot$m|grep rp.*txt);
			echo cp $mot $basepath/$n/$m.txt
			cp $mot $basepath/$n/$m.txt
			echo /scratch/braindata/eglerean/motfmri/git/merge_it.sh $ide 1 553 $basepath/$n/$m.nii  			
			/scratch/braindata/eglerean/motfmri/git/merge_it.sh $ide 1 553 $basepath/$n/$m.nii
			if [ $? -ne 0 ]; then
				echo "error. stopping"
				exit 1	
			fi  
			echo gunzip $basepath/$n/$m.nii.gz		
			gunzip $basepath/$n/$m.nii.gz
		done
		echo rm -r temp/*
		rm -r temp/*
	fi
done
