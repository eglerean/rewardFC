clear all
close all

addpath('/m/nbe/scratch/braindata/shared/toolboxes/spm12/')

files=dir('*.mat');

%% add here code to go through each roi, conver to nii, split left right, save in final nii file with each integer for each cluster, and then call bramila_makeRoiStruct

rois=zeros(91,109,91,length(files)*2);

counter=1;
maskL=zeros(81,109,91);
maskL(1:45,:,:)=1;
maskR=ones(91,109,91)-maskL;

for f=1:length(files);
	label=strrep(files(f).name,'_roi.mat','');
    	mars_rois2img(files(f).name,'temp_roi.nii')
	nii=load_nii('temp_roi.nii');
	% split left right
	imgL=nii.img.*maskL;
	imgR=nii.img.*maskR;

	rois(:,:,:,f*2-1)=imgL;
	rois(:,:,:,f*2)=imgR;
end
