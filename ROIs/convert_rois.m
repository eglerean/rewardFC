clear all
close all

addpath('/m/nbe/scratch/braindata/shared/toolboxes/spm12/')
addpath('./marsbar/marsbar')
addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/'))
files=dir('*.mat');

%% add here code to go through each roi, conver to nii, split left right, save in final nii file with each integer for each cluster, and then call bramila_makeRoiStruct

allrois=zeros(91,109,91,length(files)*2);

counter=1;
maskL=zeros(91,109,91);
maskL(1:45,:,:)=1;
maskR=ones(91,109,91)-maskL;

for f=1:length(files);
	label=strrep(files(f).name,'_roi.mat','');
    mars_rois2img(files(f).name,'trash/temp_roi.nii')
	nii=load_nii('trash/temp_roi.nii');
	% split left right
	imgL=nii.img.*maskL;
	imgR=nii.img.*maskR;

	allrois(:,:,:,f*2-1)=imgL;
	allrois(:,:,:,f*2)=imgR;
    alllabels{f*2-1,1}=[label '_L'];
    alllabels{f*2,1}=[label '_R'];
end

% get rid of voxels that have more than one roi

temp=sum(allrois,4);

mask=ones(size(temp));
mask(find(temp>1))=0;

% combine all rois as indexed single volume
roivol=zeros(size(temp));
for c=1:size(allrois,4);
    roivol=roivol+c*mask.*allrois(:,:,:,c);
end
save_nii(make_nii(roivol,[2 2 2]),'roimask.nii');

cfg=[];
cfg.roimask='roimask.nii';
cfg.labels=alllabels;
cfg.imgsize=[91 109 91];
rois = bramila_makeRoiStruct(cfg);
save ../rewardRois rois