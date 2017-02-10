clear all
close all
addpath(genpath('bramila/'))

if(0)
%% testing regressors
TR=2;
T=430;
HRF=bramila_hrf(TR);
regts=zeros(T,4);
stickregts=zeros(T,4);
for c=1:4 % four conditions
    subj=1;
    runid=1;
    temp=load(['/m/nbe/scratch/braindata/eglerean/food/data/1/regressors/' num2str(subj) '-' num2str(runid) '-' num2str(c) '.txt']);
    % from temp to time series
    tempregts=zeros(T,1);
    for t=1:length(temp)
       tempregts(round(temp(t)/TR))=1; % stick function for the event
    end
    stickregts(:,c)=tempregts;
    
    tempregts=conv(tempregts,HRF);
    regts(:,c)=tempregts(1:T);
end
regtsID=zeros(T,1);
for t=1:T
    vec=regts(t,:);
    temp=find(max(vec)==vec);
    win=temp(end);
    regtsID(t)=win; % this will be used to finf Times Of Interest (TOI)
end

    error('stop')
end

%%


% load rois
load('rewardRois.mat')
R=length(rois);
ids=find(triu(ones(R),1)); 
% for each subject

subjbasepath='/home/scratch/eglerean/food/dataout/';
NS=5;
Nruns=2;
T=430;
TR=2;

% these were for block design
% HDL=5; % haemodynamic lag in seconds
% BLlen = 0; % number of seconds used from baseline, before trial starts


% extract roi time series  
OVERWRITE=1;
allFD=zeros(T,Nruns,NS); % for storing framewise displacement
alllens=[];
for s = 1:NS
	for runid=1:Nruns
        % file to store the roi time series
		infile=[subjbasepath num2str(s) '/' num2str(runid) '_rois.mat'];
        
		if(exist(infile)~=2 || OVERWRITE == 1) % if we are here, we need to do preprocessing
			disp(['Creating file ' infile])
			tempinfile=[subjbasepath num2str(s) '/' num2str(runid) '/epi_preprocessed_FSLMNI.nii']
			cfg=[];
			cfg.StdTemplate='/home/VSSHP/glereane/code/bramila/external/MNI152_T1_2mm_brain_mask.nii';
			cfg.fileID=[num2str(s) '-' num2str(runid)];
			cfg.infile=tempinfile;
			cfg.TR=TR;
			cfg.mask=ones([91 109 91]);
			% add motion regression
			motionfile=[subjbasepath num2str(s) '/' num2str(runid) '/rp.txt'];
			disp(motionfile)
			cfg.motionparam=motionfile;
			cfg.write=0;
			cfg.motion_reg_type='volterra';
			cfg.detrend_type='Savitzky-Golay';
			outpath=[subjbasepath num2str(s) '/' num2str(runid) '/mot/'];
			mkdir(outpath);
			cfg.outpath=outpath;
			cfg.inpath=outpath;
			dtdata=bramila_detrend(cfg);
			cfg.vol=dtdata;
			cfg=bramila_motionregress(cfg);
			save_nii(make_nii(cfg.vol),[outpath num2str(runid) '.nii'])
			% filtering
			cfg= bramila_filter(cfg);
			
			cfgtemp.vol=cfg.vol;
			cfgtemp.rois=rois;
			roits=bramila_roiextract(cfgtemp);
			save(infile,'roits');	
			error('stop here please');
		else
			disp(['File ' infile ' exists.'])
			load(infile); % variable roits
		end
		roits=zscore(roits);
        
        % now load regressor time series
        
		reg=load(['/triton/becs/scratch/braindata/eglerean/motfmri/motfmri/regressor_files_mitmot/' num2str(s) '-' num2str(runid) '.txt']);
		[aa bb]=corr(reg(:,1),[1:size(reg,1)]')
        
		len=ceil((min(reg(find(reg(:,4)>0),3))+BLlen)/TR); % for all the non zero trials, find the shortest interval in TRs 
														% Note the +BLlen for the BLlen seconds of baseline before trial starts
		ends=ceil((HDL+reg(:,2)+reg(:,3))/1.7); % for each trial compute the ending time plus haemodyn delay in TRs
		windows=[ends-len ends]; % windows in TRs
		
		motionfile=[subjbasepath num2str(s) '/' num2str(runid) '.txt'];
        
        
		cfgtemp=[];
		cfgtemp.motionparam=motionfile;
		cfgtemp.prepro_suite='spm';
		FD=bramila_framewiseDisplacement(cfgtemp);
		allFD(:,s-3)=FD;
		alllens(end+1)=len;
		outlabels={'mit','mot'};
		for mm=1:2 % 1 == mit, 2 == mot
			for cc=[0 2 4 6]
				if(cc<=4)
					rows=intersect(find(reg(:,1)==mm),find(reg(:,4)==cc));
				else
					rows=intersect(find(reg(:,1)==mm),find(reg(:,4)>0));
				end
				nets=zeros(length(ids),length(rows));
				mfd=zeros(length(rows),1);
				me=zeros(length(rows),R);
				toi=[];
				for r=1:length(rows)
					net=corr(roits([windows(rows(r),1):windows(rows(r),2)],:));

					me(r,:)=mean(roits([windows(rows(r),1):windows(rows(r),2)],:));
					toi=[toi windows(rows(r),1):windows(rows(r),2)];
					nets(:,r)=net(ids);
					mfd(r)=mean(FD(windows(rows(r),1):windows(rows(r),2)));
				end
				disp(num2str(length(toi)))
				avgnet=corr(roits(toi,:));
				avgnet_par=partialcorr(roits(toi,:));
				timemask=zeros(size(roits,1),1);
				timemask(toi)=1;
				timemask=repmat(timemask,1,R);
				mfdM=mean(FD(toi));
				if(1)
					% don't correct the average
					avgnet_M=corr(timemask.*roits);
					avgnet_par_M=partialcorr(timemask.*roits);
				else
					% correct the average
					temproits=roits;
					for r=1:R					
						temproits(find(timemask(:,1)==0),r)=mean(temproits(find(timemask(:,1)==1),r));
						avgnet_M=corr(temproits);
						avgnet_par_M=partialcorr(temproits);
					end
				end
				
				% compute null distr vor avg and max
				subjperms=zeros(50000,1);
				X=fft(zscore(roits));
				XA=abs(X);
				XP=angle(X);
				for i=1:50000
					% resynth all rois
					xperm=zeros(size(roits));
					for r=1:R
						XPtemp=XP(randperm(length(XP)),r);
						xperm(:,r)=zscore(real(XA(:,r).*XPtemp));	
					end
					surrnet=corr(timemask.*xperm);
					temp=surrnet(ids);
					subjperms(i)=max(temp);
				end

				disp(['Saving ' subjbasepath num2str(s) '/' num2str(runid) '_' outlabels{mm} '_' num2str(cc) '_nets.mat']);
				save([subjbasepath num2str(s) '/' num2str(runid) '_' outlabels{mm} '_' num2str(cc) '_nets.mat'],'me','nets','ids','mfd','avgnet','avgnet_par','toi','avgnet_M','avgnet_par_M','mfdM','subjperms')
				%if(any(mfd>0.5/len)) disp(num2str(mfd)); end
				clear nets
				clear mfd
				clear toi
				clear avgnet
				clear avgnet_par
				clear avgnet_M
				clear avgnet_par_M
				clear me
			end

		end 


	end
end	
save FCsession allFD alllens rois R ids NS Nruns
