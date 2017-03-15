clear all
close all
% this is an alternative version based on Fair et al 2007 Neuroimage dx.doi.org/10.1016/j.neuroimage.2006.11.051

addpath(genpath('bramila/'))



% load rois
load('rewardRois.mat')
R=length(rois);
ids=find(triu(ones(R),1));
% for each subject

subjbasepath='/scratch/eglerean/food/dataout/';
%subjbasepath='/m/nbe/scratch/braindata/eglerean/food/dataout/';
NS=34;
Nruns=2;
T=430;
TR=2;

% these were for block design
% HDL=5; % haemodynamic lag in seconds
% BLlen = 0; % number of seconds used from baseline, before trial starts
HRF=bramila_hrf(TR);


% extract roi time series
OVERWRITE=0;
allFD=zeros(T,Nruns,NS); % for storing framewise displacement
alllens=[];
for s = 1:NS
    for runid=1:Nruns
        % file to store the roi time series
        infile=[subjbasepath '/' num2str(s) '/' num2str(runid) '/roists.mat'];
        
        if(exist(infile)~=2 || OVERWRITE == 1) % if we are here, we need to do preprocessing
           	disp(['>>> Performing preprocessing for functional connectivity']);
			disp(['Creating file ' infile])
            tempinfile=[subjbasepath num2str(s) '/' num2str(runid) '/epi_preprocessed_FSLMNI.nii']
            cfg=[];
            cfg.StdTemplate='/home/glereane/code/bramila/external/MNI152_T1_2mm_brain_mask.nii';
            %cfg.StdTemplate='/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/external/MNI152_T1_2mm_brain_mask.nii';
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
			disp(['    Detrending...']);
            dtdata=bramila_detrend(cfg);
            cfg.vol=dtdata;
			disp(['    Regressing motion params'])
            cfg=bramila_motionregress(cfg);
            save_nii(make_nii(cfg.vol),[outpath num2str(runid) '.nii'])
            % filtering
			disp(['    Temporal filtering'])
            cfg= bramila_filter(cfg);
            
            cfgtemp.vol=cfg.vol;
            cfgtemp.rois=rois;
	    cfgtemp.usemean=1;
			disp(['    Extracting ROIs timeseries'])
            roits=bramila_roiextract(cfgtemp);
            save(infile,'roits');
        else
            disp(['File ' infile ' exists, loading it.'])
            load(infile); % variable roits
        end
        roits=zscore(roits);
        
        % now load regressor time series
		
		disp(['>>> Loading regressors'])
        for c=1:4 % four conditions
            
            temp=load([subjbasepath '/' num2str(s) '/' num2str(runid) '/' num2str(runid) '-' num2str(c) '.txt']);
            regs{c}=temp;
		end
		cfg=[];
		cfg.y=roits;
		cfg.events=regs;
		cfg.TR=TR;
		cfg.hrf=HRF;

        betas=bramila_betaseries(cfg);
        
        disp(['>>> Calculating Framewise Displacement'])
        motionfile=[subjbasepath num2str(s) '/' num2str(runid) '/rp.txt'];        
        cfgtemp=[];
        cfgtemp.motionparam=motionfile;
        cfgtemp.prepro_suite='spm';
        FD=bramila_framewiseDisplacement(cfgtemp);
        allFD(:,s,runid)=FD;

		disp(['>>> Calculating connectivity matrices'])
        outlabels={'appetising'
            'bland'
            'cars'
            'fixations' };
        for cc=1:2 % only conditions 1 and 2
          	% regress the other conditions and use residual to build network
			betats=betas{cc};
			% for this subject and this run there is now a matrix with beta time series as events and columns = num of rois
			
            avgnet_M=corr(betats);
            avgnet_par_M=partialcorr(betats);
            
            % compute null distr for avg and max
			disp(['>>> Compute individual null distribution using surrogate roi time series and the max statistics'])
            subjperms=zeros(10000,1);
            X=fft(zscore(betats));
            XA=abs(X);
            XP=angle(X);
            for i=1:10000
				if(mod(i-1,1000)==0) fprintf(['..' num2str(i)]); end
                % resynth all rois
                xperm=zeros(size(betats));
                for r=1:R
                    XPtemp=XP(randperm(length(XP)),r);
                    xperm(:,r)=zscore(real(ifft(XA(:,r).*XPtemp)));
                end
                surrnet=corr(xperm);
                temp=surrnet(ids);
                subjperms(i)=max(temp);
            end
            
            disp(['Saving ' subjbasepath num2str(s) '/run' num2str(runid) '_cond' num2str(cc) '_nets_betaseries.mat']);
            save([subjbasepath num2str(s) '/run' num2str(runid) '_cond'  num2str(cc) '_nets_betaseries.mat'],'avgnet_M','avgnet_par_M','subjperms')
            %if(any(mfd>0.5/len)) disp(num2str(mfd)); end
            clear avgnet_M
            clear avgnet_par_M
			clear subjperms
        end % end of cycle that stores the networks
        
        
        
    end % end of cycle for the runs
end % end of cycle for the subjects
save([subjbasepath '/FCsession_betaseries.mat'], 'allFD',  'rois', 'R', 'ids', 'NS', 'Nruns');
