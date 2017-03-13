%% 	Runs analysis using windowed correlations for each stimulus block
%	First level analysis is done for each run for each subject. Three runs are combined as averages
%	Second level analysis also as group average with permutations
%	Option LINKS = 1: performs contrasts on the edges of the network
%	Option LINKS = 0: performs contrasts on the mean amplitude of the block (i.e. GLM)


clear all
close all
addpath(genpath('bramila/'))


subjbasepath='/m/nbe/scratch/braindata/eglerean/food/dataout/';
load([ subjbasepath '/FCsession.mat']) % variables
NS=34;
Nruns=2;
Ncond=2; % number of conditions
NITER=50000;
LINKS=1;
if(LINKS == 1)
    NC= length(ids); % number of comparisons
else
    NC = length(rois); % not sure if this is implemented
end
grouptvals=zeros(NC,1);

outlabels={'appetising'
    'bland'
    'cars'
    'fixations' };

%% 1st level analysis, loading data
all_data=zeros(Ncond,Nruns,NS,length(ids)); % all network data
all_subjperms=zeros(Ncond,Nruns,NS,NITER); % all network permutations
all_mFD=zeros(Ncond,Nruns,NS);
for s=1:NS
    for r=1:Nruns
        disp(['Subject ' num2str(s) ' Run ' num2str(r)])
        for cc=1:Ncond % 1 == mit, 2 == mot
            
            temp=load([subjbasepath '/'  num2str(s) '/run' num2str(r) '_cond' num2str(cc) '_nets.mat']);
            net=atanh(temp.avgnet_M);
            all_subjperms(cc,r,s,:)=atanh(temp.subjperms);
            all_data(cc,r,s,:)=net(ids);
            all_mFD(cc,r,s)=mean(temp.mfdM); % it's already computing the mean
        end
    end
end


% one sample main effects (average connecitivity)
allMaxTH=zeros(Ncond,Nruns);
for cc=1:Ncond
    for r=1:Nruns
        tempdata=squeeze(all_data(cc,r,:,:)); % subj x links
        tempdata=tempdata'; % 1 subj per column, they are already atanh
        
        tempperms=squeeze(all_subjperms(cc,r,:,:)); % get all fake r values
        group_perms=mean(tempperms,1);
        stats.mean=mean(tempdata,2);
        stats.maxTH=prctile(group_perms,95);
        allMaxTH(cc,r)=stats.maxTH;
        tempmask=double(stats.mean>stats.maxTH(1));
        net=zeros(R);
        net(ids)=tempmask.*stats.mean;
        net=net+net';
        net=tanh(net);
        templabel=['Avgerage ' outlabels{cc} ' ' num2str(r) '-net'];
        disp([templabel ' ' num2str(sum(tempmask))]);
        figure
        bramila_plotConn(net,rois,templabel,[-1 1],1);
        templabel=['Avg_' outlabels{cc} '_' num2str(r) '_net'];
        saveas(gcf,['pngs/' templabel '.png' ])
        save(['mats/' templabel '.mat'],'net','rois')
    end
end
save allMaxTH allMaxTH


% paired contrasts
for r=1:Nruns
   
    dataA=squeeze(all_data(1,r,:,:))';
    dataB=squeeze(all_data(2,r,:,:))';
    
    % regress framewise displ
    tempdata=[dataA dataB];
    tempdata=tempdata';
    temp_mFD=[squeeze(all_mFD(1,r,:));squeeze(all_mFD(2,r,:))];
    disp('regressing FD')
    %for i=1:NC
    %    [aa bb res]=regress(tempdata(:,i),[temp_mFD  ones(length(temp_mFD),1)]);
    %    tempdata(:,i)=res;
    %end
    
    tempdata=tempdata';
    design=[ones(1,size(dataA,2)) 2*ones(1,size(dataB,2))];
    disp('computing permutations')
    [runstats]=bramila_ttest2_np([tempdata],design,5000);
    grouptvals(:,r)=runstats.tvals;
    grouppvals(:,r)=min(runstats.pvals,[],2);
end


if(LINKS == 1)
    % plot connectivity matrices
    mask=double(grouppvals<0.05) + double(grouppvals>0.95);
    tempgrouptvals=mask.*grouptvals;
    
    for r=1:Nruns
        tempP=grouppvals(:,r);
        tempFDR=mafdr(tempP,'BHFDR','true');
        tempmask=double(tempFDR<=0.05);
        net=zeros(R);
        temptvals=tempgrouptvals(:,r).*tempmask;
        %net(ids)=tempgrouptvals(:,co);
        net(ids)=temptvals;
        uu=unique(abs(temptvals));
        if(length(uu)>1)
            uu(2)
        end
        net=net+net';
        bramila_plotConn(net,rois,'Appetising minus bland',[-4 4],6);

        saveas(gcf,['pngs/avg_M_AppVsBla.png' ])
        save(['mats/avg_M_AppVsBla.mat'],'net','rois')
        
    end
    
else
    % display summary table
    ppp=grouppvals<0.05/12;
    for ll=1:12;disp([sprintf('%-*s', 5, rois(ll).label) ' ' num2str(ppp(ll,:))]);end
    
end

