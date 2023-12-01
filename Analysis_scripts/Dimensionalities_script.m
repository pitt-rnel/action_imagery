%% Dimensionality (Figure 2f,g and Figure 5f,g)
[FSCREE,QAT,QIT,S_allaction,S_allimagery,DS,QT,S_BOTH] = deal(cell(2,1)); 
dimen = @(x) find(cumsum(x)> 90,1,'first'); 
tdim_thresh = 1; % variance threshold for dimensionality estimation
nboot = 1000; 
rng(1); 

[Rdyn, AAv,SAv,IIv,SIv,IAv,AIv,...
    AAv_f,SAv_f,IIv_f,SIv_f,IAv_f,AIv_f,...
    TVARS,FVARS,TAv,TIv,TS,AIDZ] = deal(cell(2,1)); 
for z = 1:2
    LiHis = {1:size(ZsplitAction{z}{1},1)/2, (size(ZsplitAction{z}{1},1)/2+1):size(ZsplitAction{z}{1},1)}; 
    if z == 2; LiHis = cellfun(@(x) x(20:end),LiHis,'uni',0); end
    lh = [LiHis{1},LiHis{2}]'; 
    
    fAA_sess = cell2mat(reshape(cellfun(@(x) (x(LiHis{2},:)-x(LiHis{1},:)),ZsplitAction{z}(:,1),'uni',0),1,1,[])); 
    fII_sess = cell2mat(reshape(cellfun(@(x) (x(LiHis{2},:)-x(LiHis{1},:)),ZsplitImagery{z}(:,2),'uni',0),1,1,[]));  
    fSA_sess = cell2mat(reshape(cellfun(@(x) (x(LiHis{2},:)-x(LiHis{1},:)),ZsplitShared{z}(:,1),'uni',0),1,1,[]));  
    fSI_sess = cell2mat(reshape(cellfun(@(x) (x(LiHis{2},:)-x(LiHis{1},:)),ZsplitShared{z}(:,2),'uni',0),1,1,[]));  
    fAI_sess = cell2mat(reshape(cellfun(@(x) (x(LiHis{2},:)-x(LiHis{1},:)),ZsplitAction{z}(:,2),'uni',0),1,1,[])); 
    fIA_sess = cell2mat(reshape(cellfun(@(x) (x(LiHis{2},:)-x(LiHis{1},:)),ZsplitImagery{z}(:,1),'uni',0),1,1,[])); 
    
    tAA_sess = cell2mat(reshape(ZsplitAction{z}(:,1),1,1,[])); tAA_sess = tAA_sess(lh,:,:); 
    tII_sess = cell2mat(reshape(ZsplitImagery{z}(:,2),1,1,[])); tII_sess = tII_sess(lh,:,:); 
    tSA_sess = cell2mat(reshape(ZsplitShared{z}(:,1),1,1,[])); tSA_sess = tSA_sess(lh,:,:); 
    tSI_sess = cell2mat(reshape(ZsplitShared{z}(:,2),1,1,[])); tSI_sess = tSI_sess(lh,:,:); 
    tAI_sess = cell2mat(reshape(ZsplitAction{z}(:,2),1,1,[])); tAI_sess = tAI_sess(lh,:,:); 
    tIA_sess = cell2mat(reshape(ZsplitImagery{z}(:,1),1,1,[])); tIA_sess = tIA_sess(lh,:,:); 
    tT_sess = cat(1,cat(2,tAA_sess,tSA_sess,tIA_sess),cat(2,tAI_sess,tSI_sess,tII_sess)); 
    
    fACT_sess = cat(2,fAA_sess,fSA_sess,fIA_sess); 
    fIM_sess = cat(2,fAI_sess,fSI_sess,fII_sess); 
    tACT_sess = cat(2,cell2mat(reshape(cellfun(@(x) x(lh,:),ZsplitAction{z}(:,1),'uni',0),1,1,[])),...
                      cell2mat(reshape(cellfun(@(x) x(lh,:),ZsplitShared{z}(:,1),'uni',0),1,1,[])),...
                      cell2mat(reshape(cellfun(@(x) x(lh,:),ZsplitImagery{z}(:,1),'uni',0),1,1,[]))); 
    tIM_sess = cat(2,cell2mat(reshape(cellfun(@(x) x(lh,:),ZsplitAction{z}(:,2),'uni',0),1,1,[])),...
                     cell2mat(reshape(cellfun(@(x) x(lh,:),ZsplitShared{z}(:,2),'uni',0),1,1,[])),...
                     cell2mat(reshape(cellfun(@(x) x(lh,:),ZsplitImagery{z}(:,2),'uni',0),1,1,[]))); 
    
    [AID] = deal(cell(length(ActionImageryData.(sprintf('P%d',z+1))),1)); 
    pcn = size(Subspaces{z}{1},1); 
    for i = 1:length(AID)
        for j = 1:2
            AID{i}.Action{1,j} = cellfun(@(x) x(1:pcn,LiHis{1}),ActionImageryData.(sprintf('P%d',z+1)){i}.Action.PC_combined{j},'uni',0); 
            AID{i}.Imagery{1,j} = cellfun(@(x) x(1:pcn,LiHis{1}),ActionImageryData.(sprintf('P%d',z+1)){i}.Imagery.PC_combined{j},'uni',0); 
        end
    end
    AIDZ{z} = AID;     
    
    ActIm = cellfun(@(x) [x.Action,x.Imagery],AID,'uni',0); ActIm = vertcat(ActIm{:}); 
    bA = cell(size(ActIm)); 
    for i = 1:numel(bA)
        [~,bA{i}] = bootstrp(nboot,[],ActIm{i}); 
    end
    nT = size(AID{1}.Action{1}{1},2); Aind = 1:(nT*2); Iind = (nT*2+1):4*nT; 
    
    [AAv{z},SAv{z},IIv{z},SIv{z},IAv{z},AIv{z},...
        AAv_f{z},SAv_f{z},IIv_f{z},SIv_f{z},IAv_f{z},AIv_f{z},...
        TAv{z},TIv{z}] = deal(NaN(nboot,pcn)); 
    [TVARS{z},FVARS{z}] = deal(NaN(nboot,2)); 
    for i = 1:nboot
        clc; fprintf('%d/%d (%d/%d)\n',z,2,i,nboot); 
        Act_b = cell(size(ActIm)); 
        for j = 1:numel(ActIm)
            Act_b{j} = ActIm{j}(bA{j}(:,i)); 
        end
        Act_b_means = cellfun(@(x) nanmean(cell2mat(reshape(x,1,1,[])),3),Act_b,'uni',0); 
        Act_b_means = cellfun(@(x) x - nanmean(x(:,BIND{z}),2),Act_b_means,'uni',0);
        Act_b_cats = cell(size(Act_b_means,1),1); 
        for j = 1:size(Act_b_means,1)
            Act_b_cats{j} = cell2mat(Act_b_means(j,:))'; 
        end
%         [~,TS{z}{i},~,mu] = generalProcrustes(Act_b_cats,1:10,0);
        [~,~,~,mu] = generalProcrustes(Act_b_cats,1:10,0);
        muAI = {mu(Aind,:);mu(Iind,:)}; 
        muAIf = cellfun(@(x) x(((end/2)+1):end,:)-x(1:end/2,:),muAI,'uni',0); 
        
        Qboot = DySO(muAI,99,'do_plot',false); 
        %%% TOTAL
        [VSC1,~,va] = varimax_sparsity(muAI{1}*Qboot.unique.C1,[1:2,(1:2)+nT-1]); 
        [VSC2,~,vi] = varimax_sparsity(muAI{2}*Qboot.unique.C2,[1:2,(1:2)+nT-1]); 
        VSS1 = varimax_sparsity(muAI{1}*Qboot.shared,[1:2,(1:2)+nT-1]); 
        VSS2 = varimax_sparsity(muAI{2}*Qboot.shared,[1:2,(1:2)+nT-1]); 
        VST1 = varimax_sparsity(muAI{1},[1:2,(1:2)+nT-1]); 
        VST2 = varimax_sparsity(muAI{2},[1:2,(1:2)+nT-1]); 

        TVars = cellfun(@(x) sum(nanvar(x)),muAI); 
        TVARS{z}(i,:) = TVars; 
        
        TAv{z}(i,1:size(muAI{1},2)) = nanvar(muAI{1}*VST1);
        AAv{z}(i,1:size(Qboot.unique.C1,2)) = nanvar(muAI{1}*Qboot.unique.C1*VSC1);
        SAv{z}(i,1:size(Qboot.shared,2)) = nanvar(muAI{1}*Qboot.shared*VSS1);
        IAv{z}(i,1:size(Qboot.unique.C2,2)) = nanvar(muAI{1}*Qboot.unique.C2);
        
        TIv{z}(i,1:size(muAI{2},2)) = nanvar(muAI{2}*VST2);
        IIv{z}(i,1:size(Qboot.unique.C2,2)) = nanvar(muAI{2}*Qboot.unique.C2*VSC2);
        SIv{z}(i,1:size(Qboot.shared,2)) = nanvar(muAI{2}*Qboot.shared*VSS2);
        AIv{z}(i,1:size(Qboot.unique.C1,2)) = nanvar(muAI{2}*Qboot.unique.C1);
      
        %%%% FORCE 
        VSC1f = varimax_sparsity(muAIf{1}*Qboot.unique.C1,1:2); 
        VSC2f = varimax_sparsity(muAIf{2}*Qboot.unique.C2,1:2); 
        VSS1f = varimax_sparsity(muAIf{1}*Qboot.shared,1:2); 
        VSS2f = varimax_sparsity(muAIf{2}*Qboot.shared,1:2); 

        FVars = cellfun(@(x) sum(nanvar(x)),muAIf); 
        FVARS{z}(i,:) = FVars; 
        
        AAv_f{z}(i,1:size(Qboot.unique.C1,2)) = nanvar(muAIf{1}*Qboot.unique.C1*VSC1f);
        SAv_f{z}(i,1:size(Qboot.shared,2)) = nanvar(muAIf{1}*Qboot.shared*VSS1f);
        IAv_f{z}(i,1:size(Qboot.unique.C2,2)) = nanvar(muAIf{1}*Qboot.unique.C2);
        
        IIv_f{z}(i,1:size(Qboot.unique.C2,2)) = nanvar(muAIf{2}*Qboot.unique.C2*VSC2f); 
        SIv_f{z}(i,1:size(Qboot.shared,2)) = nanvar(muAIf{2}*Qboot.shared*VSS2f);
        AIv_f{z}(i,1:size(Qboot.unique.C1,2)) = nanvar(muAIf{2}*Qboot.unique.C1);
      
        
    end
end
%%
threshrange = [0.001,0.02]; 
threshs = linspace(threshrange(1),threshrange(2),1000); 
Trange = NaN(2,6,length(threshs)); 
Frange = NaN(2,4,length(threshs)); 
thresh_1p = find(threshs>=0.01,1,'first'); 

vImagery_during_Action = cellfun(@(x) nansum(x,2),IAv,'uni',0); 
vAction_during_Imagery = cellfun(@(x) nansum(x,2),AIv,'uni',0); 

VS = [AAv, SAv, IIv, SIv, TAv, TIv]; 
VSp(1,[1 2 5]) = cellfun(@(x) (x./repmat(TVARS{1}(:,1),1,size(x,2)))*100,VS(1,[1 2 5]),'uni',0); 
VSp(2,[1 2 5]) = cellfun(@(x) (x./repmat(TVARS{2}(:,1),1,size(x,2)))*100,VS(2,[1 2 5]),'uni',0); 
VSp(1,[3 4 6]) = cellfun(@(x) (x./repmat(TVARS{1}(:,2),1,size(x,2)))*100,VS(1,[3 4 6]),'uni',0); 
VSp(2,[3 4 6]) = cellfun(@(x) (x./repmat(TVARS{2}(:,2),1,size(x,2)))*100,VS(2,[3 4 6]),'uni',0); 

VSF = [AAv_f, SAv_f, IIv_f, SIv_f]; 
VSFp(1,[1 2]) = cellfun(@(x) (x./repmat(TVARS{1}(:,1),1,size(x,2)))*100,VSF(1,[1 2]),'uni',0); 
VSFp(2,[1 2]) = cellfun(@(x) (x./repmat(TVARS{2}(:,1),1,size(x,2)))*100,VSF(2,[1 2]),'uni',0); 
VSFp(1,[3 4]) = cellfun(@(x) (x./repmat(TVARS{1}(:,2),1,size(x,2)))*100,VSF(1,[3 4]),'uni',0); 
VSFp(2,[3 4]) = cellfun(@(x) (x./repmat(TVARS{2}(:,2),1,size(x,2)))*100,VSF(2,[3 4]),'uni',0); 

for q = 1:length(threshs)
    thresh = threshs(q); 
    
    VSover(1,[1 2 5]) = cellfun(@(x) (x./repmat(TVARS{1}(:,1),1,size(x,2))) > thresh ,VS(1,[1 2 5]),'uni',0); 
    VSover(2,[1 2 5]) = cellfun(@(x) (x./repmat(TVARS{2}(:,1),1,size(x,2))) > thresh ,VS(2,[1 2 5]),'uni',0); 
    VSover(1,[3 4 6]) = cellfun(@(x) (x./repmat(TVARS{1}(:,2),1,size(x,2))) > thresh ,VS(1,[3 4 6]),'uni',0); 
    VSover(2,[3 4 6]) = cellfun(@(x) (x./repmat(TVARS{2}(:,2),1,size(x,2))) > thresh ,VS(2,[3 4 6]),'uni',0); 

    VSFover(1,[1 2]) = cellfun(@(x) (x./repmat(TVARS{1}(:,1),1,size(x,2)))>thresh ,VSF(1,[1 2]),'uni',0); 
    VSFover(2,[1 2]) = cellfun(@(x) (x./repmat(TVARS{2}(:,1),1,size(x,2)))>thresh ,VSF(2,[1 2]),'uni',0); 
    VSFover(1,[3 4]) = cellfun(@(x) (x./repmat(TVARS{1}(:,2),1,size(x,2)))>thresh ,VSF(1,[3 4]),'uni',0); 
    VSFover(2,[3 4]) = cellfun(@(x) (x./repmat(TVARS{2}(:,2),1,size(x,2)))>thresh ,VSF(2,[3 4]),'uni',0); 

    Trange(:,:,q) = cellfun(@(x) sum(mean(x)==1), VSover);
    Frange(:,:,q) = cellfun(@(x) sum(mean(x)==1), VSFover); 
        
end
for z = 1:2
    DS{z,:}.Total = mat2cell( squeeze(Trange(z,:,:))', size(Trange,3), ones(size(Trange,2),1)); 
    DS{z,:}.Force = mat2cell( squeeze(Frange(z,:,:))', size(Frange,3), ones(size(Frange,2),1)); 
end

TDimensionalities = Trange(:,:,thresh_1p); 
FDimensionalities = Frange(:,:,thresh_1p); 


%% Stats
% for each bootstrap run, calculate difference in dimensionalities
threshvec = ones(nboot,1)*0.01; 
threshvecf = ones(nboot,1)*0.01; 
[VSo,VSFo] = deal(cell(2,6)); 
VSo(1,[1 2 5]) = cellfun(@(x) (x./repmat(TVARS{1}(:,1),1,size(x,2)))*100 ,VS(1,[1 2 5]),'uni',0); 
VSo(2,[1 2 5]) = cellfun(@(x) (x./repmat(TVARS{2}(:,1),1,size(x,2)))*100 ,VS(2,[1 2 5]),'uni',0); 
VSo(1,[3 4 6]) = cellfun(@(x) (x./repmat(TVARS{1}(:,2),1,size(x,2)))*100 ,VS(1,[3 4 6]),'uni',0); 
VSo(2,[3 4 6]) = cellfun(@(x) (x./repmat(TVARS{2}(:,2),1,size(x,2)))*100 ,VS(2,[3 4 6]),'uni',0); 

[VSover_rand,VSFover_rand] = deal(cell(2,6)); 
VSover_rand(1,[1 2 5]) = cellfun(@(x) sum((x./repmat(TVARS{1}(:,1),1,size(x,2))) > repmat(threshvec,1,size(x,2)),2),VS(1,[1 2 5]),'uni',0); 
VSover_rand(2,[1 2 5]) = cellfun(@(x) sum((x./repmat(TVARS{2}(:,1),1,size(x,2))) > repmat(threshvec,1,size(x,2)),2),VS(2,[1 2 5]),'uni',0); 
VSover_rand(1,[3 4 6]) = cellfun(@(x) sum((x./repmat(TVARS{1}(:,2),1,size(x,2))) > repmat(threshvec,1,size(x,2)),2),VS(1,[3 4 6]),'uni',0); 
VSover_rand(2,[3 4 6]) = cellfun(@(x) sum((x./repmat(TVARS{2}(:,2),1,size(x,2))) > repmat(threshvec,1,size(x,2)),2),VS(2,[3 4 6]),'uni',0); 

VSFover_rand(1,[1 2]) = cellfun(@(x) sum((x./repmat(TVARS{1}(:,1),1,size(x,2)))>repmat(threshvecf,1,size(x,2)),2),VSF(1,[1 2]),'uni',0); 
VSFover_rand(2,[1 2]) = cellfun(@(x) sum((x./repmat(TVARS{2}(:,1),1,size(x,2)))>repmat(threshvecf,1,size(x,2)),2),VSF(2,[1 2]),'uni',0); 
VSFover_rand(1,[3 4]) = cellfun(@(x) sum((x./repmat(TVARS{1}(:,2),1,size(x,2)))>repmat(threshvecf,1,size(x,2)),2),VSF(1,[3 4]),'uni',0); 
VSFover_rand(2,[3 4]) = cellfun(@(x) sum((x./repmat(TVARS{2}(:,2),1,size(x,2)))>repmat(threshvecf,1,size(x,2)),2),VSF(2,[3 4]),'uni',0); 

Aunique_v_Iunique = {VSover_rand{1,1}-VSover_rand{1,3}; VSover_rand{2,1}-VSover_rand{2,3}}; 
Aunique_v_Shared = {VSover_rand{1,1}-VSover_rand{1,2}; VSover_rand{2,1}-VSover_rand{2,2}}; 
Iunique_v_Shared = {VSover_rand{1,3}-VSover_rand{1,4}; VSover_rand{2,3}-VSover_rand{2,4}}; 
Usum_v_Atotal = {VSover_rand{1,5}-VSover_rand{1,1}-VSover_rand{1,2}; VSover_rand{2,5}-VSover_rand{2,1}-VSover_rand{2,2}};  
Usum_v_Itotal = {VSover_rand{1,6}-VSover_rand{1,3}-VSover_rand{1,4}; VSover_rand{2,6}-VSover_rand{2,3}-VSover_rand{2,4}};  

boot_p = @(x) min([mean(x<=0) mean(x>=0)]); 

p_vals_AvI = cellfun(boot_p,Aunique_v_Iunique);       
p_vals_AvS = cellfun(boot_p,Aunique_v_Shared);
p_vals_IvS = cellfun(boot_p,Iunique_v_Shared); 
p_vals_UvA = cellfun(boot_p,Usum_v_Atotal);
p_vals_UvI = cellfun(boot_p,Usum_v_Itotal); 
                  
Aunique_v_Shared_F = {VSFover_rand{1,1}-VSFover_rand{1,2}; VSFover_rand{2,1}-VSFover_rand{2,2}}; 
Iunique_v_Shared_F = {VSFover_rand{1,3}-VSFover_rand{1,4}; VSFover_rand{2,3}-VSFover_rand{2,4}}; 

p_vals_AvS_F = cellfun(boot_p,Aunique_v_Shared_F);
p_vals_IvS_F = cellfun(boot_p,Iunique_v_Shared_F); 


%% Plotting
figure('Position',[277 142 1348 821]); hold on; 
s1 = subplot(2,1,1); hold on; 
s2 = subplot(2,1,2); hold on; 
total_vex = [dA, dSA, dI, dSI, dTA, dTI]; 
force_vex = [dAf, dSAf, dIf, dSIf]; 
for z = 1:2            
    [vxs,vys] = deal(cell(1,7)); 
    axes(s1); 
    mrk = {'.','o',25,7}; 
    for i = 1:6

        [vys{i},vxs{i}] = linehist(0:0.5:20,DS{z}.Total{i},'normalize',1,'suppress_plot'); 
        vxs{i} = -vxs{i}/max(vxs{i})*0.2; 
    end
    
    % plotting
    zoff = (z-1.5)*.4; 
    patch(vxs{5}+1+zoff,vys{5},[0 0 0],'EdgeColor','none','FaceAlpha',0.5);
    patch(vxs{2}+2+zoff,vys{2},clrs(3,:),'EdgeColor','none','FaceAlpha',0.5);
    patch(vxs{1}+3+zoff,vys{1},clrs(1,:),'EdgeColor','none','FaceAlpha',0.5);
    patch(vxs{6}+5+zoff,vys{6},[.5 .5 .5],'EdgeColor','none','FaceAlpha',0.5);
    patch(vxs{4}+6+zoff,vys{4},clrs(3,:),'EdgeColor','none','FaceAlpha',0.5);
    patch(vxs{3}+7+zoff,vys{3},clrs(2,:),'EdgeColor','none','FaceAlpha',0.5);
    
    plot(1+zoff,TDimensionalities(z,5),mrk{z},'Color','k','MarkerSize',mrk{z+2}); 
    plot(2+zoff,TDimensionalities(z,2),mrk{z},'Color',clrs(3,:),'MarkerSize',mrk{z+2}); 
    plot(3+zoff,TDimensionalities(z,1),mrk{z},'Color',clrs(1,:),'MarkerSize',mrk{z+2}); 
    plot(5+zoff,TDimensionalities(z,6),mrk{z},'Color',[.5 .5 .5],'MarkerSize',mrk{z+2}); 
    plot(6+zoff,TDimensionalities(z,4),mrk{z},'Color',clrs(3,:),'MarkerSize',mrk{z+2}); 
    plot(7+zoff,TDimensionalities(z,3),mrk{z},'Color',clrs(2,:),'MarkerSize',mrk{z+2}); 

    ylabel('Total Dimensionality'); set(gca,'XTick',[1 2 3 5 6 7]); xticklabels({'action','shared','unique','imagery','shared','unique'}); 
    ylim([0 16]); xlim([0 8]); 
    
    
    [vxs,vys] = deal(cell(1,6)); 
    axes(s2); 
    for i = 1:4
        [vys{i},vxs{i}] = linehist(0:0.5:12,DS{z}.Force{i},'normalize',1,'suppress_plot'); 
        vxs{i} = -vxs{i}/max(vxs{i})*0.2; 
    end
    patch(vxs{2}+1+zoff,vys{2},clrs(3,:),'EdgeColor','none','FaceAlpha',0.5);
    patch(vxs{1}+2+zoff,vys{1},clrs(1,:),'EdgeColor','none','FaceAlpha',0.5);
    patch(vxs{4}+4+zoff,vys{4},clrs(3,:),'EdgeColor','none','FaceAlpha',0.5);
    patch(vxs{3}+5+zoff,vys{3},clrs(2,:),'EdgeColor','none','FaceAlpha',0.5);

    plot(1+zoff,FDimensionalities(z,2),mrk{z},'Color',clrs(3,:),'MarkerSize',mrk{z+2}); 
    plot(2+zoff,FDimensionalities(z,1),mrk{z},'Color',clrs(1,:),'MarkerSize',mrk{z+2}); 
    plot(4+zoff,FDimensionalities(z,4),mrk{z},'Color',clrs(3,:),'MarkerSize',mrk{z+2});
    plot(5+zoff,FDimensionalities(z,3),mrk{z},'Color',clrs(2,:),'MarkerSize',mrk{z+2}); 
    
    ylabel('Force Dimensionality'); set(gca,'XTick',[1 2 4 5]); xticklabels({'shared','unique','shared','unique'}); 
    ylim([0 8]); xlim([0 8]); 
end
clean_plot;

%%
EstDimns = TDimensionalities;

%%% Example fit from shared to unique action
X = ZOav{1}(:,1:EstDimns(1,1)); X2 = X-nanmean(X); 
Y = ZOSHav{1}(:,1:EstDimns(1,2)); Y2 = Y-nanmean(Y); 
g = ~any(isnan([X2,Y2]),2); 
b = Y2(g,:)\X2(g,:); 
rec = Y2*b; 
resid = X2-rec; 
Linds = 1:(size(X,1)/2); Hinds = (size(X,1)/2+1):size(X,1); 
figure; hold on; 
clrs_dim = nanmean(cat(3,clrs,ones(size(clrs))),3); 
for i = 1:size(rec,2)
    s1 = subplot(size(rec,2),3,i*3-2); hold on; if i==1; title('Action Unique'); end
    plot(X2(Linds,i),'Color',clrs_dim(1,:)); 
    plot(X2(Hinds,i),'Color',clrs(1,:)); yl = ylim; 
    
    s2 = subplot(size(rec,2),3,i*3-1); hold on; if i==1; title('Fit from Shared'); end
    plot(rec(Linds,i),'Color',clrs_dim(3,:)); 
    plot(rec(Hinds,i),'Color',clrs(3,:)); yl2 = ylim;  
    
    s3 = subplot(size(rec,2),3,i*3); hold on; if i==1; title('Residual'); end
    plot(resid(Linds,i),'Color',[.5 .5 .5]); 
    plot(resid(Hinds,i),'Color',[0 0 0]);  yl3 = ylim; 
    
    zl = [-.5 1]; %[min([yl,yl2,yl3]),max([yl,yl2,yl3])]; 
    set(s1,'YLim',zl); set(s2,'YLim',zl); set(s3,'YLim',zl); 
end
clean_plot; 
fprintf('%.1f%% unique fit from shared\n',(1 - sum(nanvar(resid))./sum(nanvar(X2)))*100); 
    
% Example fit from unique to shared
b = X2(g,:)\Y2(g,:); 
rec = X2*b; 
resid = Y2-rec; 
Linds = 1:(size(X,1)/2); Hinds = (size(X,1)/2+1):size(X,1); 
figure; hold on; 
for i = 1:size(rec,2)
    s1 = subplot(size(rec,2),3,i*3-2); hold on; if i==1; title('Shared'); end
    plot(Y2(Linds,i),'Color',clrs_dim(3,:)); 
    plot(Y2(Hinds,i),'Color',clrs(3,:)); yl = ylim; 
    
    s2 = subplot(size(rec,2),3,i*3-1); hold on; if i==1; title('Fit from Unique'); end
    plot(rec(Linds,i),'Color',clrs_dim(1,:)); 
    plot(rec(Hinds,i),'Color',clrs(1,:)); yl2 = ylim;  
    
    s3 = subplot(size(rec,2),3,i*3); hold on; if i==1; title('Residual'); end
    plot(resid(Linds,i),'Color',[.5 .5 .5]); 
    plot(resid(Hinds,i),'Color',[0 0 0]);  yl3 = ylim; 
    
    zl = [-1 2.5];%[min([yl,yl2,yl3]),max([yl,yl2,yl3])]; 
    set(s1,'YLim',zl); set(s2,'YLim',zl); set(s3,'YLim',zl); 
end
clean_plot;  
fprintf('%.1f%% shared fit from unique\n',(1 - sum(nanvar(resid))./sum(nanvar(Y2)))*100); 

%%
nboot = 1000; 
[pvar,pvex] = deal(cell(1,2)); 
cellavfun = @(x) nanmean(cell2mat(reshape(x,1,1,[])),3); 
for z = 1:2
    for i = 1:length(AIDZ{z})
       
        Xact = cellfun(@(Y) cellfun(@(x) x'*Subspaces{z}{i,1}(:,1:EstDimns(z,1)),Y,'uni',0),AIDZ{z}{i}.Action,'uni',0);
        Xact_sh = cellfun(@(Y) cellfun(@(x) x'*Subspaces{z}{i,3}(:,1:EstDimns(z,2)),Y,'uni',0),AIDZ{z}{i}.Action,'uni',0);
        
        Xim = cellfun(@(Y) cellfun(@(x) x'*Subspaces{z}{i,2}(:,1:EstDimns(z,3)),Y,'uni',0),AIDZ{z}{i}.Imagery,'uni',0);
        Xim_sh = cellfun(@(Y) cellfun(@(x) x'*Subspaces{z}{i,3}(:,1:EstDimns(z,4)),Y,'uni',0),AIDZ{z}{i}.Imagery,'uni',0);
        
        nb = cellfun(@(x) 1:length(x),[Xact,Xim],'uni',0); 
        booti = cell(1,numel(nb)); 
        for j = 1:numel(nb)
            [~,booti{j}] = bootstrp(nboot,[],nb{j}); 
            booti{j} = [booti{j} nb{j}']; % Tack on full (non boot)
        end
        
        for j = 1:(nboot+1)
            clc; fprintf('%d/%d (%d/%d) - %d/%d\n',z,2,i,length(AIDZ{z}),j,nboot); 

            Xaa = [cellavfun(Xact{1}(booti{1}(:,j))); cellavfun(Xact{2}(booti{2}(:,j)))]; 
            Xas = [cellavfun(Xact_sh{1}(booti{1}(:,j))); cellavfun(Xact_sh{2}(booti{2}(:,j)))]; 
           
            Xii = [cellavfun(Xim{1}(booti{3}(:,j))); cellavfun(Xim{2}(booti{4}(:,j)))]; 
            Xis = [cellavfun(Xim_sh{1}(booti{3}(:,j))); cellavfun(Xim_sh{2}(booti{4}(:,j)))]; 
            
            g = ~any(isnan([Xaa Xas Xii Xis]),2); 
            
            Xaa = Xaa-nanmean(Xaa); 
            Xas = Xas-nanmean(Xas); 
            Xii = Xii-nanmean(Xii); 
            Xis = Xis-nanmean(Xis); 
            
            baa = Xaa(g,:)\Xas(g,:); raa = Xaa*baa; 
            bas = Xas(g,:)\Xaa(g,:); ras = Xas*bas; 
            bii = Xii(g,:)\Xis(g,:); rii = Xii*bii; 
            bis = Xis(g,:)\Xii(g,:); ris = Xis*bis; 
            
            if j==(nboot+1)
                pvar{z}(i,1) = 1-sum(nanvar(ras-Xaa))./sum(nanvar(Xaa)); 
                pvar{z}(i,2) = 1-sum(nanvar(raa-Xas))./sum(nanvar(Xas)); 
                pvar{z}(i,3) = 1-sum(nanvar(ris-Xii))./sum(nanvar(Xii)); 
                pvar{z}(i,4) = 1-sum(nanvar(rii-Xis))./sum(nanvar(Xis)); 
            else
                pvex{z}{i,1}(j,:) = 1-sum(nanvar(ras-Xaa))./sum(nanvar(Xaa)); 
                pvex{z}{i,2}(j,:) = 1-sum(nanvar(raa-Xas))./sum(nanvar(Xas)); 
                pvex{z}{i,3}(j,:) = 1-sum(nanvar(ris-Xii))./sum(nanvar(Xii)); 
                pvex{z}{i,4}(j,:) = 1-sum(nanvar(rii-Xis))./sum(nanvar(Xis)); 
            end
            
        end
    end
            
end
%
pvar = cellfun(@(x) x*100,pvar,'uni',0); 
pvex_bnds = cellfun(@(y) cellfun(@(x) prctile(x,[2.5 97.5])*100,y,'uni',0),pvex,'uni',0); 
mrkr = {'.','o'}; 
figure; hold on; 
for z = 1:2
    for j = 1:size(pvex_bnds{z},1)
        plot(pvar{z}(j,1)*[1 1],pvex_bnds{z}{j,2},'Color',[.5 .5 .5]); 
        plot(pvex_bnds{z}{j,1},pvar{z}(j,2)*[1 1],'Color',[.5 .5 .5]); 

        plot(pvar{z}(j,3)*[1 1],pvex_bnds{z}{j,4},'Color',[.5 .5 .5]); 
        plot(pvex_bnds{z}{j,3},pvar{z}(j,4)*[1 1],'Color',[.5 .5 .5]); 
        
    end
    plot(pvar{z}(:,1),pvar{z}(:,2),mrkr{z},'Color',clrs(1,:)); 
    plot(pvar{z}(:,3),pvar{z}(:,4),mrkr{z},'Color',clrs(2,:)); 
end
plot([40 100],[40 100],'k:');
axis square; 
xlabel('% Unique subspace variance explained by shared'); 
ylabel('% Shared subspace variance explained by unique'); 
clean_plot; 
