%% FORCE-dependent activity
nU = 4; nS = 3; 
[LH_A,LH_I,LH_SA,LH_SI,VAF_force] = deal(cell(2,1)); 
for z = 1:2
    LiHis = {1:size(ZsplitAction{z}{1},1)/2, (size(ZsplitAction{z}{1},1)/2+1):size(ZsplitAction{z}{1},1)}; 
    if z == 2; LiHis = cellfun(@(x) x(20:end),LiHis,'uni',0); end
    Ad = ZOav{z}(LiHis{2},:)-ZOav{z}(LiHis{1},:); 
    Id = ZCav{z}(LiHis{2},:)-ZCav{z}(LiHis{1},:);
    SAd = ZOSHav{z}(LiHis{2},:)-ZOSHav{z}(LiHis{1},:);
    SId = ZCSHav{z}(LiHis{2},:)-ZCSHav{z}(LiHis{1},:); 
    ATd = [Ad, SAd]; 
    ITd = [Id, SId]; 
    
    [QdAT] = varimax_sparsity(ATd,1:10); 
    [QdIT] = varimax_sparsity(ITd,1:10); 
    
    fdim_thresh = sum(nanvar([ZOav{z},ZOSHav{z},ZCOav{z}]))*tdim_thresh/sum(nanvar([Ad,SAd])); 
    
    [pcA_F,~] = varimax_sparsity(Ad,1:10); 
    [pcSA_F,~] = varimax_sparsity(SAd,1:10); 
    pcI_F = pcA_F; 
    pcSI_F = pcSA_F; 
    for i = 1:size(ZsplitAction{z},1)
        LH_A{z}(i,:) = {ZsplitAction{z}{i,1}(LiHis{1},:)*pcA_F, ZsplitAction{z}{i,1}(LiHis{2},:)*pcA_F};
        LH_I{z}(i,:) = {ZsplitImagery{z}{i,2}(LiHis{1},:)*pcI_F, ZsplitImagery{z}{i,2}(LiHis{2},:)*pcI_F};
        LH_SA{z}(i,:) = {ZsplitShared{z}{i,1}(LiHis{1},:)*pcSA_F, ZsplitShared{z}{i,1}(LiHis{2},:)*pcSA_F};
        LH_SI{z}(i,:) = {ZsplitShared{z}{i,2}(LiHis{1},:)*pcSI_F, ZsplitShared{z}{i,2}(LiHis{2},:)*pcSI_F};
    end
    VAF_force{z}.action_unique = nanvar(Ad*pcA_F)./sum(nanvar(ATd))*100; 
    VAF_force{z}.action_shared = nanvar(SAd*pcSA_F)./sum(nanvar(ATd))*100; 
    VAF_force{z}.imagery_unique = nanvar(Id*pcI_F)./sum(nanvar(ITd))*100; 
    VAF_force{z}.imagery_shared = nanvar(SId*pcSI_F)./sum(nanvar(ITd))*100; 
end


LH_Ad = LH_A; LH_Id = LH_I; LH_SAd = LH_SA; LH_SId = LH_SI; 
for z = 1:2
    for i = 1:size(LH_A{z},1)
        LH_Ad{z}(i,:) = cellfun(@(x) x-nanmean(cell2mat(reshape(LH_Ad{z}(i,:),1,1,[])),3),LH_Ad{z}(i,:),'uni',0); 
        LH_Id{z}(i,:) = cellfun(@(x) x-nanmean(cell2mat(reshape(LH_Id{z}(i,:),1,1,[])),3),LH_Id{z}(i,:),'uni',0); 
        LH_SAd{z}(i,:) = cellfun(@(x) x-nanmean(cell2mat(reshape(LH_SAd{z}(i,:),1,1,[])),3),LH_SAd{z}(i,:),'uni',0); 
        LH_SId{z}(i,:) = cellfun(@(x) x-nanmean(cell2mat(reshape(LH_SId{z}(i,:),1,1,[])),3),LH_SId{z}(i,:),'uni',0); 
    end
end

for z = 1:2
    figure; hold on; 
    t = (1:size(LH_A{z}{1},1))'; 
    sp = cell(1,2); 
%     subplot(4,2,z); hold on; title('Action'); 
    yls = []; 
    for i = 1:nU
        sp{1}{i} = subplot(nU+nS,2,i*2-1); hold on; if i==1; title(sprintf('P%d',z+1)); end 
        cellfun(@(x) patchline(t,x(:,i),'EdgeColor',clrs_dim(1,:),'EdgeAlpha',0.2),LH_Ad{z}(:,1)); 
        plot(nanmean(cell2mat(cellfun(@(x) x(:,i),LH_Ad{z}(:,1)','uni',0)),2),'Color',clrs_dim(1,:),'LineWidth',1); 
        cellfun(@(x) patchline(t,x(:,i),'EdgeColor',clrs(1,:),'EdgeAlpha',0.2),LH_Ad{z}(:,2)); 
        plot(nanmean(cell2mat(cellfun(@(x) x(:,i),LH_Ad{z}(:,2)','uni',0)),2),'Color',clrs(1,:),'LineWidth',1); 
        yls = [yls ylim]; 
    end
    for i = 1:nU; set(sp{1}{i},'YLim',[min(yls),max(yls)]); text(sp{1}{i},1,0,sprintf('%.1f%%',VAF_force{z}.action_unique(i))); end
    yls = [] ;
    for i = 1:nU
        sp{2}{i} = subplot(nU+nS,2,i*2); hold on; 
        cellfun(@(x) patchline(t,x(:,i),'EdgeColor',clrs_dim(2,:),'EdgeAlpha',0.2),LH_Id{z}(:,1)); 
        plot(nanmean(cell2mat(cellfun(@(x) x(:,i),LH_Id{z}(:,1)','uni',0)),2),'Color',clrs_dim(2,:),'LineWidth',1); 
        cellfun(@(x) patchline(t,x(:,i),'EdgeColor',clrs(2,:),'EdgeAlpha',0.2),LH_Id{z}(:,2)); 
        plot(nanmean(cell2mat(cellfun(@(x) x(:,i),LH_Id{z}(:,2)','uni',0)),2),'Color',clrs(2,:),'LineWidth',1); 
        yls = [yls ylim]; 
    end
    for i = 1:nU; set(sp{2}{i},'YLim',[min(yls),max(yls)]); text(sp{2}{i},1,0,sprintf('%.1f%%',VAF_force{z}.imagery_unique(i))); end

    yls = []; 
    for i = 1:nS
        sp{1}{i} = subplot(nU+nS,2,i*2-1+2*nU); hold on; 
        cellfun(@(x) patchline(t,x(:,i),'EdgeColor',clrs_dim(3,:),'EdgeAlpha',0.2),LH_SAd{z}(:,1)); 
        plot(nanmean(cell2mat(cellfun(@(x) x(:,i),LH_SAd{z}(:,1)','uni',0)),2),'Color',clrs_dim(3,:),'LineWidth',1); 
        cellfun(@(x) patchline(t,x(:,i),'EdgeColor',clrs(3,:),'EdgeAlpha',0.2),LH_SAd{z}(:,2)); 
        plot(nanmean(cell2mat(cellfun(@(x) x(:,i),LH_SAd{z}(:,2)','uni',0)),2),'Color',clrs(3,:),'LineWidth',1); 
        yls = [yls ylim]; 
    end
    for i = 1:nS; set(sp{1}{i},'YLim',[min(yls),max(yls)]); text(sp{1}{i},1,0,sprintf('%.1f%%',VAF_force{z}.action_shared(i))); end
    for i = 1:nS
        sp{2}{i} = subplot(nU+nS,2,i*2+2*nU); hold on; 
        cellfun(@(x) patchline(t,x(:,i),'EdgeColor',clrs_dim(3,:),'EdgeAlpha',0.2),LH_SId{z}(:,1)); 
        plot(nanmean(cell2mat(cellfun(@(x) x(:,i),LH_SId{z}(:,1)','uni',0)),2),'Color',clrs_dim(3,:),'LineWidth',1); 
        cellfun(@(x) patchline(t,x(:,i),'EdgeColor',clrs(3,:),'EdgeAlpha',0.2),LH_SId{z}(:,2)); 
        plot(nanmean(cell2mat(cellfun(@(x) x(:,i),LH_SId{z}(:,2)','uni',0)),2),'Color',clrs(3,:),'LineWidth',1); 
        yls = [yls ylim]; 
    end
    for i = 1:nS; set(sp{2}{i},'YLim',[min(yls),max(yls)]); text(sp{2}{i},1,0,sprintf('%.1f%%',VAF_force{z}.imagery_shared(i))); end
    clean_plot;
end
