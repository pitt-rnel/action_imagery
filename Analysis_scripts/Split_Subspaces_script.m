% Load ActionImageryData 
subjects = {'P2','P3'}; 
clr = lines(50); 
bline = 1:10; 
[ACAT,Zcat,Tcat,Rcat,MScat,ZsplitAction,ZsplitImagery,ZsplitShared,Zbad,...
    ZOav,ZCav,ZOSHav,ZCSHav,ZOCav,ZCOav,Fbad,Vbars,Subspaces,Q_I2A,Q_S2A,BIND,Qsub] = deal(cell(1,2)); 
for zz = 1:length(subjects)
    sub = subjects{zz}; 
    
    all_ai = cellfun(@(x) [x.Action.PC_combined,x.Imagery.PC_combined],ActionImageryData.(sub),'uni',0); 
    for i = 1:length(all_ai)
        if zz==1; bind = 1:10; else; bind = 11:20; end
        condav = cellfun(@(x) nanmean(cell2mat(reshape(x,1,1,[])),3),all_ai{i},'uni',0); 
        condav = cellfun(@(x) x - repmat(nanmean(x(:,bind),2),1,size(x,2)),condav,'uni',0); 
        
        ACAT{zz}{i,:} = cell2mat(condav)'; 
    end
    minlen = min(cellfun(@(x) size(x,2),ACAT{zz})); 
    ACAT{zz} = cellfun(@(x) x(:,1:minlen),ACAT{zz},'uni',0);
    BIND{zz} = bind; 
    
    %%% Do Procrustes alignment across sessions %%%%%%%%%%%%%%%%%%%%%%%%%%
    [Zcat{zz},Tcat{zz},Rcat{zz},MScat{zz}] = generalProcrustes(ACAT{zz},bind,0); 
    aind = 1:(size(MScat{zz},1)/2); iind = (size(MScat{zz},1)/2+1):size(MScat{zz},1); 
    
    %%% Do subspace splitting/identification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Qai] = DySO({MScat{zz}(aind,:);MScat{zz}(iind,:)},99,'do_plot',false); 
    pcsh = pca(MScat{zz}*Qai.shared); Qai.shared = Qai.shared*pcsh; 
    pcua = pca(MScat{zz}(aind,:)*Qai.unique.C1); Qai.unique.C1 = Qai.unique.C1*pcua; 
    pcui = pca(MScat{zz}(iind,:)*Qai.unique.C2); Qai.unique.C2 = Qai.unique.C2*pcui; 
    Qsub{zz} = Qai; 
    
    P_av.cond1.subspaces = cellfun(@(x) MScat{zz}(aind,:)*x,{Qai.unique.C1, Qai.unique.C2, Qai.shared},'uni',0);
    P_av.cond2.subspaces = cellfun(@(x) MScat{zz}(iind,:)*x,{Qai.unique.C1, Qai.unique.C2, Qai.shared},'uni',0);
    SSav = {Qai.unique.C1, Qai.unique.C2, Qai.shared};
    
    [Qa,Sa,va] = varimax_sparsity(P_av.cond1.subspaces{1},bline); 
    [Qi,Si,vi] = varimax_sparsity(P_av.cond2.subspaces{2},bline); 
    [Qs,Ss,vs] = varimax_sparsity([P_av.cond1.subspaces{3};P_av.cond2.subspaces{3}],bline); 
    SSav{1} = SSav{1}*Qa; 
    SSav{2} = SSav{2}*Qi; 
    SSav{3} = SSav{3}*Qs; 
    
    n_ai = min([sum(va>1), size(Si,2)]); % Dimensionality
    
    %%% Align imagery-unique subspace responses to action-unique %%%%%%%%%%
    Q_I2A{zz} = Orthogonal_Alignment(Si,Sa(:,1:n_ai),'Covariance'); 
    SSav{2} = SSav{2}*Q_I2A{zz}; 
    
    n_as = min([sum(va>1), size(Ss,2)]); 

    Subspaces{zz} = [cellfun(@(x) x.T*SSav{1},Tcat{zz},'uni',0),...
                     cellfun(@(x) x.T*SSav{2},Tcat{zz},'uni',0),...
                     cellfun(@(x) x.T*SSav{3},Tcat{zz},'uni',0)]; 
    
    ZsplitAction{zz} = [cellfun(@(x) x(aind,:)*SSav{1},Zcat{zz},'uni',0), ...
                   cellfun(@(x) x(iind,:)*SSav{1},Zcat{zz},'uni',0)];
               
    ZsplitImagery{zz} = [cellfun(@(x) x(aind,:)*SSav{2},Zcat{zz},'uni',0), ...
                   cellfun(@(x) x(iind,:)*SSav{2},Zcat{zz},'uni',0)];
               
    ZsplitShared{zz} = [cellfun(@(x) x(aind,:)*SSav{3},Zcat{zz},'uni',0), ...
                   cellfun(@(x) x(iind,:)*SSav{3},Zcat{zz},'uni',0)];

    Vbars{zz}.Action = cellfun(@(x) nanvar(x),[ZsplitAction{zz}(:,1),ZsplitImagery{zz}(:,1),ZsplitShared{zz}(:,1)],'uni',0);            
    Vbars{zz}.Imagery = cellfun(@(x) nanvar(x),[ZsplitAction{zz}(:,2),ZsplitImagery{zz}(:,2),ZsplitShared{zz}(:,2)],'uni',0);      

    ZOav{zz} = mean(cell2mat(reshape(ZsplitAction{zz}(:,1),1,1,[])),3);
    ZCav{zz} = mean(cell2mat(reshape(ZsplitImagery{zz}(:,2),1,1,[])),3);        
    ZOCav{zz} = mean(cell2mat(reshape(ZsplitAction{zz}(:,2),1,1,[])),3);       
    ZCOav{zz} = mean(cell2mat(reshape(ZsplitImagery{zz}(:,1),1,1,[])),3);
    ZOSHav{zz} = mean(cell2mat(reshape(ZsplitShared{zz}(:,1),1,1,[])),3);
    ZCSHav{zz} = mean(cell2mat(reshape(ZsplitShared{zz}(:,2),1,1,[])),3);
end

%% Procrustes alignment plots (Supp Fig 3)
numrand = 1000; 
totnull = cell(1,2); 
for z = 1:2
    newZj = cell(1,1,length(Zcat{z})); 
    randcorrs = cell(1,numrand); 
    for i = 1:numrand
        clc; fprintf('%d/%d (%d/%d)\n',z,2,i,numrand); 
        for j = 1:length(Zcat{z})
            [us,~,~] = svd(rand(size(Zcat{z}{j},2))*2-1,'econ'); 
            newZj{j} = ACAT{z}{j}*us; 
        end
        meanzj = nanmean(cell2mat(newZj),3); 
        [~,ord] = sortrows(nanvar(meanzj)','descend');
        meanzj = meanzj(ord,:); 
        g = ~any(isnan(meanzj),2); 
        
        for j = 1:length(Zcat{z})
            randcorrs{i}(:,j) = abs(diag(corr(newZj{j}(g,:),meanzj(g,:)))); 
        end
    end
    totnull{z} = cell2mat(randcorrs); 
    
end
bnd99_1 = prctile(totnull{1}(:),99); 
bnd99_2 = prctile(totnull{2}(:),99); 

figure; hold on; 
subplot(2,1,1); hold on;
cellfun(@(x) plot(x,'.-'),Rcat{1}); 
legend(cellfun(@(x) sprintf('session %d',x),num2cell(1:length(Rcat{1})),'uni',0)); 
xl = xlim; 
plot(xl,[1 1]*bnd99_1,'Color',[.5 .5 .5],'HandleVisibility','off'); 
ylim([0 1]); 
xlabel('dimensions'); ylabel('Correlation with mean');

subplot(2,1,2); hold on;
cellfun(@(x) plot(x,'.-'),Rcat{2}); 
legend(cellfun(@(x) sprintf('session %d',x),num2cell(1:length(Rcat{2})),'uni',0)); 
xl = xlim; 
plot(xl,[1 1]*bnd99_1,'Color',[.5 .5 .5],'HandleVisibility','off'); 
ylim([0 1]); 
xlabel('dimensions'); ylabel('Correlation with mean'); 
       
clean_plot; 

%% Subspace overlap and variance scatters (With shuffle control)
% [CONTROL] Attempt with shuffled
clr = lines(50); 
bline = 1:10; 
[ACATshuff,Zcatshuff,Tcatshuff,Rcatshuff,MScatshuff,ZsplitActionShuff,...
    ZsplitImageryShuff,ZsplitSharedShuff,VbarsShuff,SubspacesShuff,Q_I2AShuff,Q_S2AShuff,BINDShuff,QsubShuff] = deal(cell(1,2)); 
for zz = 1:2
    sub = subjects{zz}; 
    
    all_ai = cellfun(@(x) [x.Action.PC_combined,x.Imagery.PC_combined],ActionImageryData.(sub),'uni',0);
    
    % Shuffle conditions together
    for i = 1:length(all_ai)
        aicomb_L = [all_ai{i}{1};all_ai{i}{3}]; 
        aicomb_H = [all_ai{i}{2};all_ai{i}{4}]; 
        Lrand = randperm(length(aicomb_L)); 
        Hrand = randperm(length(aicomb_H)); 
        Lrand_ind1 = Lrand(1:length(all_ai{i}{1})); 
        Lrand_ind2 = Lrand((length(all_ai{i}{1})+1):end); 
        Hrand_ind1 = Hrand(1:length(all_ai{i}{2})); 
        Hrand_ind2 = Hrand((length(all_ai{i}{2})+1):end); 
        
        all_ai{i}{1} = aicomb_L(Lrand_ind1); 
        all_ai{i}{2} = aicomb_H(Hrand_ind1); 
        all_ai{i}{3} = aicomb_L(Lrand_ind2); 
        all_ai{i}{4} = aicomb_H(Hrand_ind2); 
    end
        
    for i = 1:length(all_ai)
        if zz==1; bind = 1:10; else; bind = 11:20; end
        condav = cellfun(@(x) nanmean(cell2mat(reshape(x,1,1,[])),3),all_ai{i},'uni',0); 
        condav = cellfun(@(x) x - repmat(nanmean(x(:,bind),2),1,size(x,2)),condav,'uni',0); 
        
        ACATshuff{zz}{i,:} = cell2mat(condav)'; 
    end
    minlen = min(cellfun(@(x) size(x,2),ACATshuff{zz})); 
    ACATshuff{zz} = cellfun(@(x) x(:,1:minlen),ACATshuff{zz},'uni',0);
    BINDShuff{zz} = bind; 
    [Zcatshuff{zz},Tcatshuff{zz},Rcatshuff{zz},MScatshuff{zz}] = generalProcrustes(ACATshuff{zz},bind,1); 
    aind = 1:(size(MScatshuff{zz},1)/2); iind = (size(MScatshuff{zz},1)/2+1):size(MScatshuff{zz},1); 

    [Qai] = DySO({MScatshuff{zz}(aind,:); MScatshuff{zz}(iind,:)},99,'do_plot',false); QsubShuff{zz} = Qai; 
    P_av.cond1.subspaces = cellfun(@(x) MScatshuff{zz}(aind,:)*x,{Qai.unique.C1, Qai.unique.C2, Qai.shared},'uni',0);
    P_av.cond2.subspaces = cellfun(@(x) MScatshuff{zz}(iind,:)*x,{Qai.unique.C1, Qai.unique.C2, Qai.shared},'uni',0);
    SSav = {Qai.unique.C1, Qai.unique.C2, Qai.shared};

    [Qa,Sa,va] = varimax_sparsity(P_av.cond1.subspaces{1},bline); 
    [Qi,Si,vi] = varimax_sparsity(P_av.cond2.subspaces{2},bline); 
    [Qs,Ss,vs] = varimax_sparsity([P_av.cond1.subspaces{3};P_av.cond2.subspaces{3}],bline); 
    SSav{1} = SSav{1}*Qa; 
    SSav{2} = SSav{2}*Qi; 
    SSav{3} = SSav{3}*Qs; 
    
    n_ai = min([sum(va>1), size(Si,2)]); 
    Q_I2AShuff{zz} = Orthogonal_Alignment(Si,Sa(:,1:n_ai),'Covariance'); 

    SSav{2} = SSav{2}*Q_I2AShuff{zz}; 
    

    ZsplitActionShuff{zz} = [cellfun(@(x) x(aind,:)*SSav{1},Zcatshuff{zz},'uni',0), ...
                   cellfun(@(x) x(iind,:)*SSav{1},Zcatshuff{zz},'uni',0)];
               
    ZsplitImageryShuff{zz} = [cellfun(@(x) x(aind,:)*SSav{2},Zcatshuff{zz},'uni',0), ...
                   cellfun(@(x) x(iind,:)*SSav{2},Zcatshuff{zz},'uni',0)];
               
    ZsplitSharedShuff{zz} = [cellfun(@(x) x(aind,:)*SSav{3},Zcatshuff{zz},'uni',0), ...
                   cellfun(@(x) x(iind,:)*SSav{3},Zcatshuff{zz},'uni',0)];

    VbarsShuff{zz}.Action = cellfun(@(x) nanvar(x),[ZsplitActionShuff{zz}(:,1),ZsplitImageryShuff{zz}(:,1),ZsplitSharedShuff{zz}(:,1)],'uni',0);            
    VbarsShuff{zz}.Imagery = cellfun(@(x) nanvar(x),[ZsplitActionShuff{zz}(:,2),ZsplitImageryShuff{zz}(:,2),ZsplitSharedShuff{zz}(:,2)],'uni',0);      
    

end



[sub_overlap,sub_overlap_shuff, sub_overlap_shuffs, sub_overlap_LH, overlap_p] = deal(cell(2,1)); 
bootnum = 1000; 
for z = 1:2
    sub = subjects{z}; 
    all_ai = cellfun(@(x) [x.Action.PC_combined,x.Imagery.PC_combined],ActionImageryData.(sub),'uni',0);
    
    for i = 1:size(ZsplitAction{z},1)
        
        O = [ZsplitAction{z}{i,1}, ZsplitImagery{z}{i,1}, ZsplitShared{z}{i,1}]; 
        C = [ZsplitAction{z}{i,2}, ZsplitImagery{z}{i,2}, ZsplitShared{z}{i,2}]; 

        pcO = pca(O); 
        pcC = pca(C); 

        vOO = nanvar(O*pcO)./sum(nanvar(O))*100; 
        vCO = nanvar(C*pcO)./sum(nanvar(C))*100; 
        vCC = nanvar(C*pcC)./sum(nanvar(C))*100; 
        vOC = nanvar(O*pcC)./sum(nanvar(O))*100; 
        
        dO = find(cumsum(vOO)>99,1,'first'); 
        dC = find(cumsum(vCC)>99,1,'first'); 
        
        sub_overlap{z}(i,1) = sum(vCO(1:dO))./sum(vCC(1:dO)); 
        sub_overlap{z}(i,2) = sum(vOC(1:dC))./sum(vOO(1:dC)); 
        
        %%%%%%%%%%%%%%%%%%%%
        aicomb_L = [all_ai{i}{1};all_ai{i}{3}]; % Low forces
        aicomb_H = [all_ai{i}{2};all_ai{i}{4}]; % High forces
        for k = 1:bootnum
            clc; fprintf('%d/%d (%d/%d) - %d/%d\n',z,2,i,size(ZsplitAction{z},1),k,bootnum); 
            Lrand = randperm(length(aicomb_L)); 
            Hrand = randperm(length(aicomb_H)); 
            Lrand_ind1 = Lrand(1:length(all_ai{i}{1})); 
            Lrand_ind2 = Lrand((length(all_ai{i}{1})+1):end); 
            Hrand_ind1 = Hrand(1:length(all_ai{i}{2})); 
            Hrand_ind2 = Hrand((length(all_ai{i}{2})+1):end); 

            all_ai{i}{1} = aicomb_L(Lrand_ind1); 
            all_ai{i}{2} = aicomb_H(Hrand_ind1); 
            all_ai{i}{3} = aicomb_L(Lrand_ind2); 
            all_ai{i}{4} = aicomb_H(Hrand_ind2); 
        
            Shuff1 = [nanmean(cell2mat(reshape(all_ai{i}{1},1,1,[])),3), ...
                      nanmean(cell2mat(reshape(all_ai{i}{2},1,1,[])),3)]'; 
                  
            Shuff2 = [nanmean(cell2mat(reshape(all_ai{i}{3},1,1,[])),3), ...
                      nanmean(cell2mat(reshape(all_ai{i}{4},1,1,[])),3)]'; 
        
            pcShuff1 = pca(Shuff1); 
            pcShuff2 = pca(Shuff2); 
            
            vS1S1 = nanvar(Shuff1*pcShuff1)./sum(nanvar(Shuff1))*100; 
            vS2S1 = nanvar(Shuff2*pcShuff1)./sum(nanvar(Shuff2))*100; 
            vS2S2 = nanvar(Shuff2*pcShuff2)./sum(nanvar(Shuff2))*100; 
            vS1S2 = nanvar(Shuff1*pcShuff2)./sum(nanvar(Shuff1))*100; 
            
            dshuff1 = sum(cumsum(vS1S1)<=99)+1; 
            dshuff2 = sum(cumsum(vS2S2)<=99)+1; 

            sub_overlap_shuffs{z}{i,1}(k,1) = sum(vS2S1(1:dshuff1))./sum(vS2S2(1:dshuff1)); 
            sub_overlap_shuffs{z}{i,1}(k,2) = sum(vS1S2(1:dshuff2))./sum(vS1S1(1:dshuff2)); 
            
        end
        
        overlap_p{z}(i,1) = sum(sub_overlap_shuffs{z}{i}(:) < nanmean(sub_overlap{z}(i,:))); 

        Shuff1 = [ZsplitActionShuff{z}{i,1}, ZsplitImageryShuff{z}{i,1}, ZsplitSharedShuff{z}{i,1}]; 
        Shuff2 = [ZsplitActionShuff{z}{i,2}, ZsplitImageryShuff{z}{i,2}, ZsplitSharedShuff{z}{i,2}]; 

        pcShuff1 = pca(Shuff1); 
        pcShuff2 = pca(Shuff2); 

        vS1S1 = nanvar(Shuff1*pcShuff1)./sum(nanvar(Shuff1))*100; 
        vS2S1 = nanvar(Shuff2*pcShuff1)./sum(nanvar(Shuff2))*100; 
        vS2S2 = nanvar(Shuff2*pcShuff2)./sum(nanvar(Shuff2))*100; 
        vS1S2 = nanvar(Shuff1*pcShuff2)./sum(nanvar(Shuff1))*100; 

        dshuff1 = find(cumsum(vS1S1)>99,1,'first'); 
        dshuff2 = find(cumsum(vS2S2)>99,1,'first'); 

        sub_overlap_shuff{z}(i,1) = sum(vS2S1(1:dshuff1))./sum(vS2S2(1:dshuff1)); 
        sub_overlap_shuff{z}(i,2) = sum(vS1S2(1:dshuff2))./sum(vS1S1(1:dshuff2));  
    end
end

%% Do Subspace Overlap and subspace variance plots (with shuffle control)

figure('Position',[127 528 1706 437]); hold on; subplot(1,3,1); hold on; title('Subspace Overlap');
sover = {sub_overlap{1}(:), sub_overlap{2}(:)}; 
[b1,b2] = deal(zeros(1,2)); 
[b1(1),b1(2)] = boot_bounds(1000,@nanmean,sover{1},2.5,97.5); 
[b2(1),b2(2)] = boot_bounds(1000,@nanmean,sover{2},2.5,97.5); 

bar(1,nanmean(sover{1}),'EdgeColor','none'); 
bar(2,nanmean(sover{2}),'EdgeColor','none'); 
plot([1 1],b1,'k','LineWidth',2); 
plot([2 2],b2,'k','LineWidth',2); 

sover_shuff = {reshape(cell2mat(sub_overlap_shuffs{1}),[],1), reshape(cell2mat(sub_overlap_shuffs{2}),[],1)}; 

b1 = prctile(sover_shuff{1},[2.5,97.5]); 
b2 = prctile(sover_shuff{2},[2.5,97.5]); 

bar(5,nanmean(sover_shuff{1}),'EdgeColor','none','FaceColor',[.5 .5 .5]); 
bar(6,nanmean(sover_shuff{2}),'EdgeColor','none','FaceColor',[.5 .5 .5]); 
plot([5 5],b1,'k','LineWidth',2); 
plot([6 6],b2,'k','LineWidth',2); 
ylim([0 1]); 

ylabel('Alignment Index'); 
set(gca,'XTick',[1.5 5.5]); 
xticklabels({'action/imagery', 'shuffled'}); 
clean_plot; 

mrkr = {'.','o';15,5}; 

clrs = lines(10); clrs(3,:) = []; 
s1 = subplot(1,3,2); hold on; title('actual');
s2 = subplot(1,3,3); hold on; title('shuffled'); 
mrk = {'.','o';15,5}; 
for zz = 1:2
    a_var = cellfun(@sum,Vbars{zz}.Action); 
    i_var = cellfun(@sum,Vbars{zz}.Imagery); 
    
    a_ve = a_var./repmat(nansum(a_var,2),1,size(a_var,2))*100; 
    i_ve = i_var./repmat(nansum(i_var,2),1,size(i_var,2))*100; 
    for i = 1:3
        plot(s1,a_ve(:,i),i_ve(:,i),mrk{1,zz},'MarkerSize',mrk{2,zz},'Color',clrs(i,:)); 
    end
    
    o_var = cellfun(@sum,VbarsShuff{zz}.Action); 
    e_var = cellfun(@sum,VbarsShuff{zz}.Imagery); 
    
    o_ve = o_var./repmat(nansum(o_var,2),1,size(o_var,2))*100; 
    e_ve = e_var./repmat(nansum(e_var,2),1,size(e_var,2))*100; 
    for i = 1:3
        plot(s2,o_ve(:,i),e_ve(:,i),mrk{1,zz},'MarkerSize',mrk{2,zz},'Color',clrs(i,:)); 
    end
end
plot(s1,[0 100],[0 100],':','Color',[.5 .5 .5]); 
plot(s2,[0 100],[0 100],':','Color',[.5 .5 .5]); 
set(s1,'XLim',[0 100],'YLim',[0 100]); axis(s1,'square'); 
set(s2,'XLim',[0 100],'YLim',[0 100]); axis(s2,'square'); 
xlabel(s1,'% Action Variance'); ylabel(s1,'% Imagery Variance'); 
xlabel(s2,'% Shuff 1 Variance'); ylabel(s2,'% Shuff 2 Variance'); 
clean_plot; 


%% Example split histogram
[vO_inO,vC_inO,vO_inC,vC_inC,vO_inQ,vC_inQ] = deal(cell(2,1)); 
for z = 1:2
    Fac = nanmean(cell2mat(reshape(Zcat{z},1,1,[])),3); 
    actind = 1:(size(Fac,1)/2); imind = (size(Fac,1)/2+1):size(Fac,1); 
    allO_facs = Fac(actind,:); 
    allC_facs = Fac(imind,:); 
    
    pco = pca(allO_facs); 
    vO_inO{z} = nanvar(allO_facs*pco)./nansum(nanvar(allO_facs))*100; 
    vC_inO{z} = nanvar(allC_facs*pco)./nansum(nanvar(allC_facs))*100; 
    
    pcc = pca(allC_facs); 
    vO_inC{z} = nanvar(allO_facs*pcc)./nansum(nanvar(allO_facs))*100; 
    vC_inC{z} = nanvar(allC_facs*pcc)./nansum(nanvar(allC_facs))*100; 
    
    Q_us = [Qsub{z}.unique.all{1}, Qsub{z}.shared]; 
    
    vO_inQ{z} = nanvar(allO_facs*Q_us)./nansum(nanvar(allO_facs))*100; 
    vC_inQ{z} = nanvar(allC_facs*Q_us)./nansum(nanvar(allC_facs))*100; 
end

figure('Position',[89 541 1701 304]); 
subplot(1,3,1); hold on; title('Action PCs'); 
bar([vO_inO{1};vC_inO{1}]',1,'EdgeColor','none'); 
xlabel('PC #'); 
ylabel('% variance'); 
legend({'action','imagery'}); 
subplot(1,3,2); hold on; title('Imagery PCs'); 
bar([vO_inC{1};vC_inC{1}]',1,'EdgeColor','none'); 
xlabel('PC #'); 

subplot(1,3,3); hold on; title('Split Subspaces'); 
bar([vO_inQ{1};vC_inQ{1}]',1,'EdgeColor','none'); 
xlabel('Dimension #'); 
clean_plot


%% Plot Action, Imagery, and Shared components
nancorr = @(x,y) corr(x(~any(isnan([x y]),2)),y(~any(isnan([x y]),2))); 
for zz = 1:2
    LiHis = {1:size(ZsplitAction{zz}{1},1)/2, (size(ZsplitAction{zz}{1},1)/2+1):size(ZsplitAction{zz}{1},1)}; 
    if z == 2; LiHis = cellfun(@(x) x(20:end),LiHis,'uni',0); end
    Li = LiHis{1}; Hi = LiHis{2}; 
    
    figure('Position',[409 287 922 675]); hold on; 
    n_ai = 5;%min([size(ZsplitAction{zz}{1},2),size(ZsplitImagery{zz}{1,2},2)]); 
    for i = 1:n_ai
        s1 = subplot(n_ai,4,i*4-3); hold on; if i==1; title(sprintf('P%d: Action',zz+1)); end
        cellfun(@(x) plot(x(Li,i),'Color',nanmean([1 1 1;clr(i,:)]),'LineWidth',1),ZsplitAction{zz}(:,1)); 
        cellfun(@(x) plot(x(Hi,i),'Color',clr(i,:),'LineWidth',1),ZsplitAction{zz}(:,1)); 
        plot(ZOav{zz}(Li,i),'Color',nanmean([1 1 1;clr(i,:)]),'LineWidth',2); 
        plot(ZOav{zz}(Hi,i),'Color',clr(i,:),'LineWidth',2)
        yl1 = ylim; 
        xl = xlim; 
        text(xl(end),mean(yl1),sprintf('%.2f',nancorr(ZOav{zz}(:,i),ZCav{zz}(:,i)))); 
        title(sprintf('P%d A_{A}: %.1f%%',zz+1,nanvar(ZOav{zz}(:,i))./sum(nanvar([ZOav{zz},ZOSHav{zz},ZCOav{zz}]))*100)); 
        
        s1b = subplot(n_ai,4,i*4-2); hold on; 
        cellfun(@(x) plot(x(Li,i),'Color',nanmean([1 1 1;clr(i,:)]),'LineWidth',1),ZsplitAction{zz}(:,2)); 
        cellfun(@(x) plot(x(Hi,i),'Color',clr(i,:),'LineWidth',1),ZsplitAction{zz}(:,2)); 
        plot(ZOCav{zz}(Li,i),'Color',nanmean([1 1 1;clr(i,:)]),'LineWidth',2); 
        plot(ZOCav{zz}(Hi,i),'Color',clr(i,:),'LineWidth',2)
        yl2 = ylim; zl = [min([yl1,yl2]), max([yl1,yl2])]; 
        ylim(zl); xlim(xl); set(s1,'YLim',zl);  
        title(sprintf('P%d A_{I}: %.1f%%',zz+1,nanvar(ZOCav{zz}(:,i))./sum(nanvar([ZCav{zz},ZCSHav{zz},ZCOav{zz}]))*100)); 
        
        s2 = subplot(n_ai,4,i*4); hold on; 
        cellfun(@(x) plot(x(Li,i),'Color',nanmean([1 1 1;clr(i,:)]),'LineWidth',1),ZsplitImagery{zz}(:,2)); 
        cellfun(@(x) plot(x(Hi,i),'Color',clr(i,:),'LineWidth',1),ZsplitImagery{zz}(:,2)); 
        plot(ZCav{zz}(Li,i),'Color',nanmean([1 1 1;clr(i,:)]),'LineWidth',2); 
        plot(ZCav{zz}(Hi,i),'Color',clr(i,:),'LineWidth',2)
        yl3 = ylim; ylim(yl3); 
        title(sprintf('P%d I_{I}: %.1f%%',zz+1,nanvar(ZCav{zz}(:,i))./sum(nanvar([ZCav{zz},ZCSHav{zz},ZOCav{zz}]))*100)); 
        
        s2b = subplot(n_ai,4,i*4-1); hold on; 
        cellfun(@(x) plot(x(Li,i),'Color',nanmean([1 1 1;clr(i,:)]),'LineWidth',1),ZsplitImagery{zz}(:,1)); 
        cellfun(@(x) plot(x(Hi,i),'Color',clr(i,:),'LineWidth',1),ZsplitImagery{zz}(:,1)); 
        plot(ZCOav{zz}(Li,i),'Color',nanmean([1 1 1;clr(i,:)]),'LineWidth',2); 
        plot(ZCOav{zz}(Hi,i),'Color',clr(i,:),'LineWidth',2)
        yl4 = ylim; zl = [min([yl3,yl4]), max([yl3,yl4])]; 
        ylim(zl); xlim(xl);  
        set(s2,'YLim',zl);  
        title(sprintf('P%d I_{A}: %.1f%%',zz,nanvar(ZCOav{zz}(:,i))./sum(nanvar([ZOav{zz},ZOSHav{zz},ZCOav{zz}]))*100)); 
%         zl = [min([yl1,yl2]),max([yl1 yl2])]; set(s1,'YLim',zl); set(s2,'YLim',zl); 
        
    end
    clean_plot;

    figure('Position',[409 287 922 675]); hold on; 
    na = 5;%size(ZsplitShared{zz}{1},2);
    for i = 1:na
        s1 = subplot(na,2,i*2-1); hold on; %if i==1; title(sprintf('P%d: Shared',zz)); end
        cellfun(@(x) plot(x(Li,i),'Color',nanmean([1 1 1;clr(i,:)]),'LineWidth',1),ZsplitShared{zz}(:,1)); 
        cellfun(@(x) plot(x(Hi,i),'Color',clr(i,:),'LineWidth',1),ZsplitShared{zz}(:,1)); 
        plot(ZOSHav{zz}(Li,i),'Color',nanmean([1 1 1;clr(i,:)]),'LineWidth',2); 
        plot(ZOSHav{zz}(Hi,i),'Color',clr(i,:),'LineWidth',2)
        yl1 = ylim; xl = xlim; 
        text(xl(end),mean(yl1),sprintf('%.2f',nancorr(ZOSHav{zz}(:,i),ZCSHav{zz}(:,i)))); 
        title(sprintf('P%d S_{A}: %.1f%%',zz+1,nanvar(ZOSHav{zz}(:,i))./sum(nanvar([ZOav{zz},ZOSHav{zz},ZCOav{zz}]))*100)); 
        
        s2 = subplot(na,2,i*2); hold on; 
        cellfun(@(x) plot(x(Li,i),'--','Color',nanmean([1 1 1;clr(i,:)]),'LineWidth',1),ZsplitShared{zz}(:,2)); 
        cellfun(@(x) plot(x(Hi,i),'--','Color',clr(i,:),'LineWidth',1),ZsplitShared{zz}(:,2)); 
        plot(ZCSHav{zz}(Li,i),'Color',nanmean([1 1 1;clr(i,:)]),'LineWidth',2); 
        plot(ZCSHav{zz}(Hi,i),'Color',clr(i,:),'LineWidth',2)
        yl2 = ylim; 
        zl = [min([yl1,yl2]),max([yl1 yl2])]; set(s1,'YLim',zl); set(s2,'YLim',zl); 
        title(sprintf('P%d S_{I}: %.1f%%',zz+1,nanvar(ZCSHav{zz}(:,i))./sum(nanvar([ZCav{zz},ZCSHav{zz},ZOCav{zz}]))*100)); 
    end
    clean_plot;
end