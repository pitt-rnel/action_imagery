nboot = 1000; 
cost_method = 'VAF'; 
EstDimns = round(TDimensionalities);
Rdyn = cell(1,2); 
for z = 1:2
    LiHis = {1:size(ZsplitAction{z}{1},1)/2, (size(ZsplitAction{z}{1},1)/2+1):size(ZsplitAction{z}{1},1)}; 
    if z == 2; LiHis = cellfun(@(x) x(11:end),LiHis,'uni',0); end
    lh = [LiHis{1},LiHis{2}]'; 
    
    [AID] = deal(cell(length(ActionImageryData.(sprintf('P%d',z+1))),1)); 
    pcn = size(Subspaces{z}{1},1); 
    for i = 1:length(AID)
        for j = 1:2
            AID{i}.Action{1,j} = cellfun(@(x) x(1:pcn,LiHis{1}),ActionImageryData.(sprintf('P%d',z+1)){i}.Action.PC_combined{j},'uni',0); 
            AID{i}.Imagery{1,j} = cellfun(@(x) x(1:pcn,LiHis{1}),ActionImageryData.(sprintf('P%d',z+1)){i}.Imagery.PC_combined{j},'uni',0); 
        end
    end

    ActIm = cellfun(@(x) [x.Action,x.Imagery],AID,'uni',0); ActIm = vertcat(ActIm{:}); 
    bA = cell(size(ActIm)); 
    for i = 1:numel(bA)
        [~,bA{i}] = bootstrp(nboot,[],ActIm{i}); 
    end
    nT = size(AID{1}.Action{1}{1},2); Aind = 1:(nT*2); Iind = (nT*2+1):4*nT; 
    
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
        
        ZOboot = muAI{1}*Qboot.unique.C1*VSC1; 
        ZOSHboot = muAI{1}*Qboot.shared*VSS1; 
        ZIboot = muAI{2}*Qboot.unique.C2*VSC2; 
        ZISHboot = muAI{2}*Qboot.shared*VSS2; 
        
        %%%%
        [~,~,ro] = ReverseFit(ZOboot(:,1:EstDimns(z,1)),ZOSHboot(:,1:EstDimns(z,2)),cost_method);
        [~,~,ros] = ReverseFit(ZOSHboot(:,1:EstDimns(z,2)),ZOboot(:,1:EstDimns(z,1)),cost_method);
        [~,~,ri] = ReverseFit(ZIboot(:,1:EstDimns(z,3)),ZISHboot(:,1:EstDimns(z,4)),cost_method);
        [~,~,ris] = ReverseFit(ZISHboot(:,1:EstDimns(z,4)),ZIboot(:,1:EstDimns(z,3)),cost_method);
        
        Rdyn{z}(i,:) = {ro, ros, ri, ris}; 
    end
end

%%
figure('Position',[239 510 1366 291]); hold on; 
s1 = subplot(1,4,1); hold on; title('Action-unique'); 
plot(fliplr(-(1:length(RS{1}.A))),RS{1}.A,'.-','Color',clrs(1,:)); 
plot(fliplr(-(1:length(RS{2}.A))),RS{2}.A,'.-','Color',clrs(2,:)); 
ylim([0 1]); 
xl1 = xlim; 
ylabel('Max Corr'); 
xlabel('Component # (sorted)'); 

s1b = subplot(1,4,2); hold on; title('Action Shared'); 
plot(fliplr(-(1:length(RS{1}.SA))),RS{1}.SA,'.-','Color',clrs(1,:)); 
plot(fliplr(-(1:length(RS{2}.SA))),RS{2}.SA,'.-','Color',clrs(2,:)); 
ylim([0 1]); 
xl1b = xlim; 

s2 = subplot(1,4,3); hold on; title('Imagery-unique'); 
plot(fliplr(-(1:length(RS{1}.I))),RS{1}.I,'.-','Color',clrs(1,:)); 
plot(fliplr(-(1:length(RS{2}.I))),RS{2}.I,'.-','Color',clrs(2,:)); 
ylim([0 1]); 
xl2 = xlim; 

s2b = subplot(1,4,4); hold on; title('Imagery Shared'); 
plot(fliplr(-(1:length(RS{1}.SI))),RS{1}.SI,'.-','Color',clrs(1,:)); 
plot(fliplr(-(1:length(RS{2}.SI))),RS{2}.SI,'.-','Color',clrs(2,:)); 
ylim([0 1]); 
xl2b = xlim; 

zl = [min([xl1(1),xl1b(1),xl2(1),xl2b(1)]),0]; 
set(s1,'XLim',zl); 
set(s1b,'XLim',zl); 
set(s2,'XLim',zl); 
set(s2b,'XLim',zl); 

clean_plot; 
