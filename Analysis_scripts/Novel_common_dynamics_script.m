%% Identify novel->common components for Action, Imagery, and Shared subspaces
clrs = lines(10); clrs(3,:) = []; 
clrs_dim = nanmean(cat(3,clrs,ones(size(clrs))),3); 
cost_method = 'VAF'; 
[UY,UY_S,UYI,UYI_S,RS,VE,resvar,Bfits] = deal(cell(1,2)); 
EstDimns = round(TDimensionalities);

for z = 1:2
    [Uy,Bx,r,vaf,Zy,Zx,cost] = ReverseFit(ZOav{z}(:,1:EstDimns(z,1)),ZOSHav{z}(:,1:EstDimns(z,2)),cost_method);
    UY{z} = Uy; 
    actx = ZOav{z}(:,1:EstDimns(z,1))*Uy; actx = actx - nanmean(actx); 
    fity = ZOSHav{z}(:,1:EstDimns(z,2))*Bx; fity = fity-nanmean(fity); 
    res = actx-fity; 
    resvar{z}.action_unique = sum(nanvar(res))./sum(nanvar(actx))*100; 
    
    for i = 1:size(Bx,2); Bx(:,i) = Bx(:,i)./norm(Bx(:,i)); end
    Bfits{z}.StoA = Bx; 
    
    [Uys,Bxs,rs,vafs,~,~,costs] = ReverseFit(ZOSHav{z}(:,1:EstDimns(z,2)),ZOav{z}(:,1:EstDimns(z,1)),cost_method);
    UY_S{z} = Uys; 
    actx = ZOSHav{z}(:,1:EstDimns(z,2))*Uys; actx = actx - nanmean(actx); 
    fity = ZOav{z}(:,1:EstDimns(z,1))*Bxs; fity = fity-nanmean(fity); 
    res = actx-fity; 
    resvar{z}.action_shared = sum(nanvar(res))./sum(nanvar(actx))*100; 
    for i = 1:size(Bxs,2); Bxs(:,i) = Bxs(:,i)./norm(Bxs(:,i)); end
    Bfits{z}.AtoS = Bxs; 
    
    dmax = size(Uy,2); 
    figure('Position',[250 72 700 880]); hold on; 
    
    RS{z}.A = r'; 
    RS{z}.SA = rs';
    VE{z}.A = {nanvar(ZOav{z}(:,1:EstDimns(z,1))*Uy)./sum(nanvar(ZOav{z}(:,1:EstDimns(z,1))))*100,...
               nanvar(ZOSHav{z}(:,1:EstDimns(z,2))*Bx)./sum(nanvar(ZOSHav{z}(:,1:EstDimns(z,2))))*100}; 
    VE{z}.SA = {nanvar(ZOSHav{z}(:,1:EstDimns(z,2))*Uys)./sum(nanvar(ZOSHav{z}(:,1:EstDimns(z,2))))*100,...
               nanvar(ZOav{z}(:,1:EstDimns(z,1))*Bxs)./sum(nanvar(ZOav{z}(:,1:EstDimns(z,1))))*100}; 

    sps = cell(dmax,2); 
    for i = 1:dmax
        sps{i,1} = subplot(dmax+1,2,i*2-1); hold on; if i==1; title(sprintf('P%d',z+1)); end
        plot(ZOav{z}(Li,1:EstDimns(z,1))*Uy(:,i),'Color',clrs_dim(1,:));
        plot(ZOav{z}(Hi,1:EstDimns(z,1))*Uy(:,i),'Color',clrs(1,:));

        sps{i,2} = subplot(dmax+1,2,i*2); hold on; title(sprintf('r = %.2f',r(i))); 
        plot(ZOSHav{z}(Li,1:EstDimns(z,2))*Bx(:,i),'Color',clrs_dim(3,:));
        plot(ZOSHav{z}(Hi,1:EstDimns(z,2))*Bx(:,i),'Color',clrs(3,:));

    end

    for i = 1:size(ZsplitAction{z},1)
        sess_a = ZsplitAction{z}{i,1}(:,1:EstDimns(z,1))*Uy; 
        sess_s = ZsplitShared{z}{i,1}(:,1:EstDimns(z,2))*Bx; 
        fig = gcf; 

        for j = 1:size(sps,1)
            fig.CurrentAxes = sps{j,1}; 
            patchline(1:length(Li),sess_a(Li,j),'EdgeColor',clrs_dim(1,:),'EdgeAlpha',0.25); 
            patchline(1:length(Hi),sess_a(Hi,j),'EdgeColor',clrs(1,:),'EdgeAlpha',0.25); 

            fig.CurrentAxes = sps{j,2}; 
            patchline(1:length(Li),sess_s(Li,j),'EdgeColor',clrs_dim(3,:),'EdgeAlpha',0.25); 
            patchline(1:length(Hi),sess_s(Hi,j),'EdgeColor',clrs(3,:),'EdgeAlpha',0.25); 
        end
    end
    clean_plot
    
    dmax = size(Uys,2); 
    figure('Position',[250 72 700 880]); hold on; 
    
    sps = cell(dmax,2); 
    for i = 1:dmax
        sps{i,2} = subplot(dmax+1,2,i*2-1); hold on; title(sprintf('r = %.2f',rs(i))); 
        plot(ZOav{z}(Li,1:EstDimns(z,1))*Bxs(:,i),'Color',clrs_dim(1,:));
        plot(ZOav{z}(Hi,1:EstDimns(z,1))*Bxs(:,i),'Color',clrs(1,:));

        sps{i,1} = subplot(dmax+1,2,i*2); hold on; if i==1; title(sprintf('P%d',z+1)); end
        plot(ZOSHav{z}(Li,1:EstDimns(z,2))*Uys(:,i),'Color',clrs_dim(3,:));
        plot(ZOSHav{z}(Hi,1:EstDimns(z,2))*Uys(:,i),'Color',clrs(3,:));

    end

    for i = 1:size(ZsplitAction{z},1)
        sess_a = ZsplitAction{z}{i,1}(:,1:EstDimns(z,1))*Bxs; 
        sess_s = ZsplitShared{z}{i,1}(:,1:EstDimns(z,2))*Uys; 
        fig = gcf; 

        for j = 1:size(sps,1)
            fig.CurrentAxes = sps{j,2}; 
            patchline(1:length(Li),sess_a(Li,j),'EdgeColor',clrs_dim(1,:),'EdgeAlpha',0.25); 
            patchline(1:length(Hi),sess_a(Hi,j),'EdgeColor',clrs(1,:),'EdgeAlpha',0.25); 

            fig.CurrentAxes = sps{j,1}; 
            patchline(1:length(Li),sess_s(Li,j),'EdgeColor',clrs_dim(3,:),'EdgeAlpha',0.25); 
            patchline(1:length(Hi),sess_s(Hi,j),'EdgeColor',clrs(3,:),'EdgeAlpha',0.25); 
        end
    end
    clean_plot
end

% Imagery
for z = 1:2
    [Uy,Bx,r,vaf,Zy,Zx,cost] = ReverseFit(ZCav{z}(:,1:EstDimns(z,3)),ZCSHav{z}(:,1:EstDimns(z,4)),cost_method);
    UYI{z} = Uy; 
    actx = ZCav{z}(:,1:EstDimns(z,3))*Uy; actx = actx - nanmean(actx); 
    fity = ZCSHav{z}(:,1:EstDimns(z,4))*Bx; fity = fity-nanmean(fity); 
    res = actx-fity; 
    resvar{z}.imagery_unique = sum(nanvar(res))./sum(nanvar(actx))*100; 
%     r = cost; 
    for i = 1:size(Bx,2); Bx(:,i) = Bx(:,i)./norm(Bx(:,i)); end
    
    [Uys,Bxs,rs,vafs,~,~,costs] = ReverseFit(ZCSHav{z}(:,1:EstDimns(z,4)),ZCav{z}(:,1:EstDimns(z,3)),cost_method);
    UYI_S{z} =Uys; 
    actx = ZCSHav{z}(:,1:EstDimns(z,4))*Uys; actx = actx - nanmean(actx); 
    fity = ZCav{z}(:,1:EstDimns(z,3))*Bxs; fity = fity-nanmean(fity); 
    res = actx-fity; 
    resvar{z}.imagery_shared = sum(nanvar(res))./sum(nanvar(actx))*100; 
%     rs = costs; 
    for i = 1:size(Bxs,2); Bxs(:,i) = Bxs(:,i)./norm(Bxs(:,i)); end
    
    dmax = size(Uy,2); 
    figure('Position',[250 72 700 880]); hold on; 
    
    RS{z}.I = r'; 
    RS{z}.SI = rs';
    VE{z}.I = {nanvar(ZCav{z}(:,1:EstDimns(z,3))*Uy)./sum(nanvar(ZCav{z}(:,1:EstDimns(z,3))))*100,...
               nanvar(ZCSHav{z}(:,1:EstDimns(z,4))*Bx)./sum(nanvar(ZCSHav{z}(:,1:EstDimns(z,4))))*100}; 
    VE{z}.SI = {nanvar(ZCSHav{z}(:,1:EstDimns(z,4))*Uys)./sum(nanvar(ZCSHav{z}(:,1:EstDimns(z,4))))*100,...
               nanvar(ZCav{z}(:,1:EstDimns(z,3))*Bxs)./sum(nanvar(ZCav{z}(:,1:EstDimns(z,3))))*100}; 

    sps = cell(dmax,2); 
    for i = 1:dmax
        sps{i,1} = subplot(dmax+1,2,i*2-1); hold on; if i==1; title(sprintf('P%d',z+1)); end
        plot(ZCav{z}(Li,1:EstDimns(z,3))*Uy(:,i),'Color',clrs_dim(2,:));
        plot(ZCav{z}(Hi,1:EstDimns(z,3))*Uy(:,i),'Color',clrs(2,:));

        sps{i,2} = subplot(dmax+1,2,i*2); hold on; title(sprintf('r = %.2f',r(i))); 
        plot(ZCSHav{z}(Li,1:EstDimns(z,4))*Bx(:,i),'Color',clrs_dim(3,:));
        plot(ZCSHav{z}(Hi,1:EstDimns(z,4))*Bx(:,i),'Color',clrs(3,:));

    end

    for i = 1:size(ZsplitImagery{z},1)
        sess_i = ZsplitImagery{z}{i,2}(:,1:EstDimns(z,3))*Uy; 
        sess_s = ZsplitShared{z}{i,2}(:,1:EstDimns(z,4))*Bx; 
        fig = gcf; 

        for j = 1:size(sps,1)
            fig.CurrentAxes = sps{j,1}; 
            patchline(1:length(Li),sess_i(Li,j),'EdgeColor',clrs_dim(2,:),'EdgeAlpha',0.25); 
            patchline(1:length(Hi),sess_i(Hi,j),'EdgeColor',clrs(2,:),'EdgeAlpha',0.25); 

            fig.CurrentAxes = sps{j,2}; 
            patchline(1:length(Li),sess_s(Li,j),'EdgeColor',clrs_dim(3,:),'EdgeAlpha',0.25); 
            patchline(1:length(Hi),sess_s(Hi,j),'EdgeColor',clrs(3,:),'EdgeAlpha',0.25); 
        end
    end
    clean_plot
    
    dmax = size(Uys,2); 
    figure('Position',[250 72 700 880]); hold on; 
  
    sps = cell(dmax,2); 
    for i = 1:dmax
        sps{i,2} = subplot(dmax+1,2,i*2-1); hold on; title(sprintf('r = %.2f',rs(i))); 
        plot(ZCav{z}(Li,1:EstDimns(z,3))*Bxs(:,i),'Color',clrs_dim(2,:));
        plot(ZCav{z}(Hi,1:EstDimns(z,3))*Bxs(:,i),'Color',clrs(2,:));

        sps{i,1} = subplot(dmax+1,2,i*2); hold on; if i==1; title(sprintf('P%d',z+1)); end
        plot(ZCSHav{z}(Li,1:EstDimns(z,4))*Uys(:,i),'Color',clrs_dim(3,:));
        plot(ZCSHav{z}(Hi,1:EstDimns(z,4))*Uys(:,i),'Color',clrs(3,:));

    end

    for i = 1:size(ZsplitImagery{z},1)
        sess_i = ZsplitImagery{z}{i,2}(:,1:EstDimns(z,3))*Bxs; 
        sess_s = ZsplitShared{z}{i,2}(:,1:EstDimns(z,4))*Uys; 
        fig = gcf; 

        for j = 1:size(sps,1)
            fig.CurrentAxes = sps{j,2}; 
            patchline(1:length(Li),sess_i(Li,j),'EdgeColor',clrs_dim(2,:),'EdgeAlpha',0.25); 
            patchline(1:length(Hi),sess_i(Hi,j),'EdgeColor',clrs(2,:),'EdgeAlpha',0.25); 

            fig.CurrentAxes = sps{j,1}; 
            patchline(1:length(Li),sess_s(Li,j),'EdgeColor',clrs_dim(3,:),'EdgeAlpha',0.25); 
            patchline(1:length(Hi),sess_s(Hi,j),'EdgeColor',clrs(3,:),'EdgeAlpha',0.25); 
        end
    end
    clean_plot
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
