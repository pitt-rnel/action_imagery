%% Monte Carlo correlations (Figure 3a)
numrun = 10000; 
[oc_r,sh_r,oc_fit,opv,nulldcorr,maxdcorr,veOC] = deal(cell(1,2)); 
[ZOav,ZCav,ZOSHav,ZCSHav] = deal(cell(1,2)); 
for z = 1:2
    ZOav{z} = mean(cell2mat(reshape(ZsplitAction{z}(:,1),1,1,[])),3);
           
    ZCav{z} = mean(cell2mat(reshape(ZsplitImagery{z}(:,2),1,1,[])),3);
           
    ZOSHav{z} = mean(cell2mat(reshape(ZsplitShared{z}(:,1),1,1,[])),3);
    
    ZCSHav{z} = mean(cell2mat(reshape(ZsplitShared{z}(:,2),1,1,[])),3);
end
nancorr = @(x,y) corr(x(~any(isnan([x,y]),2),:),y(~any(isnan([x,y]),2),:)); 
for z = 1:2
    for i = 1:numrun
        clc; fprintf('%d/%d (%d/%d)\n',z,2,i,numrun);
        
        [pco,~,~,~,evo] = pca([ZOav{z};ZCav{z}]); 
        od = find(cumsum(evo)>95,1,'first'); 

        [U,~,~] = svd(rand(od,1)*2-1); 
        
        odim = ZOav{z}*pco(:,1:od)*U(:,1); 
        ndim = ZOav{z}*pco(:,1:od)*U(:,2:end); 
        cdim = ZCav{z}*pco(:,1:od)*U(:,1); 
        cnull = ZCav{z}*pco(:,1:od)*U(:,2:end); 
        cfull = ZCav{z}*pco(:,1:od)*U(:,1:end); 
        
        cnull_cent = cnull - nanmean(cnull); 
        cfull_cent = cfull - nanmean(cfull); 
        odim_cent = odim - nanmean(odim); 
        g = ~any(isnan(cnull),2);
        
        b = cnull_cent(g,:)\odim_cent(g); Q = b./norm(b); 
        nulldcorr{z}(i,:) = nancorr(cnull*Q,odim); 
        
        b = cfull_cent(g,:)\odim_cent(g); Q = b./norm(b); 
        maxdcorr{z}(i,:) = nancorr(cfull*Q,odim);
        
        oc_r{z}(i,:) = nancorr(odim,cdim); 
        opv{z}(i,:) = nanvar(odim)./sum(nanvar([ZOav{z},ZOSHav{z}])); 
        
        veOC{z}(i,1) = nanvar(ZOav{z}*pco(:,1:od)*U(:,1))./sum(nanvar([ZOav{z},ZOSHav{z}]))*100; 
        veOC{z}(i,2) = nanvar(ZCav{z}*pco(:,1:od)*U(:,1))./sum(nanvar([ZCav{z},ZCSHav{z}]))*100; 
 
        shd = size(ZOSHav{z},2); 
        [U,~,~] = svd(rand(shd,1)*2-1); 
        
        oshdim = ZOSHav{z}*U(:,1); 
        cshdim = ZCSHav{z}*U(:,1); 
        
        sh_r{z}(i,:) = nancorr(oshdim,cshdim); 
    end
end

%
figure; hold on; 
subplot(2,1,1); hold on; title('shared corr'); 
[xs,ns] = linehist(0:.025:1,sh_r{1},'normalize',100); 
[xs2,ns2] = linehist(0:.025:1,sh_r{2},'normalize',100); 
yl = ylim; 
plot(nanmedian(sh_r{1})*[1 1],yl,'Color',clrs(1,:)); 
plot(nanmedian(sh_r{2})*[1 1],yl,'Color',clrs(2,:)); 
text(nanmedian(sh_r{1}),yl(2)*1.1,sprintf('%.2f',nanmedian(sh_r{1}))); 
text(nanmedian(sh_r{2}),yl(2)*1.1,sprintf('%.2f',nanmedian(sh_r{2}))); 

subplot(2,1,2); hold on; title('action-imagery corr'); 
[x,n] = linehist(0:.025:1,oc_r{1},'normalize',100); 
[x2,n2] = linehist(0:.025:1,oc_r{2},'normalize',100); 
yl = ylim; 
plot(nanmedian(oc_r{1})*[1 1],yl,'Color',clrs(1,:)); 
plot(nanmedian(oc_r{2})*[1 1],yl,'Color',clrs(2,:)); 
text(nanmedian(oc_r{1}),yl(2)*1.1,sprintf('%.2f',nanmedian(oc_r{1}))); 
text(nanmedian(oc_r{2}),yl(2)*1.1,sprintf('%.2f',nanmedian(oc_r{2}))); 

linehist(0:.025:1,nulldcorr{1},'Color',[.7 .7 .7],'normalize',100); 
linehist(0:.025:1,nulldcorr{2},'Color',[.7 .7 .7],'LineType',':','normalize',100); 

plot(nanmedian(nulldcorr{1})*[1 1],yl,'-','Color',[.5 .5 .5]); 
plot(nanmedian(nulldcorr{2})*[1 1],yl,':','Color',[.5 .5 .5]); 

clean_plot; 

r_corr_bounds_sh = cellfun(@(x) prctile(x,[2.5 97.5]),sh_r,'uni',0); 
r_corr_bounds_un = cellfun(@(x) prctile(x,[2.5 97.5]),oc_r,'uni',0); 

[~,p_dimrem(1)] = ttest(oc_r{1},nulldcorr{1}); 
[~,p_dimrem(2)] = ttest(oc_r{2},nulldcorr{2}); 

