%% Single Channel modulations
mrkrs = {nanmean([1 1 1 ; clrs(1,:)]),clrs(1,:); ...
         nanmean([1 1 1 ; clrs(2,:)]),clrs(2,:)}; 
fh = figure('Position',[544 443 1067 528.4000]); hold on; 
ratioLH = deal(cell(2,1)); 
zls = {[0 150],[0 100]};
clrs = lines(10); clrs2 = nanmean(cat(3,ones(10,3),clrs),3); 
sp = cell(1,2); 
for z = 1:2 
    sp{z} = subplot(1,2,z); hold on; 
    sub = subjects{z}; 
    MATO = cellfun(@(x) x.Action.Motor_raw,ActionImageryData.(sub),'uni',0); MATO = vertcat(MATO{:,:}); 
    MATC = cellfun(@(x) x.Imagery.Motor_raw,ActionImageryData.(sub),'uni',0); MATC = vertcat(MATC{:,:}); 

    matosL = cell2mat(cellfun(@(x) mean(cell2mat(reshape(x,1,1,[])),3)',MATO(:,1)','uni',0)); 
    matcsL = cell2mat(cellfun(@(x) mean(cell2mat(reshape(x,1,1,[])),3)',MATC(:,1)','uni',0)); 
    matosH = cell2mat(cellfun(@(x) mean(cell2mat(reshape(x,1,1,[])),3)',MATO(:,2)','uni',0)); 
    matcsH = cell2mat(cellfun(@(x) mean(cell2mat(reshape(x,1,1,[])),3)',MATC(:,2)','uni',0)); 
    
    OC_L = [diff(prctile(matosL*50,[5 95]))', diff(prctile(matcsL*50,[5 95]))']; 
    OC_H = [diff(prctile(matosH*50,[5 95]))', diff(prctile(matcsH*50,[5 95]))']; 
    
    OCALLmod = diff(prctile([matosL;matcsL;matosH;matcsH]*50,[5 95])); % times 50 to go from bin count -> FR (Hz)
    OCALLmu = nanmean([matosL;matcsL;matosH;matcsH])*50; 
    
    badL = OC_L(:,1)==0; 
    badH = OC_H(:,1)==0; 
    
    OC_L(badL,:) = []; 
    OC_H(badH,:) = []; 
    matosL(:,badL) = []; matcsL(:,badL) = []; 
    matosH(:,badH) = []; matcsH(:,badH) = []; 
    
    [~,~,ratioLH{z}(:,1)] = boot_bounds(10000,@(X) X(:,1)\X(:,2),OC_L,2.5,97.5); 
    [~,~,ratioLH{z}(:,2)] = boot_bounds(10000,@(X) X(:,1)\X(:,2),OC_H,2.5,97.5); 
    
    slope_bnd = prctile(ratioLH{z},[2.5 97.5]); 
    slope_av = [OC_L(:,1)\OC_L(:,2), OC_H(:,1)\OC_H(:,2)];  
    
    plot(zls{z},slope_av(1)*zls{z},'Color',mrkrs{z,1},'LineWidth',2); 
    patch([zls{z} fliplr(zls{z})],[0, slope_bnd(1,1)*zls{z}(2), slope_bnd(2,1)*zls{z}(2), 0],mrkrs{z,1},'EdgeColor','none','FaceAlpha',0.2); 
    
    plot(zls{z},slope_av(2)*zls{z},'Color',mrkrs{z,2},'LineWidth',2); 
    patch([zls{z} fliplr(zls{z})],[0, slope_bnd(1,2)*zls{z}(2), slope_bnd(2,2)*zls{z}(2), 0],mrkrs{z,2},'EdgeColor','none','FaceAlpha',0.2); 
    
    title(sprintf('%.0f +/- %.0f    %.0f +/- %.0f',100*slope_av(1),100*slope_bnd(2,1)-100*slope_av(1),100*slope_av(2),100*slope_bnd(2,2)-100*slope_av(2))); 
    
    plot(OC_L(:,1),OC_L(:,2),'.','MarkerSize',15,'Color',mrkrs{z,1}); 
    plot(OC_H(:,1),OC_H(:,2),'.','MarkerSize',15,'Color',mrkrs{z,2}); 
    
    plot(zls{z},zls{z},'k:'); xlim(zls{z}); ylim(zls{z}); axis square; 
    xlabel('action modulation (Hz)'); ylabel('imagery modulation (Hz)'); 
   
end
clean_plot; 