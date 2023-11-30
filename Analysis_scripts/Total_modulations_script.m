%% total modulations script (Figure 1)
figure; hold on; 
mrk = {'-',':'}; 
for i = 1:2
    sub = subjects{i}; 
    MATO = cellfun(@(x) x.Action.Motor_raw,ActionImageryData.(sub),'uni',0); MATO = vertcat(MATO{:,:}); 
    MATC = cellfun(@(x) x.Imagery.Motor_raw,ActionImageryData.(sub),'uni',0); MATC = vertcat(MATC{:,:}); 

    actavs = cellfun(@(x) mean(cell2mat(reshape(x,1,1,[])),3)'*50,MATO,'uni',0); 
    modavs = cellfun(@(x) mean(x - repmat(nanmean(x(find(~isnan(x(:,1)),10,'first'),:)),size(x,1),1),2),actavs,'uni',0); 
    s1 = subplot(2,1,1); hold on; title('Motor'); 

    cellfun(@(x) plot(x,mrk{i},'Color',clrs(1,:)),modavs(:,1)); 
    cellfun(@(x) plot(x,mrk{i},'Color',clrs(2,:)),modavs(:,2)); 
    
end
ylim_motO = ylim;

for i = 1:2
    sub = subjects{i}; 
    MATO = cellfun(@(x) x.Action.Motor_raw,ActionImageryData.(sub),'uni',0); MATO = vertcat(MATO{:,:}); 
    MATC = cellfun(@(x) x.Imagery.Motor_raw,ActionImageryData.(sub),'uni',0); MATC = vertcat(MATC{:,:}); 
    
    actavs = cellfun(@(x) mean(cell2mat(reshape(x,1,1,[])),3)'*50,MATC,'uni',0); 
    modavs = cellfun(@(x) mean(x - repmat(nanmean(x(find(~isnan(x(:,1)),10,'first'),:)),size(x,1),1),2),actavs,'uni',0); 
    s2 = subplot(2,1,2); hold on; 
    
    cellfun(@(x) plot(x,mrk{i},'Color',clrs(1,:)),modavs(:,1)); 
    cellfun(@(x) plot(x,mrk{i},'Color',clrs(2,:)),modavs(:,2)); 
    
end
ylim_motC = ylim; 
ylim_mot = [min([ylim_motO,ylim_motC]),max([ylim_motO,ylim_motC])]; 
ylim(s1,ylim_mot); 
ylim(s2,ylim_mot); 
clean_plot; 
clean_plot(s1,'timebar');
clean_plot(s2,'timebar'); 

