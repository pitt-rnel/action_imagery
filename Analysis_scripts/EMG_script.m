%% EMG
[FORCEC,EMGC,EMGC_av,FORCEC_av] = deal(cell(2,1)) ; 

for zz  = 1:length(FORCE)
    sub = subjects{zz}; 
    for z = 1:length(FORCE{zz})
        for i = 1:2
            FORCEC{zz,:}{z,:}.Action{i} = cellfun(@(x) x - nanmean(x(BIND{zz})),FORCE{zz}{z}.Action{i},'uni',0); 
            FORCEC{zz,:}{z,:}.Imagery{i} = cellfun(@(x) x - nanmean(x(BIND{zz})),FORCE{zz}{z}.Imagery{i},'uni',0); 
            
            FORCEC_av{zz,:}{z,:}.Action{i} = nanmean(cell2mat(reshape(FORCEC{zz}{z}.Action{i},1,1,[])),3); 
            FORCEC_av{zz,:}{z,:}.Imagery{i} = nanmean(cell2mat(reshape(FORCEC{zz}{z}.Imagery{i},1,1,[])),3); 
            
            if isfield(ActionImageryData.(sub){z}.Action,'EMG')
                EMGC{zz,:}{z,:}.Action{i} = cellfun(@(x) sqrt(sum(x.^2)),ActionImageryData.(sub){z}.Action.EMG{i},'uni',0); 
                EMGC{zz,:}{z,:}.Imagery{i} = cellfun(@(x) sqrt(sum(x.^2)),ActionImageryData.(sub){z}.Imagery.EMG{i},'uni',0); 

                EMGC_av{zz,:}{z,:}.Action{i} = nanmean(cell2mat(reshape(EMGC{zz}{z}.Action{i},1,1,[])),3); 
                EMGC_av{zz,:}{z,:}.Imagery{i} = nanmean(cell2mat(reshape(EMGC{zz}{z}.Imagery{i},1,1,[])),3); 
            end
            
        end
    end
end

figure; hold on; 
subplot(2,1,1); hold on; title('EMG'); 
E1 = EMGC{1}{1}.Action;
F1 = FORCEC{1}{1}.Action;
cellfun(@(x) plot(x,'Color',nanmean([1 1 1; clrs(1,:)])),E1{1}); 
cellfun(@(x) plot(x,'Color',clrs(1,:)),E1{2}); 
xlim([0 400]); 
subplot(2,1,2); hold on; title('Force'); 
cellfun(@(x) plot(x,'Color',[.5 .5 .5]),F1{1}); 
cellfun(@(x) plot(x,'Color',[0 0 0]),F1{2}); 
xlim([0 400]); 
clean_plot; 

figure; hold on; 
clo = lines(10); 
peaklag = []; 
for i = 1:length(EMGC{1})
    if ~isempty(EMGC{1}{i})
        E1 = EMGC{1}{i}.Action;
        F1 = FORCEC{1}{i}.Action;
        E1all = cell2mat(vertcat(E1{:})')'; 
        F1all = cell2mat(vertcat(F1{:})')'; 
        g = find(~any(isnan([E1all,F1all]),2)); 
        [r,lags] = xcorr(E1all(g),F1all(g),50,'coeff'); 
        
        plot(lags*.02*1000,r,'Color',clo(i,:)); 

        peaklag = [peaklag, lags(find(r==max(r),1,'first'))*.02*1000]; 
    end
end
yl = ylim; 
plot([0 0],yl,'k:'); ylim(yl); 
xlabel('lag (ms)'); 
ylabel('correlation'); 
plot(mean(peaklag)*[1 1],[0 yl(2)],'Color',clo(1,:)); 
t = text(mean(peaklag),yl(1),sprintf('offset = %.0f',mean(peaklag)),'HorizontalAlignment','left','VerticalAlignment','bottom'); 
t.Rotation = 90; 
clean_plot; 


yl = [0 3000]; 
zz = 1; % only P2 has EMG recordings
figure('Position',[342 518 1084 475]); hold on; 
for i = 1:length(EMGC{zz})
    if ~isempty(EMGC{zz}{i})
        subplot(length(EMGC{zz}),2,(i-1)*2+1); hold on; title(sprintf('EMG P%d (%d)',zz+1,i));
        for j = 1:length(EMGC{zz}{i}.Action{1})
        	plot(EMGC{zz}{i}.Action{1}{j},'Color',nanmean([1 1 1;clrs(i,:)]));
        end
        for j = 1:length(EMGC{zz}{i}.Action{2}) 
            plot(EMGC{zz}{i}.Action{2}{j},'Color',clrs(i,:));
        end
        for j = 1:length(EMGC{zz}{i}.Imagery{1})
        	plot((1:length(EMGC{zz}{i}.Imagery{1}{j})) + length(EMGC{zz}{i}.Imagery{1}{j}) + 25, ...
                EMGC{zz}{i}.Imagery{1}{j},'Color',nanmean([1 1 1;clrs(i,:)]));
        end
        for j = 1:length(EMGC{zz}{i}.Imagery{2}) 
            plot((1:length(EMGC{zz}{i}.Imagery{2}{j})) + length(EMGC{zz}{i}.Imagery{2}{j}) + 25, ...
                EMGC{zz}{i}.Imagery{2}{j},'Color',clrs(i,:));
        end
        ylim(yl); set(gca,'YTick',yl); 
    end
    ylabel('EMG magnitude'); 
end

yl = [-10 80]; 
for i = 1:length(FORCEC{zz})
    if ~isempty(EMGC{zz}{i})
        subplot(length(FORCEC{zz}),2,(i-1)*2+2); hold on; title(sprintf('Force P%d (%d)',zz+1,i));
        for j = 1:length(FORCEC{zz}{i}.Action{1})
        	plot(FORCEC{zz}{i}.Action{1}{j},'Color',nanmean([1 1 1;clrs(i,:)]));
        end
        for j = 1:length(FORCEC{zz}{i}.Action{2}) 
            plot(FORCEC{zz}{i}.Action{2}{j},'Color',clrs(i,:));
        end
        for j = 1:length(FORCEC{zz}{i}.Imagery{1})
        	plot((1:length(FORCEC{zz}{i}.Imagery{1}{j})) + length(FORCEC{zz}{i}.Imagery{1}{j}) + 25, ...
                FORCEC{zz}{i}.Imagery{1}{j},'Color',nanmean([1 1 1;clrs(i,:)]));
        end
        for j = 1:length(FORCEC{zz}{i}.Imagery{2}) 
            plot((1:length(FORCEC{zz}{i}.Imagery{2}{j})) + length(FORCEC{zz}{i}.Imagery{2}{j}) + 25, ...
                FORCEC{zz}{i}.Imagery{2}{j},'Color',clrs(i,:));
        end
        ylim(yl); set(gca,'YTick',[0 80]); 
    end
    ylabel('Force (N)'); 
end
ylim(yl); 
clean_plot; 
