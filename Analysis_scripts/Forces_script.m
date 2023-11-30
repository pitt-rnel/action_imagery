%% Forces (Supp Fig 1)
[FORCE] = deal(cell(2,1)); 
for zz = 1:2
    sub = subjects{zz}; 
    for i = 1:length(ActionImageryData.(sub))
        FORCE{zz}{i,:}.Action = ActionImageryData.(sub){i}.Action.Force; 
        FORCE{zz}{i,:}.Imagery = ActionImageryData.(sub){i}.Imagery.Force; 
    end
end
FORCEperc = FORCE;
for zz = 1:2
    sub = subjects{zz}; 
    MVCS = cellfun(@(x) x.Action.MVC_bounds,ActionImageryData.(sub),'uni',0); 
    for z = 1:length(FORCE{zz})
        mvcfunc = @(x,mvc) (x-mvc(1))/(mvc(2)-mvc(1))*100; 
        FORCEperc{zz}{z}.Action{1} = cellfun(@(x) mvcfunc(x,MVCS{z}),FORCEperc{zz}{z}.Action{1},'uni',0); 
        FORCEperc{zz}{z}.Action{2} = cellfun(@(x) mvcfunc(x,MVCS{z}),FORCEperc{zz}{z}.Action{2},'uni',0); 
        FORCEperc{zz}{z}.Imagery{1} = cellfun(@(x) mvcfunc(x,MVCS{z}),FORCEperc{zz}{z}.Imagery{1},'uni',0); 
        FORCEperc{zz}{z}.Imagery{2} = cellfun(@(x) mvcfunc(x,MVCS{z}),FORCEperc{zz}{z}.Imagery{2},'uni',0); 
    end
end

[FORCEC,FORCEC_av] = deal(cell(2,1)); 
for zz  = 1:length(FORCE)
    for z = 1:length(FORCE{zz})
        for i = 1:2
            FORCEC{zz,:}{z,:}.Action{i} = cellfun(@(x) x - nanmean(x(BIND{zz})),FORCEperc{zz}{z}.Action{i},'uni',0); 
            FORCEC{zz,:}{z,:}.Imagery{i} = cellfun(@(x) x - nanmean(x(BIND{zz})),FORCEperc{zz}{z}.Imagery{i},'uni',0); 
            
            FORCEC_av{zz,:}{z,:}.Action{i} = nanmean(cell2mat(reshape(FORCEC{zz}{z}.Action{i},1,1,[])),3); 
            FORCEC_av{zz,:}{z,:}.Imagery{i} = nanmean(cell2mat(reshape(FORCEC{zz}{z}.Imagery{i},1,1,[])),3); 
        end
    end
end

figure('Position',[342 518 1084 475]); hold on; 
yl = [0 60]; 
for zz = 1:2
    subplot(2,2,(zz-1)*2+1); hold on; title('Action (averages)');
    cellfun(@(x) plot(x.Action{1},'Color',nanmean([1 1 1;clrs(1,:)])),FORCEC_av{zz});
    cellfun(@(x) plot(x.Action{2},'Color',clrs(1,:)),FORCEC_av{zz}); 
    ylim(yl); 

    subplot(2,2,(zz-1)*2+2); hold on; title('Imagery (averages)'); 
    cellfun(@(x) plot(x.Imagery{1},'Color',nanmean([1 1 1;clrs(1,:)])),FORCEC_av{zz});
    cellfun(@(x) plot(x.Imagery{2},'Color',clrs(1,:)),FORCEC_av{zz});
    ylim(yl); 
end
clean_plot; 

yl = [-10 80]; 
for zz = 1:2
    figure('Position',[342 518 1084 475]); hold on; 
    for i = 1:length(FORCEC{zz})
        subplot(length(FORCEC{zz}),2,(i-1)*2+1); hold on; title(sprintf('Action P%d (%d)',zz+1,i));
        cellfun(@(x) plot(x,'Color',nanmean([1 1 1;clrs(i,:)])),FORCEC{zz}{i}.Action{1});
        cellfun(@(x) plot(x,'Color',clrs(i,:)),FORCEC{zz}{i}.Action{2}); 
        ylim(yl); set(gca,'YTick',[0 80]); 
    end
    

    for i = 1:length(FORCEC{zz})
        subplot(length(FORCEC{zz}),2,(i-1)*2+2); hold on; title(sprintf('Imagery P%d (%d)',zz+1,i));
        cellfun(@(x) plot(x,'Color',nanmean([1 1 1;clrs(i,:)])),FORCEC{zz}{i}.Imagery{1});
        cellfun(@(x) plot(x,'Color',clrs(i,:)),FORCEC{zz}{i}.Imagery{2}); 
        ylim(yl); set(gca,'YTick',[0 80]); 
    end
    ylim(yl); 
    clean_plot; 
end

%% Check forces during imagery
[FIz,FAz,peakFAz,peakFIz,pfs] = deal(cell(1,2)); 
for z = 1:2
    FI = cellfun(@(x) vertcat(x.Imagery{:}),FORCEperc{z},'uni',0); 
    FIz{z} = vertcat(FI{:}); 
    peakFIz{z} = cellfun(@(x) max(x),FIz{z}); 

    FA = cellfun(@(x) vertcat(x.Action{:}),FORCEperc{z},'uni',0); 
    FAz{z} = vertcat(FA{:}); 
    peakFAz{z} = cellfun(@(x) max(x),FAz{z}); 
end

