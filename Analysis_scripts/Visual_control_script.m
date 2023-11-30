%% Visual Control
[ZA,FA,FFA,ZFA,Viscue] = deal(cell(2,1)); 
emptzer = @(x) max([x,0]); 
for zz = 1:2
    sub = subjects{zz}; 
    Viscue{zz} = cellfun(@(x) x.Action.HighForceCueTimes,ActionImageryData.(sub),'uni',0); 
    Zsub = cellfun(@(x) x.Action.PC_combined{2},ActionImageryData.(sub),'uni',0); 
    ZFsub = cellfun(@(x) x.Action.PC_combined_extended{2},ActionImageryData.(sub),'uni',0); 
    FA{zz} = cellfun(@(x) x.Action.Force{2},ActionImageryData.(sub),'uni',0); 
    FFA{zz} = cellfun(@(x) x.Action.Force_extended{2},ActionImageryData.(sub),'uni',0); 
    
    npc = size(Subspaces{zz}{1,1},1); 
    for i = 1:length(Zsub)

        ZA{zz}{i,1} = cellfun(@(x) x(1:npc,:)'*Subspaces{zz}{i,1},Zsub{i},'uni',0); 
        ZA{zz}{i,2} = cellfun(@(x) x(1:npc,:)'*Subspaces{zz}{i,2},Zsub{i},'uni',0); 
        ZA{zz}{i,3} = cellfun(@(x) x(1:npc,:)'*Subspaces{zz}{i,3},Zsub{i},'uni',0); 
        
        ZFA{zz}{i,1} = cellfun(@(x) x(1:npc,:)'*Subspaces{zz}{i,1},ZFsub{i},'uni',0); 
        ZFA{zz}{i,2} = cellfun(@(x) x(1:npc,:)'*Subspaces{zz}{i,2},ZFsub{i},'uni',0); 
        ZFA{zz}{i,3} = cellfun(@(x) x(1:npc,:)'*Subspaces{zz}{i,3},ZFsub{i},'uni',0); 
        
        for j = 1:length(ZA{zz}{i,1})
            removenan = any(isnan(ZFA{zz}{i,1}{j}),2); 
            ZFA{zz}{i,1}{j}(removenan,:) = []; 
            ZFA{zz}{i,2}{j}(removenan,:) = []; 
            ZFA{zz}{i,3}{j}(removenan,:) = []; 
            FFA{zz}{i}{j}(removenan) = []; 
            lastnan = find(~removenan,1,'first')-1; 
            Viscue{zz}{i}(j) = Viscue{zz}{i}(j)-lastnan; 
       
            removenan2 = any(isnan(ZA{zz}{i,1}{j}),2); 
            lastnan2 = find(~removenan2,1,'first')-1; 
        end

    end
end

% Unique space
[forceon,tranon,comps,forces] = deal(cell(2,1)); 
ss = 1; 
transfinder = @(x) find(nansum(x(round(end/2):(round(end/2)+150),:).^2)-nansum(x((round(end/2)+151):end,:).^2) ...
    == max(nansum(x(round(end/2):(round(end/2)+150),:).^2)-nansum(x((round(end/2)+151):end,:).^2))); 

tcoms = cellfun(transfinder,ZOav); % index of onset transient components for P2 and P3
for z = 1:2
    trancom = tcoms(z); 
    for j = 1:length(FFA{z}) % sessions 
        % Transient threshold
        peak = cellfun(@(x) max(x(:,trancom)),ZFA{z}{j,ss}); 
        val = cellfun(@(x) min(x(:,trancom)),ZFA{z}{j,ss}); 
        possible_threshs = linspace(mean(peak),mean(val),100); 
        num_crosses = NaN(length(possible_threshs),1); 
        Z = ZFA{z}{j,ss}'; Z = cellfun(@(x) x(~any(isnan(x),2),:),Z,'uni',0); 
        for k = 1:length(possible_threshs)
            find(~isnan(x(:,trancom)),20,'first'); 
            num_crosses(k,:) = mean(cellfun(@(x) sum(diff(x(1:find(x(:,trancom)==max(x(:,trancom)),1,'first'),trancom)>possible_threshs(k))==1) ,Z)==1); 
        end
        thresh = possible_threshs(find(num_crosses==max(num_crosses),1,'last')); 
        missout = find(cellfun(@(x) sum(diff(x(1:find(x(:,trancom)==max(x(:,trancom)),1,'first'),trancom)>thresh)==1) ,Z)~=1); 

        % force threshold
        allfon = nanmean(cell2mat(cellfun(@(x) x,FA{z}{j},'uni',0)),1)';
        peak = find(allfon==max(allfon),1,'first'); 
        possible_threshs = linspace(allfon(peak),allfon(1),100); 
        num_crosses = NaN(length(possible_threshs),1); 
        ZF = FFA{z}{j}'; ZF = cellfun(@(x) x(~isnan(x)),ZF,'uni',0); 
        for k = 1:length(possible_threshs)
            num_crosses(k,:) = mean(cellfun(@(x) sum(diff(x(1:round(end/2))>possible_threshs(k))==1),ZF)==1); 
        end
        fthresh = possible_threshs(find(num_crosses==max(num_crosses),1,'last')); 
        
        for i = 1:length(FFA{z}{j})
            forceon{z}{j,:}(i,:) = .02*(emptzer(find(FFA{z}{j}{i}>fthresh,1,'first'))-Viscue{z}{j}(i)); 
            tranon{z}{j,:}(i,:) = .02*(emptzer(find(ZFA{z}{j,ss}{i}(:,trancom)>thresh,1,'first'))-Viscue{z}{j}(i)); 
        end
        misses = [missout'; find(tranon{z}{j,:}==0)]; 
        forceon{z}{j,:}(misses) = NaN; 
        tranon{z}{j,:}(misses) = NaN; 
        
        for i = 1:length(tranon{z}{j,:})
            if ~isnan(tranon{z}{j,:}(i))
                foni =forceon{z}{j,:}(i)/.02 + Viscue{z}{j}(i); 
                if foni > 50
                    fini = round((foni-50):(foni+75)); 
                    comps{z}{j,:}{i,:} = ZFA{z}{j,ss}{i}(fini,trancom)-thresh; 
                    forces{z}{j,:}{i,:} = FFA{z}{j}{i}(fini)-fthresh; 
                else
                    fini = round(1:(foni+75)); 
                    ppad = round(50-foni); 
                    comps{z}{j,:}{i,:} = [NaN(ppad,1); ZFA{z}{j,ss}{i}(fini,trancom)]-thresh; 
                    forces{z}{j,:}{i,:} = [NaN(ppad,1)', FFA{z}{j}{i}(fini)]'-fthresh; 
                end
            end
        end
        emps = cellfun(@isempty,comps{z}{j}); 
        comps{z}{j}(emps) = []; 
        forces{z}{j}(emps) = []; 
    end
end
figure('Position',[23 503 1819 450]); 

subplot(2,3,1); hold on; title('Unique'); 
cellfun(@(x) patchline((1:length(x))*.02,x,'EdgeColor',clrs(1,:),'EdgeAlpha',0.2),vertcat(comps{1}{:})); 
yyaxis right; 
cellfun(@(x) patchline((1:length(x))*.02,x,'EdgeColor',[0 0 0],'EdgeAlpha',0.2),vertcat(forces{1}{:})); 
ylabel('Newtons'); 

subplot(2,3,4); hold on; 
cellfun(@(x) patchline((1:length(x))*.02,x,'EdgeColor',clrs(1,:),'EdgeAlpha',0.2),vertcat(comps{2}{:})); 
yyaxis right; 
cellfun(@(x) patchline((1:length(x))*.02,x,'EdgeColor',[0 0 0],'EdgeAlpha',0.2),vertcat(forces{2}{:})); 
ylabel('Newtons'); 

subplot(2,3,[2 5]); hold on; 
plot(cell2mat(tranon{1}),cell2mat(forceon{1}),'k.');
plot(cell2mat(tranon{2}),cell2mat(forceon{2}),'ko');
axis square 
xlabel('transient onset (s)'); 
ylabel('force onset (s)'); 
plot([-.5 1.5],[-.5 1.5],'k:'); 
xlim([-.5 1.5]),ylim([-.5 1.5]); 

subplot(2,3,[3 6]); hold on; 
[x1,h1] = linehist(-1:0.05:1,cell2mat(forceon{1})-cell2mat(tranon{1}),'suppress_plot'); 
[x2,h2] = linehist(-1:0.05:1,cell2mat(forceon{2})-cell2mat(tranon{2}),'suppress_plot'); 

patch(x1,h1,clrs(1,:),'FaceAlpha',0.2,'EdgeColor','none'); 
patch(x2,h2,clrs(2,:),'FaceAlpha',0.2,'EdgeColor','none'); 
xlabel('transient-force lag'); 

clean_plot;

% Shared space
[forceon,tranon,comps,forces] = deal(cell(2,1)); 
ss = 3; 
tcoms = [1 2]; % index of onset transient components for P2 and P3
for z = 1:2
    trancom = tcoms(z); 
    for j = 1:length(FFA{z}) % sessions 
        % Transient threshold
        peak = cellfun(@(x) max(x(:,trancom)),ZFA{z}{j,ss}); 
        val = cellfun(@(x) min(x(:,trancom)),ZFA{z}{j,ss}); 
        possible_threshs = linspace(mean(peak),mean(val),100); 
        num_crosses = NaN(length(possible_threshs),1); 
        Z = ZFA{z}{j,ss}'; Z = cellfun(@(x) x(~any(isnan(x),2),:),Z,'uni',0); 
        for k = 1:length(possible_threshs)
            find(~isnan(x(:,trancom)),20,'first'); 
            num_crosses(k,:) = mean(cellfun(@(x) sum(diff(x(1:find(x(:,trancom)==max(x(:,trancom)),1,'first'),trancom)>possible_threshs(k))==1) ,Z)==1); 
        end
        thresh = possible_threshs(find(num_crosses==max(num_crosses),1,'last')); 
        missout = find(cellfun(@(x) sum(diff(x(1:find(x(:,trancom)==max(x(:,trancom)),1,'first'),trancom)>thresh)==1) ,Z)~=1); 

        % force threshold
        allfon = nanmean(cell2mat(cellfun(@(x) x,FA{z}{j},'uni',0)),1)';
        peak = find(allfon==max(allfon),1,'first'); 
        possible_threshs = linspace(allfon(peak),allfon(1),100); 
        num_crosses = NaN(length(possible_threshs),1); 
        ZF = FFA{z}{j}'; ZF = cellfun(@(x) x(~isnan(x)),ZF,'uni',0); 
        for k = 1:length(possible_threshs)
            num_crosses(k,:) = mean(cellfun(@(x) sum(diff(x(1:round(end/2))>possible_threshs(k))==1),ZF)==1); 
        end
        fthresh = possible_threshs(find(num_crosses==max(num_crosses),1,'last')); 
        
        for i = 1:length(FFA{z}{j})
            forceon{z}{j,:}(i,:) = .02*(emptzer(find(FFA{z}{j}{i}>fthresh,1,'first'))-Viscue{z}{j}(i)); 
            tranon{z}{j,:}(i,:) = .02*(emptzer(find(ZFA{z}{j,ss}{i}(:,trancom)>thresh,1,'first'))-Viscue{z}{j}(i)); 
        end
        misses = [missout'; find(tranon{z}{j,:}==0)]; 
        forceon{z}{j,:}(misses) = NaN; 
        tranon{z}{j,:}(misses) = NaN; 
        
        for i = 1:length(tranon{z}{j,:})
            if ~isnan(tranon{z}{j,:}(i))
                foni =forceon{z}{j,:}(i)/.02 + Viscue{z}{j}(i); 
                if foni > 50
                    fini = round((foni-50):(foni+75)); 
                    comps{z}{j,:}{i,:} = ZFA{z}{j,ss}{i}(fini,trancom)-thresh; 
                    forces{z}{j,:}{i,:} = FFA{z}{j}{i}(fini)-fthresh; 
                else
                    fini = round(1:(foni+75)); 
                    ppad = round(50-foni); 
                    comps{z}{j,:}{i,:} = [NaN(ppad,1); ZFA{z}{j,ss}{i}(fini,trancom)]-thresh; 
                    forces{z}{j,:}{i,:} = [NaN(ppad,1)', FFA{z}{j}{i}(fini)]'-fthresh; 
                end
            end
        end
        emps = cellfun(@isempty,comps{z}{j}); 
        comps{z}{j}(emps) = []; 
        forces{z}{j}(emps) = []; 
    end
end
figure('Position',[23 503 1819 450]); 

subplot(2,3,1); hold on; title('Shared'); 
cellfun(@(x) patchline((1:length(x))*.02,x,'EdgeColor',clrs(3,:),'EdgeAlpha',0.2),vertcat(comps{1}{:})); 
yyaxis right; 
cellfun(@(x) patchline((1:length(x))*.02,x,'EdgeColor',[0 0 0],'EdgeAlpha',0.2),vertcat(forces{1}{:})); 
ylabel('Newtons'); 

subplot(2,3,4); hold on; 
cellfun(@(x) patchline((1:length(x))*.02,x,'EdgeColor',clrs(3,:),'EdgeAlpha',0.2),vertcat(comps{2}{:})); 
yyaxis right; 
cellfun(@(x) patchline((1:length(x))*.02,x,'EdgeColor',[0 0 0],'EdgeAlpha',0.2),vertcat(forces{2}{:})); 
ylabel('Newtons'); 

subplot(2,3,[2 5]); hold on; 
plot(cell2mat(tranon{1}),cell2mat(forceon{1}),'k.');
plot(cell2mat(tranon{2}),cell2mat(forceon{2}),'ko');
axis square 
xlabel('transient onset (s)'); 
ylabel('force onset (s)'); 
plot([-.5 1.5],[-.5 1.5],'k:'); 
xlim([-.5 1.5]),ylim([-.5 1.5]); 

subplot(2,3,[3 6]); hold on; 
[x1,h1] = linehist(-1:0.05:1,cell2mat(forceon{1})-cell2mat(tranon{1}),'suppress_plot'); 
[x2,h2] = linehist(-1:0.05:1,cell2mat(forceon{2})-cell2mat(tranon{2}),'suppress_plot'); 

patch(x1,h1,clrs(1,:),'FaceAlpha',0.2,'EdgeColor','none'); 
patch(x2,h2,clrs(2,:),'FaceAlpha',0.2,'EdgeColor','none'); 
xlabel('transient-force lag'); 

clean_plot;
