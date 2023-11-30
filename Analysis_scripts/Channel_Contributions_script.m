%% Channel Contributions (Supp Fig 5)
[contrib_A,contrib_I,contrib_S,contribs] = deal(cell(2,1));
catsub = vertcat(Subspaces{:}); 
for zz = 1:2
    sub = subjects{zz}; 
    for z = 1:length(ActionImageryData.(sub))
        act = [ActionImageryData.(sub){z}.Action.Motor_full,ActionImageryData.(sub){z}.Imagery.Motor_full];
        condav = cell2mat(cellfun(@(x) mean(cell2mat(reshape(x,1,1,[])),3),act,'uni',0))'; 
        mus = nanmean(condav,1); 

        lowDproj = ActionImageryData.(sub){z}.processing.Combined_ActionImagery_PC_weights(:,1:size(catsub{z,1},1)); 
        SA = lowDproj*catsub{z,1}; 
        SI = lowDproj*catsub{z,2}; 
        SS = lowDproj*catsub{z,3}; 

        for i = 1:size(condav,2) %channels
            contrib_A{z,1}(i,1) = norm(SA(i,:))*mus(i); 
            contrib_I{z,1}(i,1) = norm(SI(i,:))*mus(i); 
            contrib_S{z,1}(i,1) = norm(SS(i,:))*mus(i); 
        end
    end
    contribs{zz}.A = contrib_A; 
    contribs{zz}.I = contrib_I; 
    contribs{zz}.S = contrib_S; 
end

[pdip,pdip_log] = deal(cell(1,2)); 
for z = 1:2
    figure; hold on; 
    
    subplot(3,2,1); hold on; title('Action-unique'); 
    linehist(0:0.01:.5,cell2mat(contribs{z}.A),'normalize',1); 
    xlabel('contribution'); ylabel('proportion'); 
    pdip{z}(1) = dipTest(cell2mat(contribs{z}.A)); 
    
    subplot(3,2,3); hold on; title('Imagery-unique'); 
    linehist(0:0.01:.5,cell2mat(contribs{z}.I),'normalize',1); 
    xlabel('contribution'); ylabel('proportion'); 
    pdip{z}(2) = dipTest(cell2mat(contribs{z}.I)); 
    
    subplot(3,2,5); hold on; title('Shared'); 
    linehist(0:0.01:.5,cell2mat(contribs{z}.S),'normalize',1); 
    xlabel('contribution'); ylabel('proportion');
    pdip{z}(3) = dipTest(cell2mat(contribs{z}.S)); 
    
    
    subplot(3,2,2); hold on; title('Action-unique & Imagery-unique'); 
    linehist(-3:0.1:3,log(cell2mat(contribs{z}.A)./cell2mat(contribs{z}.I)),'normalize',1); 
    xlabel('log ratio of contribution'); ylabel('proportion'); 
    pdip_log{z}(1) = dipTest(log(cell2mat(contribs{z}.A)./cell2mat(contribs{z}.I))); 
    
    subplot(3,2,4); hold on; title('Action-unique & Shared'); 
    linehist(-3:0.1:3,log(cell2mat(contribs{z}.A)./cell2mat(contribs{z}.S)),'normalize',1); 
    xlabel('log ratio of contribution'); ylabel('proportion'); 
    pdip_log{z}(2) = dipTest(log(cell2mat(contribs{z}.A)./cell2mat(contribs{z}.S))); 
    
    subplot(3,2,6); hold on; title('Imagery-unique & Shared'); 
    linehist(-3:0.1:3,log(cell2mat(contribs{z}.I)./cell2mat(contribs{z}.S)),'normalize',1); 
    xlabel('log ratio of contribution'); ylabel('proportion'); 
    pdip_log{z}(3) = dipTest(log(cell2mat(contribs{z}.I)./cell2mat(contribs{z}.S))); 
    
    clean_plot; 
end
