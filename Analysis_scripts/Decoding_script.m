%% Weiner Cascade Decoding performance
R2fun = @(Y,Ypred) 1 - sum((Ypred(~any(isnan([Ypred,Y]),2),:)-Y(~any(isnan([Ypred,Y]),2),:)).^2)./...
    sum((Y(~any(isnan([Ypred,Y]),2),:)-nanmean(Y(~any(isnan([Ypred,Y]),2),:))).^2); 
[ZA,FA,TFA,FAU,FAS] = deal(cell(2,1)); 
nboot = 1000; 
clean2target = false; 
include_imagery = false; 
temp_offset = false; 

for zz = 1:2
   sub = subjects{zz}; 
   Zsub = cellfun(@(x) vertcat(x.Action.PC_combined{:}),ActionImageryData.(sub),'uni',0); 
   FA{zz} = cellfun(@(x) vertcat(x.Action{:}),FORCEperc{zz},'uni',0); 
   FAU{zz} = FA{zz}; 
   FAS{zz} = FA{zz}; 
    if temp_offset
        ForceTrace = nanmean(cell2mat(reshape(cellfun(@(x) cell2mat(x.Action),FORCEC_av{zz},'uni',0),1,1,[])),3); 
        XU = ZOav{zz};
        XS = ZOSHav{zz}; %[ZOav{zz},ZOSHav{zz}]; 
        offs = 0:50; 
        R2_tempoff = zeros(length(offs),3); 
        for offi = 1:length(offs)
            F = [ForceTrace((offs(offi)+1):end),repmat(ForceTrace(end),1,offs(offi))]'; 
            g = ~any(isnan([XU,XS,F]),2); 
            bo = [ones(sum(g),1) XU(g,:)]\F(g); 
            rec = [ones(size(XU,1),1) XU]*bo; 
            R2_tempoff(offi,1) = R2fun(F,rec); 
            
            bo = [ones(sum(g),1) XS(g,:)]\F(g); 
            rec = [ones(size(XS,1),1) XS]*bo; 
            R2_tempoff(offi,2) = R2fun(F,rec); 
            
            bo = [ones(sum(g),1) XU(g,:) XS(g,:)]\F(g); 
            rec = [ones(size(XS,1),1) XU XS]*bo; 
            R2_tempoff(offi,3) = R2fun(F,rec); 
        end
        f_lag = round(mean([find(R2_tempoff(:,1)==max(R2_tempoff(:,1)),1,'first'),... 
                 find(R2_tempoff(:,2)==max(R2_tempoff(:,2)),1,'first')])); 
        f_lagU = find(R2_tempoff(:,1)==max(R2_tempoff(:,1)),1,'first');      
        f_lagS = find(R2_tempoff(:,2)==max(R2_tempoff(:,2)),1,'first'); 
        for i = 1:length(FA{zz}) % sess
            for j = 1:length(FA{zz}{i}) %trials
                FA{zz}{i}{j} = [FA{zz}{i}{j}(f_lag:end),repmat(FA{zz}{i}{j}(end),1,(f_lag-1))]; 
                nn = find(~isnan(FA{zz}{i}{j})); 
                FA{zz}{i}{j} = interp1(nn,FA{zz}{i}{j}(nn),1:length(FA{zz}{i}{j})); 
                
                FAU{zz}{i}{j} = [FAU{zz}{i}{j}(f_lagU:end),repmat(FAU{zz}{i}{j}(end),1,(f_lagU-1))]; 
                nn = find(~isnan(FAU{zz}{i}{j})); 
                FAU{zz}{i}{j} = interp1(nn,FAU{zz}{i}{j}(nn),1:length(FAU{zz}{i}{j})); 
                
                FAS{zz}{i}{j} = [FAS{zz}{i}{j}(f_lagS:end),repmat(FAS{zz}{i}{j}(end),1,(f_lagS-1))]; 
                nn = find(~isnan(FAS{zz}{i}{j})); 
                FAS{zz}{i}{j} = interp1(nn,FAS{zz}{i}{j}(nn),1:length(FAS{zz}{i}{j})); 
            end
        end
             
    end
             
    npc = size(Subspaces{zz}{1,1},1); 
    for i = 1:length(Zsub)
        ZA{zz}{i,1} = cellfun(@(x) x(1:npc,:)'*Subspaces{zz}{i,1},Zsub{i},'uni',0); 
        ZA{zz}{i,2} = cellfun(@(x) x(1:npc,:)'*Subspaces{zz}{i,2},Zsub{i},'uni',0); 
        ZA{zz}{i,3} = cellfun(@(x) x(1:npc,:)'*Subspaces{zz}{i,3},Zsub{i},'uni',0); 
        
        for j = 1:size(UY{zz},2)
            na = size(UY{zz},2); 
            ZA{zz}{i,j+3} = cellfun(@(x) x(1:npc,:)'*Subspaces{zz}{i,1}(:,1:na)*UY{zz}(:,j:end),Zsub{i},'uni',0); 
%             noj = 1:size(UY{zz},2); noj(j) = []; 
%             ZA{zz}{i,j+3} = cellfun(@(x) x(1:npc,:)'*Subspaces{zz}{i,1}(:,1:na)*UY{zz}(:,noj),Zsub{i},'uni',0); 
        end
    end
end

[ULpred,SLpred,CLpred,HCLpred,YLact,R2regU,R2regS,R2regC,R2regHC,Yreg,HULpred,HSLpred,YregH,...
    R2regHU,R2regHS,R2regMU,R2regMS,R2regMC,Ytrial,Upredtrial,Spredtrial,Upreds,Spreds,SSpreds,Yacts,YUacts,YSacts] = deal(cell(2,1)); 
for zz = 1:2
    sub = subjects{zz}; 
    for q = 1:size(ZA{zz},1)
        if include_imagery
            tiz = [cell2mat(ActionImageryData.(sub){q}.Action.trial_num), cell2mat(ActionImageryData.(sub){q}.Imagery.trial_num)+max(cell2mat(ActionImageryData.(sub){q}.Action.trial_num))]; 
        else
            tiz = cell2mat(ActionImageryData.(sub){q}.Action.trial_num); % action trial orders
        end
        clc; fprintf('%d/%d - (sess %d/%d)\n',zz,2,q,size(ZA{zz},1)); 
        X1 = ZA{zz}{q,1};
        X2 = ZA{zz}{q,3};
        Xss = cell(1,size(ZA{zz},2)-3); %cell(1,3);%size(ZA{zz},2)-3); 
        for j = 1:length(Xss)
            Xss{j} = ZA{zz}{q,j+3}; 
        end
        
        dU = ceil(EstDimns(zz,1)); 
        dS = ceil(EstDimns(zz,2)); 
        
%         [dU,dS] = deal(max([dU,dS])); 
        
        X1 = cellfun(@(x) x(:,1:dU),X1,'uni',0); 
        X2 = cellfun(@(x) x(:,1:dS),X2,'uni',0); 
        X3 = cell(size(X1)); 
        for j = 1:length(X1)
            X3{j} = [X1{j},X2{j}]; 
        end

        X1 = cellfun(@(x) x-repmat(nanmean(x(1:10,:)),size(x,1),1),X1,'uni',0); 
        X2 = cellfun(@(x) x-repmat(nanmean(x(1:10,:)),size(x,1),1),X2,'uni',0); 
        X3 = cellfun(@(x) x-repmat(nanmean(x(1:10,:)),size(x,1),1),X3,'uni',0); 
        for j = 1:length(Xss)
            Xss{j} = cellfun(@(x) x-repmat(nanmean(x(1:10,:)),size(x,1),1),Xss{j},'uni',0); 
        end
        
%         X1cat = cell2mat(X1); 
%         X2cat = cell2mat(X2); 
%         g = ~any(isnan([X1cat, X2cat]),2); 
%         [A,B,r,U,V] = canoncorr(X1cat(g,:),X2cat(g,:)); 
        mind = min([dU,dS]); 
        g = ~any(isnan([ZOav{zz},ZOSHav{zz}]),2); 
        X = ZOav{zz}(g,1:dU); X = X - nanmean(X); 
        Y = ZOSHav{zz}(g,1:dS); Y = Y - nanmean(Y); 
%         
%         [A,B,~,U,V] = canoncorr(X,Y); 
        
        bU = X\Y; 
        bS = Y\X; 
        
        [pnU,~,~] = svd(bU); 
        potU = pnU(:,1:size(bU,2)); 
        nullU = pnU(:,(size(bU,2)+1):end); 
%         [~,~,v1] = svd(U');
%         vre = v1(:,(mind+1):end)*v1(:,(mind+1):end)'; 
%         
%         Xn = (X'*(v1(:,(mind+1):end)*v1(:,(mind+1):end)'))'; 
%         
%         
%         [u1,s1,uv1] = svd(Xn'); 
%         
%         v1null = (v1(:,(mind+1):end)*v1(:,(mind+1):end)')u1(:,1:mind); 
%         
%         [~,~,v2] = svd(V');
%         Yn = (ZOSHav{zz}(g,1:dS)'*(v(:,(mind+1):end)*v(:,(mind+1):end)'))'; 
%         B1 = X1cat(g,:)\X2cat(g,:); 
%         B2 = X2cat(g,:)\X1cat(g,:); 
        lxss = length(Xss); 
        for j = 1:length(X1)
            Xss{lxss+1}{j,:} = X1{j}*bU; 
%             Xss{lxss+2}{j,:} = X1{j} - X2{j}*bS;
%             Xss{lxss+3}{j,:} = X2{j} - X1{j}*bU;
        end
        

        Y = FA{zz}{q}; Y = cellfun(@(x) x',Y,'uni',0); 
        Y = cellfun(@(x) x-nanmean(x(1:10)),Y,'uni',0); 
        
        YU = FAU{zz}{q}; YU = cellfun(@(x) x',YU,'uni',0); 
        YU = cellfun(@(x) x-nanmean(x(1:10)),YU,'uni',0); 
        
        YS = FAS{zz}{q}; YS = cellfun(@(x) x',YS,'uni',0); 
        YS = cellfun(@(x) x-nanmean(x(1:10)),YS,'uni',0); 
        
        for i = 1:length(X1)
            nani = any(isnan(X1{i}),2); 
            X1{i}(nani,:) = []; 
            X2{i}(nani,:) = []; 
            X3{i}(nani,:) = []; 
            for j = 1:length(Xss)
                Xss{j}{i}(nani,:) = []; 
            end
            Y{i}(nani) = []; 
            YU{i}(nani) = []; 
            YS{i}(nani) = []; 
        end
        emptied = cellfun(@(x) isempty(x),Y); 
        X1(emptied) = []; 
        X2(emptied) = []; 
        X3(emptied) = []; 
        Y(emptied) = []; 
        YU(emptied) = []; 
        YS(emptied) = []; 
        for j = 1:length(Xss); Xss{j}(emptied) = []; end
        tizg = tiz; tizg(emptied) =[]; 
        [~,tizo] = sortrows(tizg'); 
        tiz = 1:length(Y); tiz(tizo) = 1:length(Y); 
              
        X1 = cellfun(@(x) smooth_pad(x,10),X1,'uni',0); 
        X2 = cellfun(@(x) smooth_pad(x,10),X2,'uni',0); 
        X3 = cellfun(@(x) smooth_pad(x,10),X3,'uni',0); 
        for j = 1:length(Xss) 
            Xss{j} = cellfun(@(x) smooth_pad(x,10),Xss{j},'uni',0); 
            if size(Xss{j}{1},2)>size(Xss{j}{1},1) 
                Xss{j} = cellfun(@(x) x(:),Xss{j},'uni',0); 
            end
        end
        
        Y = cellfun(@(x) smooth_pad(x,10)',Y,'uni',0); 
        YU = cellfun(@(x) smooth_pad(x,10)',YU,'uni',0); 
        YS = cellfun(@(x) smooth_pad(x,10)',YS,'uni',0); 
        
        [HoldY,HoldX1,HoldX2,HoldX3,HoldReg] = deal(cell(size(Y))); 
        for i = 1:length(Y)
            thresh = min(Y{i})+.25*range(Y{i}); 
            oni = find(Y{i}>thresh); 
            mid12 = round(prctile(oni,[25 75])); 
            midi = mid12(1):mid12(end); 
            
            HoldReg{i} = false(1,length(Y{i})); 
            HoldReg{i}(midi) = true; 

            HoldY{i} = nanmean(Y{i}(midi)); 
            HoldX1{i} = nanmean(X1{i}(midi,:)); 
            HoldX2{i} = nanmean(X2{i}(midi,:)); 
            HoldX3{i} = nanmean(X3{i}(midi,:)); 
        end
        
        tpredSS = cell(length(Xss),1);
        [tpredU,tpredS,tpredHU,tpredHS,Yq,YqU,YqS,tpredCOMB,tpredHC] = deal(cell(length(X1),1)); 
        for i = 1:length(X1)
            testi = i; 
            traini = 1:length(X1); traini(i) = []; 

            trainX1 = cell2mat(X1(traini)); 
            trainX2 = cell2mat(X2(traini)); 
            trainX3 = cell2mat(X3(traini)); 
            trainXss = cell(1,length(Xss)); 
            for j = 1:length(Xss)
                trainXss{j} = cell2mat(Xss{j}(traini)); 
            end
            trainY = cell2mat(Y(traini)); 
            trainYU = cell2mat(YU(traini));
            trainYS = cell2mat(YS(traini));

            b1 = [ones(length(trainYU),1), trainX1]\trainYU; 
            b2 = [ones(length(trainYS),1), trainX2]\trainYS; 
            b3 = [ones(length(trainY),1), trainX3]\trainY; 
            bss = cell(length(Xss),1); 
            for j = 1:length(Xss)
                bss{j} = [ones(length(trainY),1), trainXss{j}]\trainY; 
            end
            
            pred1 = [ones(length(trainYU),1), trainX1]*b1; 
            pred2 = [ones(length(trainYS),1), trainX2]*b2; 
            pred3 = [ones(length(trainY),1), trainX3]*b3; 
            predss = cell(length(Xss),1); 
            for j = 1:length(predss)
                predss{j} = [ones(length(trainY),1), trainXss{j}]*bss{j}; 
            end
            
            bb1 = polyfit(pred1,trainYU,3); 
            bb2 = polyfit(pred2,trainYS,3); 
            bb3 = polyfit(pred3,trainY,3); 
            bbss = cell(length(Xss),1); 
            for j = 1:length(Xss)
                bbss{j} = polyfit(predss{j},trainY,3);  
            end
            
            p1 = [ones(size(X1{testi},1),1), X1{testi}]*b1;
            p2 = [ones(size(X2{testi},1),1), X2{testi}]*b2; 
            p3 = [ones(size(X3{testi},1),1), X3{testi}]*b3; 
            pss = cell(length(Xss),1); 
            for j = 1:length(Xss)
                pss{j} = [ones(size(Xss{j}{testi},1),1), Xss{j}{testi}]*bss{j};
            end
             
            tpredU{tiz(i)} = [polyval(bb1,p1);NaN(50,size(pred1,2))]; 
            tpredS{tiz(i)} = [polyval(bb2,p2);NaN(50,size(pred2,2))]; 
            tpredCOMB{tiz(i)} = [polyval(bb3,p3);NaN(50,size(pred3,2))]; 
             
            for j = 1:length(Xss)
                tpredSS{j}{tiz(i),:} = [polyval(bbss{j},pss{j});NaN(50,size(predss{j},2))]; 
                tpredSS{j}{tiz(i)}(tpredSS{j}{tiz(i)}<0) = 0; 
            end
            Yq{tiz(i)} = [Y{testi};NaN(50,1)]; 
            YqU{tiz(i)} = [YU{testi};NaN(50,1)]; 
            YqS{tiz(i)} = [YS{testi};NaN(50,1)]; 
           
%             
%             tpredU{i} = [ones(size(X1{testi},1),1), X1{testi}]*b1; 
%             tpredS{i} = [ones(size(X2{testi},1),1), X2{testi}]*b2; 
            
            tpredU{tiz(i)}(tpredU{tiz(i)}<0) = 0; 
            tpredS{tiz(i)}(tpredS{tiz(i)}<0) = 0; 
            tpredCOMB{tiz(i)}(tpredCOMB{tiz(i)}<0) = 0; 

        end
        Upreds{zz}{q,:} = tpredU; 
        Spreds{zz}{q,:} = tpredS; 
        for j = 1:length(tpredSS)
            SSpreds{zz}{q,j} = tpredSS{j}; 
        end
        Yacts{zz}{q,:} = Yq; 
        YUacts{zz}{q,:} = YqU; 
        YSacts{zz}{q,:} = YqS; 

        Yreg{zz}{q} = cell2mat(Yq); 
        ULpred{zz}{q} = cell2mat(tpredU);   
        SLpred{zz}{q} = cell2mat(tpredS);
        CLpred{zz}{q} = cell2mat(tpredCOMB); 
             
    end

end

r2_uss = cell(length(Yacts),1); 
for z = 1:length(Yacts)
    
    Ucell = vertcat(Upreds{z}{:}); 
    Scell = vertcat(Spreds{z}{:});
    [SScell,ssc] = deal(cell(1,size(SSpreds{z},2))); 
    for i = 1:size(SSpreds{z},2) 
        SScell{i} = vertcat(SSpreds{z}{:,i});
    end
    Ycell = vertcat(Yacts{z}{:}); 
    YcellU = vertcat(YUacts{z}{:}); 
    YcellS = vertcat(YSacts{z}{:}); 
    
    [~,bootsam] = bootstrp(nboot,[],Ucell); 
    for i = 1:size(bootsam,2)
        clc; fprintf('%d/%d (%d/%d)\n',z,length(Yacts),i,size(bootsam,2)); 
        uc = cell2mat(Ucell(bootsam(:,i))); 
        sc = cell2mat(Scell(bootsam(:,i))); 
        yc = cell2mat(Ycell(bootsam(:,i))); 
        ycU = cell2mat(YcellU(bootsam(:,i))); 
        ycS = cell2mat(YcellS(bootsam(:,i))); 
        
        r2_uss{z}(i,1:2) = [R2fun(ycU,uc),R2fun(ycS,sc)];%,R2fun(yc,ssc)];
        for j = 1:length(SScell)
            ssc{j} = cell2mat(SScell{j}(bootsam(:,i)));
        
            r2_uss{z}(i,j+2) = R2fun(yc,ssc{j});
        end
        
    end 
end
means = cellfun(@(x) mean(x),r2_uss,'uni',0); 
bnds = cellfun(@(x) prctile(x,[2.5 97.5]),r2_uss,'uni',0); 
%%

[pvar_UY,ts] = deal(cell(2,1)); 

figure('Position',[69 673 1787 213]); hold on; 
subplot(1,6,1:5); hold on; 
plot(Yreg{1}{5},'k'); 
plot(ULpred{1}{5},'Color',clrs(1,:)); 
plot(SLpred{1}{5},'Color',clrs(3,:)); 
ah = gca; clean_plot(ah,'timebar'); 

subplot(1,6,6); hold on; 
for i = 1:2
    
    plot(repmat(1:2,2,1)+.2*i-.3,bnds{i}(:,[1 2]),'Color',[.8 .8 .8],'LineWidth',3); 
    plot((1:2)+.2*i-.3,means{i}([1 2]),'.','Color',clrs(i,:)); 

    X = ZOav{i}(:,1:EstDimns(i,1))*UY{i};
    Y = ZOSHav{i}(:,1:EstDimns(i,2));
    g = ~any(isnan([X,Y]),2); 
    for j = 1:size(X,2) 
        b = [ones(sum(g),1),Y(g,:)]\X(g,j:end); 
        rec = [ones(size(Y,1),1) Y]*b; 
        pvar_UY{i}(j) = sum(nanvar(rec))./sum(nanvar(X(:,j:end)))*100; 
    end
end
set(gca,'XTick',1:2); 
xticklabels({'Action-unique','Shared'}); 
xlim([.5 2.5]); ylim([.2 1]); 
clean_plot; 
% 
figure('Position',[342 401 1306 496]); hold on; 
subplot(1,2,1); hold on;
for i = 1:2
    plot(repmat(1:3,2,1)+.2*i-.3,bnds{i}(:,[1 2 end]),'Color',[.8 .8 .8],'LineWidth',3); 
    plot((1:3)+.2*i-.3,means{i}([1 2 end]),'.','Color',clrs(i,:)); 
end
ylim([.2 .8]); 
set(gca,'XTick',[1 2 3]); 
ylabel('Decoding R2'); 
xticklabels({'Action-unique','Shared','Action-unique (shared dyn)'}); 
clean_plot

subplot(1,2,2); hold on; 
plot(repmat(pvar_UY{1},2,1),bnds{1}(:,3:(end-1)),'Color',[.8 .8 .8],'LineWidth',3); 
plot(repmat(pvar_UY{2},2,1),bnds{2}(:,3:(end-1)),'Color',[.8 .8 .8],'LineWidth',3); 
plot(pvar_UY{1},means{1}(3:(end-1)),'.-','Color',clrs(1,:));
plot(pvar_UY{2},means{2}(3:(end-1)),'.-','Color',clrs(2,:)); 
xlabel('% Var explained by Shared'); 
ylabel('Decoding R2');
ylim([.2 .8]); 
clean_plot
% 
