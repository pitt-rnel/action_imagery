function [Subspaces,out] = DySO(Xcell,var_cutoff,varargin)
% Splits data into orthogonal subspaces containing condition-unique and 
% condition-shared activity.
%   [Subspaces] = DySO(Xcell,var_cutoff)
%
%- INPUTS -----------------------------------------------------------------
%
%   Xcell:      Cell array of matrices containing neural firing rates for
%               different conditions. In each cell, rows correspond to 
%               samples and columns correspond to neural dimensions. All 
%               cells must have the same dimensionality. It's recommended 
%               to provide activity with reduced dimensionality, but it is 
%               preferable that it be an overestimation of the true 
%               dimensionality. 
%
%   var_cutoff: Percent variance cutoff used to delineate potent and null
%               spaces. If no cutoff is given, the default value is 99.
%
%  optional inputs:    
%               'do_plot': flag to return plot (true) [default] or not (false)
%               
%               'combinations': custom list of combinations to check.
%                     Useful if there are many conditions and full
%                     enumeration of all combinations would be
%                     computationally infeasible. 
%
%- OUTPUTS ----------------------------------------------------------------
%
%   Subspaces:  Struct containing 'unique' and 'shared' fields, which
%               contain the axes for each identified subspace. Together,
%               they form a full orthonormal basis, i.e.
%                  
%                   >> Qfull = [cell2mat(Subspaces.unique.all), Subspaces.shared];
%                   >> Qfull'*Qfull
%                   ans = 
%                       1.00    0.00    0.00    0.00    ...
%                       0.00    1.00    0.00    0.00    ...
%                       0.00    0.00    1.00    0.00    ...
%                       0.00    0.00    0.00    1.00    ...
%                       ...
%                   >> norm(Qfull)  
%                   ans =
%                       1.000
%
%   out:        Struct containing projections of original data into identified
%               subspaces
%     
%- Description of method --------------------------------------------------
%   DySO (Dynamic Subspace Overlap) is based on the manifold 
%   view of neural activity, in which we can represent population activity 
%   as a point (state) within a low dimensional manifold. For different 
%   tasks, the neural state may inhabit different subspaces within the 
%   manifold. If the neural state projects in part onto an axis or subspace 
%   *only* during condition A, then we consider that axis/subspace to be
%   "A-unique". Likewise, if the neural state contains projections onto a 
%   different axis/subspace during multiple conditions (e.g. A, B, and C),
%   we consider it to be "ABC-shared".
%
%   The intuition behind the present method is as follows: 
%   Consider three conditions, A, B and C. If there is some axis/subspace
%   that is unique to A, then it must exist in the null space of B and C.
%   Likewise, a B-unique subspace exists in the null space of A and C, etc. .
%   We use PCA and a variance cutoff to identify the relevant null spaces 
%   for identifying unique activity: 
%
%       A_unique exists in {B,C}_null, B_unique in {A,C}_null, and C_unique
%       in {A,B} null
%
%   Doing so gives us the **form** of the unique activity, but since each
%   condition is calculated separately, they will not be mathematically
%   orthogonal (each basis estimation will incorporate some degree of 
%   noise, largely due to low-variance dimensions in the original data). 
%
%   However, once we have identified the profiles of the unique activity,
%   we can find an orthonormal transformation of the original space that
%   reconstructs those known profiles. We do this by minimizing sum 
%   squared error using the Manopt toolbox. The output of the optimization
%   reflects the combination of unique spaces, and we can simply take the
%   null space of this to define the shared space. 
%
%   ***The end result is an
%   orthonormal transformation of the original basis such that each axis
%   reflects either condition-unique or condition-shared activity***
%   
%
%   ----
%   created by:
%   Brian Dekleva
%   5/2/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% meansub = @(x) x - repmat(nanmean(x),size(x,1),1); 

Nc = length(Xcell); % Number of conditions
Nd = size(Xcell{1},2); % Data dimensionality
if nargin < 2 || isempty(var_cutoff)
    var_cutoff = 99; % Default % variance cutoff if none provided
end

plot_toggle_loc = find(cellfun(@(x) all(strcmp(x,'do_plot')),varargin))+1; 
custom_combs_loc = find(cellfun(@(x) all(strcmp(x,'combinations')),varargin))+1; 
if isempty(custom_combs_loc)
    list_type = 'single'; 
elseif ischar(varargin{custom_combs_loc})
    list_type = varargin{custom_combs_loc}; 
else
    list_type = 'custom'; 
end


if isempty(plot_toggle_loc)
    do_plot = true; 
else
    do_plot = varargin{plot_toggle_loc}; 
end


switch list_type
    case 'single' % Default
        combi = {(1:length(Xcell))'}; 
        comb_lens = 1; 
    case 'full' 
        comb_lens = Nc-1; % Number of combination lengths
        combi = cell(comb_lens,1); 
        for i = 1:comb_lens
            combi{i} = nchoosek(1:Nc,i); % Enumerate all combinations of length i
        end
    case 'custom'
        if iscell(varargin{custom_combs_loc})
            combi = varargin{custom_combs_loc}; 
            comb_lens = length(combi);
        else
            combi = varargin(custom_combs_loc); 
            comb_lens = 1;
        end
end

% Check inputs and remove NaNs
assert(length(unique(cellfun(@(x) size(x,2),Xcell)))==1,'All data must have same dimensionality'); 
Xcell = Xcell(:); 
for i = 1:length(Xcell)
    Xcell{i} = Xcell{i}(~any(isnan(Xcell{i}),2),:); 
end

% Concatenate all conditions 
Z = cell2mat(Xcell); %cell2mat(cellfun(@(x) meansub(x),Xcell,'uni',0)); 

% Set up optimization params
options.verbosity = 0; 
options.maxiter = 1000; 
problem = []; 
warning('off', 'manopt:getHessian:approx');

[Q_COMB,Q_CAT,C_inds,XV] = deal(cell(1,comb_lens)); 
Nuller = cell(1,comb_lens+1); Nuller{1} = eye(Nd); 
for z = 1:comb_lens % Loop through combination lengths
    
    Zz = Z*Nuller{z}; % Only check for unique dimensions not captured by earlier combination lengths

    C = size(combi{z},1); % number of combinations of length z
    C_inds{z} = cell(1,C); 
    
    nullthresh = (100-var_cutoff); 
    cflag = false; 
    
    while ~cflag
%         figure; hold on; clrs = lines(10); 
        [potXnoi,nullXnoi,uniqueX] = deal(cell(C,1)); 
        for i = 1:C % Loop through combinations of length z

            % determine dimensionality of the set of conditions not containing i
            noi = 1:Nc; noi(ismember(noi,combi{z}(i,:))) = [];   

            [pcnoi,~,~,~,evnoi] = pca(cell2mat(Xcell(noi))); 
            dimnoi = find(cumsum(evnoi)>var_cutoff,1,'first'); 
            var_in_null = sum(var(cell2mat(Xcell(noi))*pcnoi(:,(dimnoi+1:end)))); 

            [pcnoi] = pca(cell2mat(Xcell(noi))*Nuller{z}); 
            dimnoi = find(cumsum(flipud(diag(cov(cell2mat(Xcell(noi))*Nuller{z}*pcnoi))))<=var_in_null,1,'last'); % potent dimensionality of ~i
            potXnoi{i} = pcnoi(:,1:(end-dimnoi)); 
            nullXnoi{i} = pcnoi(:,(end-dimnoi+1):end); % null space of ~i

            [evXs,dimUis] = deal(cell(1,size(combi{z},2))); 
            [dimX] = deal(zeros(1,size(combi{z},2))); 
            for j = 1:size(combi{z},2)
                % determine dimensionality of condition(s) in combi{i}
%                 [~,~,~,~,evXs{j}] = pca(Xcell{combi{z}(i,j)});
                pcj = pca(Xcell{combi{z}(i,j)});
                evXs{j} = (var(Xcell{combi{z}(i,j)}*pcj)./sum(var(Xcell{combi{z}(i,j)}))*100)'; 
                dimX(j) = find(cumsum(evXs{j})>var_cutoff,1,'first'); 

                % Activity from condition(s) i projected into the null space of condition(s) 
                % ~i represents i-unique
                i_in_noi = Xcell{combi{z}(i,j)}*Nuller{z}*nullXnoi{i}; 
                
%                 pcnoi = pca(i_in_noi); 
%                 evi = Xcell{combi{z}(i,j)}*Nuller{z}*[potXnoi{i},nullXnoi{i}*pcnoi]; evi = (var(evi)./sum(var(evi))*100)'; 
%                 subplot(C,1,i); hold on; 
%                 b = bar([evnoi,evi],'EdgeColor','none'); yl = ylim; 
%                 plot(size(potXnoi{i},2)*[1 1]+.5,yl,'k:'); 
%                 b(1).FaceColor = [.5 .5 .5];
%                 b(2).FaceColor = clrs(i,:); 
%                 b(1).CData = repmat([.5 .5 .5],20,1); 
%                 b(2).CData = repmat(clrs(i,:),20,1); 
%                 make_pretty; axis off
                
                unique_i = pca(i_in_noi); % do pca on projected data
                uniqueXvar = var(i_in_noi*unique_i)./trace(cov(Xcell{combi{z}(i,j)}))*100; % explained condition var
                % Discard dimensions that explain % variance below threshold
                dimUis{j} = sum(~(cumsum(fliplr(uniqueXvar)) < nullthresh)); 

            end

            % only consider to be unique to the set if it applies to *all*
            % conditions within the set
            num_unique_to_all = min(cell2mat(dimUis)); 

            if num_unique_to_all>0

                i_in_noi = cell2mat(Xcell(combi{z}(i,:)))*Nuller{z}*nullXnoi{i}; % concatenated i conds in null of ~i
                unique_i = pca(i_in_noi); % do pca on projected data

                uniqueX{i} = nullXnoi{i}*unique_i(:,1:num_unique_to_all); % estimated low-d condition-unique activity

                C_inds{z}{i} = repmat(i,1,num_unique_to_all); % index of the combination for later record-keeping 
                XV{z}{i} = uniqueXvar(1:dimUis{j}); 
            else
                uniqueX{i} = zeros(Nd,0); % No unique, so create empty space-saver
                XV{z}{i} = []; 
            end
        end

        D_unique = cellfun(@(x) size(x,2),uniqueX); % dimensionalities of unique spaces

        cflag = sum(D_unique)<=size(Nuller{z},2); 
        nullthresh = nullthresh*1.1; 
    end
    
    uniqueI = find(D_unique>0); 
    D_unique = D_unique(uniqueI); uniqueX = uniqueX(uniqueI); 
    dimU = length(D_unique); % Number of unique condition sets of length z
    dimUtot = sum(D_unique); % Total dimensionality of all unique spaces
    
    if dimU == 0 
        fprintf('No unique spaces found of length %d\n',size(combi{z},2)); 
        Nuller{z+1} = Nuller{z}; 
        C_inds{z} = []; 
        continue
    end

%     Z_u_v = cellfun(@(x) nanvar(Zz*x)./sum(nanvar(Zz))*100,uniqueX,'uni',0); 
    
    % Get projection of all conditions into each unique space
    Z_uniques = cellfun(@(x) Zz*x,uniqueX,'uni',0); 
    

    problem.M = stiefelfactory(size(Zz,2),dimUtot); 
    problem.cost = @(Q) sum((Zz*Q - cell2mat(Z_uniques')).^2,'all');
    problem.egrad = @(Q) 2*Zz'*(Zz*Q - cell2mat(Z_uniques'));
    
    % If you want, you can check gradient with: 
    % >> figure; checkgradient(problem); 

    % Solve for the COMBINED orthonormal subspace QZ = [A_unique,B_unique,C_uniques,...]. 
    % This ensures that each resulting subspace will be orthonormal as well as orthogonal
    % to each other. 
    fprintf('fitting unique spaces of length %d...\n',z); 
    QZ = trustregions(problem,[],options); 
    
    % Flip to positive (just for easier interpretation/visualization)
    QZ = flip_positive(QZ,C_inds{z},Xcell,Nuller{z},combi{z},Nc); 

    % Create Usplit to recover activity from each condition-unique subspace.
    % i.e.
    %   Usplit{1} = [ I ; 0; 0 ]
    %   Usplit{2} = [ 0 ; I; 0 ]
    %   ...
    Usplit = getUsplit(D_unique); 
    for i = 1:length(Usplit)
        % We can now use U_split to back out the separate unique spaces, 
        % and return to the original dimensionality 
        Q_COMB{z}{i} = Nuller{z}*QZ*Usplit{i};
    end
    Q_CAT{z} = Nuller{z}*QZ; % Concatenated version 
    
    Nuller{z+1} = null(cell2mat(Q_CAT(1:z))'); % Update 'previous-unique'-null transformation

    % Remove indices of unused combinations
    C_inds{z} = cellfun(@(x) x(~isnan(x)),C_inds{z},'uni',0); 
    C_inds{z}(cellfun(@(x) isempty(x),C_inds{z})) = []; 
end 
QZ = cell2mat(Q_CAT); 

% Fill out output struct
Subspaces.unique.all = Q_CAT; 

for i = 1:length(C_inds)
    for j = 1:length(C_inds{i})
        combname = ['C' num2str(combi{i}(C_inds{i}{j}(1),:))]; 
        combname = strrep(combname,'  ','_'); 
        Subspaces.unique.(combname) = Q_COMB{i}{j}; 
        
        out.proj.unique.(combname) = cellfun(@(x) x*Q_COMB{i}{j},Xcell','uni',0); 
    end
end

% The shared space is defined as the null space of QZ. 
Subspaces.shared = null(QZ'); 
pcsh = pca(cell2mat(Xcell)*Subspaces.shared);%
Subspaces.shared = Subspaces.shared*pcsh; 
Subspaces.full = [cell2mat(Subspaces.unique.all),Subspaces.shared]; 

out.proj.shared = cellfun(@(x) x*Subspaces.shared,Xcell','uni',0); 

unique_fields = fieldnames(Subspaces.unique); unique_fields(1) = []; % remove 'all' field
out.proj.split = cell(length(Xcell),length(unique_fields)+1); 
for i = 1:length(unique_fields)
    out.proj.split(:,i) = cellfun(@(x) x*Subspaces.unique.(unique_fields{i}),Xcell,'uni',0); 
end
out.proj.split(:,end) = cellfun(@(x) x*Subspaces.shared,Xcell,'uni',0); 
fprintf('Done\n'); 

if do_plot
    %% For plotting
    
    Nu = length(unique_fields);
%     lowD_unique = cell(Nu,1); 
%     for i = 1:Nu
%         [pcI,~,~,~,varIunique] = pca(cell2mat(Xcell)*Subspaces.unique.(unique_fields{i})); 
%         dI = find(cumsum(varIunique)>var_cutoff,1,'first'); 
%         lowD_unique{i} = pcI(:,1:dI); 
%     end
    [pcS,~,~,~,varS] = pca(cell2mat(Xcell)*Subspaces.shared); 
    dS = find(cumsum(varS)>var_cutoff,1,'first'); 
    lowD_shared = pcS(:,1:dS); 

    fs = factor(Nu+1); it = 1; while length(fs)==1; fs = factor(Nu+1+it); it = it+1; end
    sp_i = [(Nu+it)/fs(end) fs(end)]; 
    clrs = lines(50); 
    figure('Position',[207 97 1419 899]); hold on; 
    sp = cell(Nu+1,1); 
    for i = 1:Nu
        spname = unique_fields{i}; spname = strrep(spname,'_','|'); spname = strrep(spname,'C',''); 
        iv = sum(nanvar(cell2mat(Xcell)*Subspaces.unique.(unique_fields{i})))./sum(nanvar(cell2mat(Xcell)))*100; 
        sp{i} = subplot(sp_i(1),sp_i(2),i); hold on; title(sprintf('unique %s (%.1f%%)',spname,iv)); 
        for j = 1:Nc
            if ismember(j,cell2mat(cellfun(@str2num,num2cell(unique_fields{i}),'uni',0)))
                plot(Xcell{j}*Subspaces.unique.(unique_fields{i}),'Color',clrs(i,:),'LineWidth',2); 
            else
                plot(Xcell{j}*Subspaces.unique.(unique_fields{i}),'Color',[.5 .5 .5],'LineWidth',1); 
            end
        end
    end

    sp{end} = subplot(sp_i(1),sp_i(2),Nu+1); hold on; 
    legcell = cell(1,size(lowD_shared,2)); 
    for i = 1:Nc
        for j = 1:size(lowD_shared,2)
            plot(Xcell{i}*Subspaces.shared*lowD_shared(:,j),'Color',clrs(j,:),'LineWidth',2); 
            if i==1; legcell{j} = sprintf('shared %d',j); end
        end
        if i==1; leg = legend(legcell); end
    end
    iv = sum(nanvar(cell2mat(Xcell)*Subspaces.shared))./sum(nanvar(cell2mat(Xcell)))*100; 
    title(sprintf('shared (%.1f%%)',iv)); 

    allys = cell2mat(cellfun(@(x) x.YLim,sp,'uni',0)); 
    zl = [min(allys(:)), max(allys(:))]; 

    for i = 1:length(sp); sp{i}.YLim = zl; end 
    if exist('make_pretty','file'); make_pretty; end
    leg.String((length(legcell)+1):end) = []; 
    
    %%%% Bar plots
    V_conds = cell(length(Xcell),length(unique_fields)+1,1); 
    for i = 1:length(Xcell)
        cvar = sum(nanvar(Xcell{i})); 
        for j = 1:length(unique_fields)
            V_conds{i,j} = nanvar(Xcell{i}*Subspaces.unique.(unique_fields{j}))./cvar*100; 
            
        end
        V_conds{i,end} = nanvar(Xcell{i}*Subspaces.shared*pcS)./cvar*100; 
    end
    
    titlenames = [unique_fields', {'shared'}]; titlenames = strrep(titlenames,'_','|');
    titlenames = cellfun(@(x) strrep(x,'C','unique '),titlenames,'uni',0); 
%     sp = cell(1,length(unique_fields)+1); 
    figure('Position',[106 448 1660 402]); hold on; 
    bar(cell2mat(V_conds)'); yl = ylim; 
    splitloc = cumsum(cellfun(@(x) size(x,2),V_conds(1,1:end-1)))+.5; 
    labloc = nanmean([0 splitloc; splitloc, Nd]); 
    
    for i = 1:(length(titlenames)-1)
        plot(splitloc(i)*[1 1],yl,'k:'); 
    end
    for i = 1:length(titlenames)
        text(labloc(i),yl(2),titlenames{i},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',12); 
    end
%     
    make_pretty; xlim([0 Nd+1]);  g = gca; set(g,'XTick',1:Nd); 
    ylabel('% Variance'); xlabel('Dimensions');

end
end
% 
function Usplit = getUsplit(D_unique)
    dimU = length(D_unique); % Number of unique conditions
    Usplit = cell(dimU,1);

    Usplit{1} = zeros(sum(D_unique),D_unique(1)); 
    Usplit{1}(1:D_unique(1),:) = eye(D_unique(1)); 
    for i = 2:dimU
        Usplit{i} = zeros(sum(D_unique),D_unique(i)); 
        splitind = (sum(D_unique(1:(i-1)))+1):(sum(D_unique(1:(i-1)))+D_unique(i)); 
        Usplit{i}(splitind,:) = eye(D_unique(i)); 
    end
end

function QZ = flip_positive(QZ,C_inds,Xcell,Nuller,combi,Nc)
    ci = cell2mat(C_inds); ci(isnan(ci)) = []; 
    for i = 1:size(QZ,2)
        dprojI = cell2mat(Xcell(combi(ci(i),:)))*Nuller*QZ(:,i); 
        noi = 1:Nc; noi(combi(ci(i),:)) = []; 
        dprojnoI = cell2mat(Xcell(noi))*Nuller*QZ(:,i); 
        
        imax = dprojI(abs(dprojI)==max(abs(dprojI))); 
        noimax = dprojnoI(abs(dprojnoI)==max(abs(dprojnoI))); 
        if sign(imax-noimax)<0
            QZ(:,i) = -QZ(:,i);
        end
    end
end

function ndims = scree_dimensionality(scree)

    N = length(scree); 
    d2_line = @(x,y) norm(cross([1 scree(1) 0] - [length(scree) scree(end) 0],[x y 0] -...
        [length(scree) scree(end) 0])) / norm([1 scree(1) 0] - [length(scree) scree(end) 0]); 

    ds = NaN(N,1); 
    for k = 1:N
        ds(k) = d2_line(k,scree(k)); 
    end

    ndims = find(ds==max(ds),1,'first'); 
end