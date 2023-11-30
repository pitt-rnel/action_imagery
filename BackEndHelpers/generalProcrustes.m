function[Zcell,transform,rs,meanshape,Xcellraw] = generalProcrustes(Xcell_in,baseline_inds,do_Varimax)
% generalProcrustes.m performs generalized Procrustes analysis on the
% D-dimensional traces in Xcell_in. Generalized Procrustes analysis is a
% method for aligning a single D-dimensional shape across multiple
% datasets. This function also performs an Varimax-based decomposition of the
% resulting aligned components. 
%
% Inputs:
%
%   Xcell_in     : S x 1 Cell array containing data matrices. Each cell
%                  should contain data matrix of size T x D. 
%
%                   e.g. for 5 datasets containing 10-dimensional traces of
%                   length 200: 
%                   
%                   >> Xcell_in
%
%                   Xcell_in = 
%
%                     5x1 cell array
%
%                       {200x10 double}
%                       {200x10 double}
%                       {200x10 double}           
%                       {200x10 double}    
%                       {200x10 double}
%                           
%   baseline_inds  : Indices corresponding to baseline activity in the
%                    data. This is used for Varimax rotation 
%                   (default = 1:10)
%
%   do_Varimax         : Flag to perform Varimax on the resulting
%                    aligned components. Setting to false will only return
%                    the standard Generalized Procrustes Analysis result
%                    (default = true)
% 
% Outputs:      
%   
%   Zcell          : Cell of same size as Xcell_in containing transformed
%                   (rotated) components
%
%   transform      : Cell array of structs containing the Procrustes
%                    transformations for each dataset. 
%                       transform{i}.T : Rotation matrix
%                       transform{i}.b : Isometric scaling factor
%                       transform{i}.c : offsets
%
%                    Zcell{i} can be reconstructed via:
%                    = Xcell_in{i}*transform{i}.T*transform{i}.b+transform{i}.c
%
%   rs             : Cell array containing Pearson correlation coefficients
%                    between a dataset's transformed components and the 
%                    Procrustes template (mean shape) 
%
%   meanshape      : Array of size Xcell_in{1} containing the mean shape
%                    obtained via Generalized Procrustes Analysis + ICA

% Check that inputs are correct 
assert(all(cellfun(@(x) size(x,1)==size(Xcell_in{1},1),Xcell_in)),'Each dataset must have same number of samples')
if ~all(cellfun(@(x) all(size(x)==size(Xcell_in{1})),Xcell_in))
    warning('Inputs are of different dimensionalities. Be careful'); 
end
% Set defaults if necessary
if nargin < 2
    baseline_inds = 1:10; 
    do_Varimax = true; 
elseif nargin < 3
    do_Varimax = true; 
end
% Initialize
XcellO = Xcell_in(:)'; 

% If different dimensionalities (they shouldn't be, but sometimes can be
% okay if you know what you're doing)
maxdim = max(cellfun(@(x) size(x,2),XcellO)); 
for i = 1:length(XcellO)
    if size(XcellO{i},2)<maxdim
        XcellO{i} = [XcellO{i},zeros(size(XcellO{i},1),maxdim-size(XcellO{i},2))]; 
    end
end

% Remove NaN and refigure baselines accordingly
bads = any(isnan(cell2mat(XcellO)),2); 
blines = false(1,length(bads)); blines(baseline_inds) = true; 
blines(bads) = []; blines = find(blines); 

% Remove bad times and mean subtract
XcellO = cellfun(@(x) x(~bads,:),XcellO,'uni',0); 
Xcell = cellfun(@(x) x-repmat(mean(x),size(x,1),1),XcellO,'uni',0); 

% Do initial Procrustes alignment to find best starting dataset
dinit = zeros(length(Xcell)); 
for i = 1:length(Xcell)
    for j = 1:length(Xcell)
        [dinit(i,j),~] = procrustes(Xcell{i},Xcell{j}); 
    end
end
bestinit = find(mean(dinit)==min(mean(dinit)),1,'first'); 

% Initialize
refshape = Xcell{bestinit}; % Starting reference shape
startshape = refshape; 
[Z] = deal(cell(1,length(Xcell))); 
[prods] = deal(NaN(1,length(Xcell))); 
d_prev = 1e4; % Stopping criteria
d_change =  inf;
while d_change > 1e-10 

    % Procrustes alignment to reference shape
    for i = 1:length(Xcell)
        [prods(i),Z{i}] = procrustes(refshape,Xcell{i}); 
    end

    % Compute mean shape
    meanshape = nanmean(cell2mat(reshape(Z,1,1,[])),3); 

    % Procrustes distance between mean shape and reference shape
    d = procrustes(refshape,meanshape); 

    refshape = meanshape./sum(var(meanshape)); % Set reference shape to mean shape
   
    % Change in procrustes distance
    d_change = d_prev - d; 
    d_prev = d; 
end
% the procedure can lead to drastic changes in overall magnitudes, so
% rescale the reference to be of equal magnitude as the original 
[~,~,tStargettoFin] = procrustes(startshape,meanshape); 
meanshape = meanshape*tStargettoFin.b; 

%% If doing Varimax
[Zcell,transform,rs] = deal(cell(length(Xcell),1)); 
if do_Varimax
    [~,shape] = varimax_sparsity(meanshape,blines); 
else
    shape = meanshape - repmat(nanmean(meanshape(blines,:)),size(meanshape,1),1);
    [pcshape] = pca(shape); 
    shape = (meanshape - nanmean(meanshape))*pcshape; 
end
% 
% rorig = NaN(length(Xcell),size(shape,2)); 
% for i = 1:length(Xcell)
%     [~,z] = procrustes(shape,XcellO{i}); 
%     rorig(i,:) = diag(corr(z,shape));
%     
% end
% [~,r_ords] = sortrows(nanmean(rorig)','descend'); 
% shape = shape(:,r_ords); 
% meanshape = meanshape(:,r_ords); 
%% Loop through datasets and find new Procrustes alignments
Xcellraw = cell(size(XcellO)); 
for i = 1:length(Xcell)
    Zcell{i} = NaN(length(bads),size(Xcell{i},2)); 
    [~,Zcell{i}(~bads,:),tb] = procrustes(shape,XcellO{i}); 

    extrabls = nanmean(Zcell{i}(baseline_inds,:)); 

    Zcell{i} = Zcell{i} - repmat(extrabls,size(Zcell{i},1),1); 

    rs{i} = diag(corr(Zcell{i}(~bads,:),shape))'; 

    transform{i} = tb; 
    transform{i}.c = NaN(length(bads),size(XcellO{i},2)); 
    transform{i}.c(~bads,:) = tb.c - repmat(extrabls,size(tb.c,1),1); 
    
    Xcellraw{i} = NaN(length(bads),size(XcellO{i},2)); 
    Xcellraw{i}(~bads,:) = XcellO{i}; 
end

meanshape = nanmean(cell2mat(reshape(Zcell,1,1,[])),3); 

end
