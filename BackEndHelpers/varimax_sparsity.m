function [Q,Xsparse,sparse_VE,baselines] = varimax_sparsity(X,baseline_indices)
% Perform varimax rotation on baseline-subtracted activity
%
%%%%%%%% INPUTS %%%%%%%%
%
%   X: Data matrix (N x D), where N is the number of samples and D is the
%      dimensionality of the input. This works best on trial-averaged
%      responses
%
%   baseline_indices: 
%                    numeric (e.g. 1:10): Indices corresponding to "baseline", which is usually
%                    pre-trial or inter-trial activity. 
%                    Leave blank to automatically fit baselines or set to 'none' for no
%                    baseline subtraction.
%
%%%%%%%% OUTPUTS %%%%%%%%
%
%   Q: Varimax rotation matrix (D x D)
%
%   Xsparse: Varimax-rotated data
%
%   sparse_VE: % Variance explained by each dimension
%
%   baselines: baseline values

if nargin > 1 && ~isempty(baseline_indices) && isnumeric(baseline_indices)
    bi = find(~any(isnan(X),2),1,'first')+baseline_indices-1; 
    baselines = nanmean(X(bi,:),1); 
    autobase = false; 

elseif nargin > 1 && ischar(baseline_indices) && strcmp(baseline_indices,'none')
    baselines = zeros(1,size(X,2)); 
    autobase = false; 
else
    autobase = true;
end
bads = any(isnan(X),2); 
Xsparse = NaN(size(X)); 
X(bads,:) = []; 

if autobase
    clc; fprintf('finding baselines...\n'); 
    mincost = @(x) -varimax_cost(rotatefactors(X-x,'Method','varimax','Maxit',1000)); 

    % baselines cannot exist outside of the data range
    lb = min(X); 
    ub = max(X); 
    
    options = optimoptions(@fmincon,'Display','off'); 
    baselines = fmincon(mincost,mean(X),[],[],[],[],lb,ub,[],options); 
    
    fprintf('finding Varimax rotation...\n'); 
    if size(X,2)==1
        Qf = 1; 
    else
        Xb = X-baselines; 
        [~,Qf] = rotatefactors(Xb(var(Xb,[],2)>0,:),'Method','varimax','Maxit',10000,'Normalize','on'); 
    end
else
    if size(X,2)==1
        Qf = 1; 
    else
        Xb = X-baselines; 
        [~,Qf] = rotatefactors(Xb(var(Xb,[],2)>0,:),'Method','varimax','Maxit',10000,'Normalize','on'); 
    end
end

% Flip to positive
Xprime = (X-baselines)*Qf; 
signs = sign(mean([min(Xprime);max(Xprime)])); 
% signs = sign(median(Xprime)); 
Qf = Qf*diag(signs); 

[~,vrank] = sortrows(nanvar(Xprime)','descend'); 
Q = Qf(:,vrank); 

Xsparse(~bads,:) = (X-baselines)*Q; 
sparse_VE = nanvar(Xsparse)./sum(nanvar(Xsparse))*100; 

end

function vcost = varimax_cost(A)
    % standardize
    h = sqrt(sum(A.^2, 2)); 
    A = bsxfun(@rdivide, A, h); 
    p = size(A,1); 
    % varimax cost function
    vcost = 1/p*sum(sum(A.^4,1)- 1/p*sum(A.^2).^2); 
end
