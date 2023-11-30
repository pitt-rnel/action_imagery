function [Uy,Bx,r,vaf,Zy,Zx,cost] = ReverseFit(Y,X,method)
switch method
    case 'VAF'
        Y = Y-nanmean(Y); 
        X = X-nanmean(X); 
end

assert(size(X,1)==size(Y,1),'X and Y must be the same length'); 

nY = size(Y,2); 
g = ~any(isnan([X,Y]),2); 

Yloop = Y(g,:); 
Xg = X(g,:); 
[nullerY,Qfit] = deal(cell(1,nY-1)); 
for q = 1:(nY-1) % dimensions

    clc; fprintf('Sorting from worst - (%d/%d)\n',q,nY); 
    [Qfit{q},~,~,costfun] = LeastAlignableComponent(Yloop,Xg,method);
 
    [U,~,~] = svd(Qfit{q}); 
    nullerY{q} = U(:,2:end); 
    
    Yloop = Yloop*nullerY{q}; 
end

%% Compile reduced space eigenvectors into full space
potfullY = NaN(nY,nY-1); % Initialize 
potfullY(:,1) = Qfit{1}; 
for i = 2:length(Qfit)
    nullerpre2 = eye(size(nullerY{1},1)); 
    for j = 1:(i-1)
        nullerpre2 = nullerpre2*nullerY{j}; 
    end
    potfullY(:,i) = nullerpre2*Qfit{i}; % Re-project
end    
[UlastY,~,~] = svd(potfullY); % SVD to bring matrix from N-1 to N dimensionality
Uy = [potfullY, UlastY(:,nY:end)]; 

Zy = Y * Uy; 
Zx = NaN(size(Zy)); 
Bx = NaN(size(X,2),size(Zy,2)); 
for i = 1:size(Zy,2)
    b = X(g,:)\Zy(g,i); 
    Zx(:,i) = X*b;
    Bx(:,i) = b; 
end

r = diag(corr(Zx(g,:),Zy(g,:))); 
vaf = 1 - nansum((Zx-Zy).^2)./nansum((Zy - repmat(nanmean(Zy),size(Zy,1),1)).^2); 
cost = NaN(size(r)); 
for i = 1:length(cost)
    cost(i) = costfun(Y(g,:),X(g,:),Uy(:,i)); 
end
clc; 
fprintf('Finished \nr = [')
fprintf(' %.2f ',r); 
fprintf(']\n'); 
fprintf('vaf = [')
fprintf(' %.2f ',vaf); 
fprintf(']\n')
end