function [Xrot,rmetric] = Orthogonal_Alignment(X,Y,metric)
% Collection of methods for finding an orthogonal transformation of X to
% align to Y. The alignment can be set to maximize correlation or covariance,
% or to minimize the residuals. 

assert(size(X,1)==size(Y,1),'X and Y must be the same length'); 

g = ~any(isnan([X,Y]),2); 

Xo = X; Yo = Y; 
X = Xo(g,:)-nanmean(Xo,1); 
Y = Yo(g,:)-nanmean(Yo,1); 

options.verbosity = 0;
options.maxiter = 1000; 
warning('off', 'manopt:getHessian:approx');
options.maxiter = 1000; 
switch metric
    
    case 'Correlation'
%         [pc,~,~,~,~] = pca(X); 
%         % X = X*pc(:,ev > 2); 
%         X = X*pc(:,1:size(Y,2)); 

        nX = size(X,2);
        nY = size(Y,2); 

        problem.M = stiefelfactory(nX,nY); 
        problem.cost = @(Q) -trace((Y'*X*Q)./sqrt((Y'*Y).*((X*Q)'*(X*Q)))); 
        problem.egrad = @(Q) -(X'*Y*(eye(size(Y,2))./sqrt((Y'*Y).*((X*Q)'*X*Q))) - ...
            (0.5*(X'*X*Q)*(((X*Q)'*Y).*eye(size(Y,2))./sqrt((Y'*Y).*((X*Q)'*X*Q)).*(Y'*Y)./((Y'*Y).*((X*Q)'*X*Q))) + ...
            0.5*(X'*X*Q)*((Y'*X*Q).*eye(size(Y,2))./sqrt((Y'*Y).*((X*Q)'*X*Q)).*(Y'*Y)./((Y'*Y).*((X*Q)'*X*Q))))); 

        
%         b = X\Y; Qstart = b./norm(b); 
        
        [Q, ~, ~, ~] = trustregions(problem, [],options);

        Q = Q*diag(sign(diag(corr(Y,X*Q)))); 
%         Xrot = pc(:,1:size(Y,2))*Q; 
        Xrot = Q; 
        rmetric = diag(corr(X*Xrot,Y)); 
        
    case 'Covariance'
        
        problem.M = stiefelfactory(size(X,2),size(Y,2)); 
        problem.cost = @(Q) -sum(diag(Y'*(X*Q)).^2); 
        problem.egrad = @(Q) -2*X'*Y*diag(diag(Y'*(X*Q)));

        [Q, ~, ~, ~] = trustregions(problem, [],options);
        
        Q = Q*diag(sign(diag(corr(Y,X*Q)))); 
        Xrot = Q; 
        rmetric = diag((X*Xrot)'*Y); 
        
    case 'Residuals'
        
        problem.M = stiefelfactory(size(X,2),size(Y,2)); 
        problem.cost = @(Q) sum((X*Q - Y).^2,'all');
        problem.egrad = @(Q) 2*X'*(X*Q - Y); 
        
        [Xrot,~,~,~] = trustregions(problem,[],options); 
        rmetric = mean((X*Xrot - Y).^2); 
end
