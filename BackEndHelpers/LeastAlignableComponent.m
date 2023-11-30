function [Qf,r,best_b,costfun] = LeastAlignableComponent(Ytarget,Xin,method)
% 
if nargin < 3
    method = 'Correlation'; 
end

g = ~any(isnan([Xin,Ytarget]),2); 
X = Xin(g,:); 
Y = Ytarget(g,:); 

options.verbosity = 0;
options.maxiter = 1000; 
warning('off', 'manopt:getHessian:approx');
options.maxiter = 1000; 

problem = []; 
problem.M = stiefelfactory(size(Y,2),1); 

switch method
    case 'Correlation'
    
        problem.cost = @(Q) ((X*(pinv(X'*X)*X'*(Y*Q)))'*Y*Q)./sqrt(((X*(pinv(X'*X)*X'*(Y*Q)))'*(X*(pinv(X'*X)*X'*(Y*Q)))).*((Y*Q)'*(Y*Q))); 
        problem.egrad = @(Q) 2/sqrt(((Q'*(Y'*X*(pinv(X'*X))*X'*(X*(pinv(X'*X))*X'*(Y*Q))))*(Q'*(Y'*(Y*Q)))))*(Y'*(X*(pinv(X'*X))*X'*(Y*Q))) - ...
            (2*(((0.5./sqrt(((Q'*(Y'*X*(pinv(X'*X))*X'*(X*(pinv(X'*X))*X'*(Y*Q))))*(Q'*(Y'*(Y*Q))))))/((Q'*(Y'*X*(pinv(X'*X))*X'*...
            (X*(pinv(X'*X))*X'*(Y*Q))))*(Q'*(Y'*(Y*Q)))))*(Q'*(Y'*(X*(pinv(X'*X))*X'*(Y*Q))))*(Q'*(Y'*(Y*Q)))*(Y'*X*(pinv(X'*X))*X'*...
            (X*(pinv(X'*X))*X'*(Y*Q))))+2*(((0.5./sqrt(((Q'*(Y'*X*(pinv(X'*X))*X'*(X*(pinv(X'*X))*X'*(Y*Q))))*(Q'*(Y'*(Y*Q))))))/...
            ((Q'*(Y'*X*(pinv(X'*X))*X'*(X*(pinv(X'*X))*X'*(Y*Q))))*(Q'*(Y'*(Y*Q)))))*(Q'*(Y'*(X*(pinv(X'*X))*X'*(Y*Q))))*(Q'*(Y'*X*(pinv(X'*X))*...
            X'*(X*(pinv(X'*X))*X'*(Y*Q))))*(Y'*(Y*Q))));

        costfun = @(Y,X,Q) ((X*(pinv(X'*X)*X'*(Y*Q)))'*Y*Q)./sqrt(((X*(pinv(X'*X)*X'*(Y*Q)))'*(X*(pinv(X'*X)*X'*(Y*Q)))).*((Y*Q)'*(Y*Q)));
    case 'Residuals'
        
        problem.cost = @(Q) -sum(((X*(pinv(X'*X)*X'*(Y*Q)))-Y*Q).^2); 
        problem.egrad = @(Q) -(2*Y'*X*pinv(X'*X)*X'*(X*pinv(X'*X)*X'*(Y*Q) - (Y*Q)) - 2*Y'*(X*pinv(X'*X)*X'*(Y*Q) - (Y*Q))); 
        
        costfun = @(Y,X,Q) -sum(((X*(pinv(X'*X)*X'*(Y*Q)))-Y*Q).^2); 
 
    case 'VAF'
        
        Y = Y-nanmean(Y); 
        X = X-nanmean(X); 
        
        problem.cost = @(Q) -sum(((X*(pinv(X'*X)*X'*(Y*Q)))-Y*Q).^2)./(sum((Y*Q).^2)); 
        problem.egrad = @(Q) -((2./sum((Y*Q).^2))*Y'*X*pinv(X'*X)*X'*(X*pinv(X'*X)*X'*(Y*Q)-(Y*Q))-(2./sum((Y*Q).^2))*Y'*(X*pinv(X'*X)*X'*(Y*Q)-(Y*Q))-1./(sum((Y*Q).^2).^2)*(((Y*Q)'*X*pinv(X'*X)*X' + (-(Y*Q))').^2)*(2*ones(size(Y,1),1))*Y'*(Y*Q)); 

        costfun = @(Y,X,Q) 1-sum(((X*(pinv(X'*X)*X'*(Y*Q)))-Y*Q).^2)./(sum((Y*Q).^2)); 
        
    case 'VAF_no_offset'
               
        problem.cost = @(Q) -sum(((X*(pinv(X'*X)*X'*(Y*Q)))-Y*Q).^2)./(sum((Y*Q).^2)); 
        problem.egrad = @(Q) -((2./sum((Y*Q).^2))*Y'*X*pinv(X'*X)*X'*(X*pinv(X'*X)*X'*(Y*Q)-(Y*Q))-(2./sum((Y*Q).^2))*Y'*(X*pinv(X'*X)*X'*(Y*Q)-(Y*Q))-1./(sum((Y*Q).^2).^2)*(((Y*Q)'*X*pinv(X'*X)*X' + (-(Y*Q))').^2)*(2*ones(size(Y,1),1))*Y'*(Y*Q)); 

        costfun = @(Y,X,Q) 1-sum(((X*(pinv(X'*X)*X'*(Y*Q)))-Y*Q).^2)./(sum((Y*Q).^2)); 
end

fprintf('finding dimension of Y that is most poorly fit by X...'); 
Qf = trustregions(problem,[],options); 

best_b = X\(Y*Qf); 
r = corr(X*best_b,Y*Qf); 
fprintf('done (r = %.2f)\n',r); 

end