function problem=setupFSolveProblem(myGrid,Index,mtx,RHSphi,RHSw,R,W0,x0,n,parameters)
% this function setup the problem for fsolve
% input: 
%   x0 is the initial guess for both phi and w
%
% -- Longfei Li

bcType=parameters.bcType;
tol=parameters.tol;
maxIter=parameters.maxIter;
% setup matrices
% bc for MTXs are  implemented in side of getMTX functions
Aphi = getMTX_phiEqn(Index,mtx,parameters);
Aw = getMTX_wEqn(Index,mtx,parameters,0.*W0);  % we don't need the PHI argument for fsolve, so pass zero 
if(bcType==3)
    % A is the augmented matrix and Q is the kernal
    Aw=removeMatrixSingularity(Aw,myGrid,Index);
end    


% setup problem for fsolve
problem.options = optimoptions('fsolve','Display','iter-detailed',...
    'Algorithm', 'trust-region-dogleg','TolX',tol,'TolFun',tol,...
    'Jacobian','on','MaxIter',maxIter); % we know how to compute Jacobiam now. So turn it on
problem.objective = @(x) fsolveFun(x,n,W0,Aphi,Aw,mtx,Index,RHSphi,RHSw,R,parameters);
problem.x0 = x0;
problem.solver = 'fsolve';

end