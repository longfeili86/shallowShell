function problem=setupFSolveProblem(myGrid,Index,mtx,RHSphi,RHSw,R,W0,PHI0,n,parameters)
% this function setup the problem for fsolve
%
% -- Longfei Li

bcType=parameters.bcType;
tol=parameters.tol;

% setup matrices
% bc for MTXs are  implemented in side of getMTX functions
Aphi = getMTX_phiEqn(Index,mtx,parameters);
Aw = getMTX_wEqn(Index,mtx,parameters,PHI0);
if(bcType==3)
    % A is the augmented matrix and Q is the kernal
    Aw=removeMatrixSingularity(Aw,myGrid,Index);
end    


% setup problem for fsolve
problem.options = optimoptions('fsolve','Display','iter-detailed',...
    'Algorithm', 'trust-region-dogleg','TolX',tol,'TolFun',tol,...
    'Jacobian','on'); % we know how to compute Jacobiam now. So turn it on
problem.objective = @(x) fsolveFun(x,n,W0,Aphi,Aw,mtx,Index,RHSphi,RHSw,R,parameters);
problem.x0 = getInitialGuess(n,W0,PHI0,Aphi,RHSphi,parameters);
problem.solver = 'fsolve';

end