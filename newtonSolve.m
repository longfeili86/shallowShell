function x=newtonSolve(n,PHI0,W0,Index,mtx,parameters,myGrid,RHSphi,RHSw,R)
% solve the coupled problem using newton iteration methods
%
% -- Longfei Li

infoPrefix = '--newtonSolve--: '; % all info displayed by this function includes this prefix


bcType=parameters.bcType;
maxIter=parameters.maxIter;
tol=parameters.tol;
useLU=parameters.useLU;
solver=parameters.solver;

isConverged=false;
numberOfLevels=4; % we need 4 levels to estimate conv rate of iteration

Nphi=1:n; % phi solutions
Nw=n+1:2*n;  % w solutions
Nlambda=2*n+1:2*n+3; % lambda solutions
nadd=0; % number of additional variables
if(bcType==3)
    nadd=3;
end
xSol=zeros(2*n+nadd,numberOfLevels); % holder for solutions
F=zeros(2*n+nadd,numberOfLevels); % holder for F vectors

Aphi = getMTX_phiEqn(Index,mtx,parameters); % Aphi does not change
Aw=getMTX_wEqn(Index,mtx,parameters,PHI0);   
if(bcType==3)
    quiet=false;
    Aw=removeMatrixSingularity(Aw,myGrid,Index,quiet); 
end

step=0;
% do this to avoid copying data for new stage
% solution at step 0 is stored in new
[prev2,prev,cur,new] = step2IterLevels(step); % we need four stages to estimate convegence rate
xSol(:,new)=getInitialGuess(n,W0,PHI0,Aphi,RHSphi,parameters);
F(:,new)=FEvaluation(xSol(:,cur),n,Aphi,Aw,RHSphi,RHSw,R,parameters);
while(~isConverged && step<=maxIter)
    tStart=tic;
    step=step+1;
    [prev2,prev,cur,new] = step2IterLevels(step); %shift index
    
    % evaluate Jacobian at current time
    Jcur=JEvaluation(xSol(:,cur),n,W0,Aphi,Aw,mtx,Index,parameters);   
    Res = -Jcur\F(:,cur);
    xSol(:,new)=Res+xSol(:,cur);
    F(:,new)=FEvaluation(xSol(:,new),n,Aphi,Aw,RHSphi,RHSw,R,parameters);

    tStep(step)=toc(tStart);
    
%     res=sqrt(sum(Res.^2));
%     fval=sqrt(sum(Fcur.^2)); % use L2 norm for now
%     fprintf('%sres=%e, fval=%e, tol=%e\n',infoPrefix,res,fval,tol);
%     if(res<tol || fval<tol)
%         isConverged=true;
%     end

    [isConverged,res,p(step)]=checkConvergence(step,tol,xSol,n,tStep(step));
end

if(~isConverged)
    fprintf('%sIteration does not converge after %i steps (maxIter=%i)\n',infoPrefix,step,maxIter);
    fprintf('%sres=%e, tol=%e\n',infoPrefix,res,tol);
    error('abort newtonSolve');
else
    fprintf('%s-----------Iteration Summary-----------\n',infoPrefix); 
    fprintf('%sIteration converges after %i steps (maxIter=%i)\n',infoPrefix,step,maxIter); 
    fprintf('%sAverage time used per iteration step is %e sec\n',infoPrefix,mean(tStep)); 
    fprintf('%sEstimated computational order of convergence p=%f\n',infoPrefix,mean(p(5:end)));    
    x=xSol(:,new); % solution
    %get fval so we can compare with newton or fsolve:
    fval=sqrt(sum(F(:,new).^2)); % use L2 norm for now
    fprintf('%sres=%e, fval=%e, tol=%e\n',infoPrefix,res,fval,tol);
    fprintf('%s---------------------------------------\n',infoPrefix); 
end



end