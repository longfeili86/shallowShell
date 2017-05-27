function [x,fval,exitflag,output]=newtonSolve(n,x0,W0,Index,mtx,parameters,myGrid,RHSphi,RHSw,R,xOld,dx)
% solve the coupled problem using newton iteration methods
% Input: x0 is the initial guess for both phi and w 
% -- Longfei Li

infoPrefix = '--newtonSolve--: '; % all info displayed by this function includes this prefix


bcType=parameters.bcType;
maxIter=parameters.maxIter;
tol=parameters.tol;
useLU=parameters.useLU;
solver=parameters.solver;
usePAC=parameters.usePAC;
relaxFactor=parameters.relaxFactor;

isConverged=false;
numberOfLevels=4; % we need 4 levels to estimate conv rate of iteration

Nphi=1:n; % phi solutions
Nw=n+1:2*n;  % w solutions
Nlambda=2*n+1:2*n+3; % lambda solutions
nadd=0; % number of additional variables
if(bcType==3)
    nadd=3;
end
if(usePAC)
    nadd=nadd+1;
end
xSol=zeros(2*n+nadd,numberOfLevels); % holder for solutions
F=zeros(2*n+nadd,numberOfLevels); % holder for F vectors

Aphi = getMTX_phiEqn(Index,mtx,parameters); % Aphi does not change
Aw=getMTX_wEqn(Index,mtx,parameters,0.*W0); % we don't need the PHI argument for newton, so pass zero   
if(bcType==3)
    quiet=false;
    Aw=removeMatrixSingularity(Aw,myGrid,Index,quiet); 
end

step=0;
% do this to avoid copying data for new stage
% solution at step 0 is stored in new
[prev2,prev,cur,new] = step2IterLevels(step); % we need four stages to estimate convegence rate
xSol(:,new)=x0;
F(:,new)=FEvaluation(xSol(:,new),n,Aphi,Aw,RHSphi,RHSw,R,Index,parameters,xOld,dx);
while(~isConverged && step<=maxIter)
    tStart=tic;
    step=step+1;
    [prev2,prev,cur,new] = step2IterLevels(step); %shift index
    
    % evaluate Jacobian at current time
    Jcur=JEvaluation(xSol(:,cur),n,W0,Aphi,Aw,mtx,Index,parameters,dx);   
    Res = -Jcur\F(:,cur);
    xSol(:,new)=relaxFactor*Res+xSol(:,cur); % relax netwon
    F(:,new)=FEvaluation(xSol(:,new),n,Aphi,Aw,RHSphi,RHSw,R,Index,parameters,xOld,dx);

    tStep(step)=toc(tStart);
    
%     res=sqrt(sum(Res.^2));
%     fval=sqrt(sum(Fcur.^2)); % use L2 norm for now
%     fprintf('%sres=%e, fval=%e, tol=%e\n',infoPrefix,res,fval,tol);
%     if(res<tol || fval<tol)
%         isConverged=true;
%     end

    [isConverged,res,p(step)]=checkConvergence(step,tol,xSol,n,tStep(step));
    if(isConverged==-9999)
        break;
    end       
end

%prepare output stuff
output.algorithm=solver;
output.iterations=step;

message=sprintf('res=%e,tol=%e, maxIter=%i\n',res,tol,maxIter);
message=sprintf('%sAverage time used per iteration step is %e sec\n',message,mean(tStep));
message=sprintf('%sEstimated computational order of convergence p=%f\n',message,mean(p(5:end)));
output.message=message;

if(isConverged==-9999)
    exitflag=-9999;
elseif(isConverged>0)
    exitflag=1;
else
    exitflag=0;
end 


x=xSol(:,new); % solution (could be unconverged)
fval=F(:,new); 



end