function x=picardSolve(n,PHI0,W0,Index,mtx,parameters,myGrid,RHSphi,RHSw,R)
% solve the coupled problem using picard iteration methods
%
% -- Longfei Li

infoPrefix = '--picardSolve--: '; % all info displayed by this function includes this prefix


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


Aphi = getMTX_phiEqn(Index,mtx,parameters); % Aphi does not change
if(useLU)
    fprintf('%suse LU decomposition for phi equation\n',infoPrefix);
    [Lphi,Uphi]=lu(Aphi,0.);
end
% Aw does not change for exPicard, so we compute it only once before the
% iteration starts
if(strcmp(solver,'exPicard'))
    Aw=getMTX_wEqn(Index,mtx,parameters,PHI0);   
    if(bcType==3)
        quiet=false;
        Aw=removeMatrixSingularity(Aw,myGrid,Index,quiet); 
    end
    if(useLU)
        fprintf('%suse LU decomposition for w equation\n',infoPrefix);    
        [Lw,Uw]=lu(Aw,0.);
    end
end



step=0;
% do this to avoid copying data for new stage
% solution at step 0 is stored in new
[prev2,prev,cur,new] = step2IterLevels(step); % we need four stages to estimate convegence rate
xSol(:,new)=getInitialGuess(n,W0,PHI0,Aphi,RHSphi,parameters);
while(~isConverged && step<=maxIter)
    % time each step:
    tStart=tic;
    step=step+1;
    [prev2,prev,cur,new] = step2IterLevels(step); %shift index
    
    % solve phi equation 
    if(useLU)
        y=Lphi\RHSphi(xSol(Nw,cur),xSol(Nphi,cur));
        xSol(Nphi,new)=Uphi\y;
    else
        xSol(Nphi,new)=Aphi\RHSphi(xSol(Nw,cur),xSol(Nphi,cur));
    end
    
    % solve w equation
    if(strcmp(solver,'imPicard')) % recompute Aw at evert step for imPicard 
        quiet=step>1; % supress repeating message from each step
        Aw=getMTX_wEqn(Index,mtx,parameters,xSol(Nphi,new),quiet);   
        if(bcType==3)
            % A is the augmented matrix and Q is the kernal
            Aw=removeMatrixSingularity(Aw,myGrid,Index,quiet); 
        end
    end
    Rw=RHSw(xSol(Nw,cur),xSol(Nphi,new)); 
    if(bcType==3)
        Rw=[Rw;R]; % additional rhs of w equation
    end
    if(useLU && strcmp(solver,'exPicard')) % we can use LU here for exPicard
        y=Lw\Rw;
        xTemp=Uw\y;
    else
        xTemp=Aw\Rw;
    end
    xSol(Nw,new)=xTemp(1:n);
    if(bcType==3)
       xSol(Nlambda,new) = xTemp(n+1:n+3); 
    end
    tStep(step)=toc(tStart);
    [isConverged,res,p(step)]=checkConvergence(step,tol,xSol,n,tStep(step));
 
end

if(~isConverged)
    fprintf('%sIteration does not converge after %i steps (maxIter=%i)\n',infoPrefix,step,maxIter);
    fprintf('%sres=%e, tol=%e\n',infoPrefix,res,tol);
    error('abort picardSolve');
else
    fprintf('%s-----------Iteration Summary-----------\n',infoPrefix); 
    fprintf('%sIteration converges after %i steps (maxIter=%i)\n',infoPrefix,step,maxIter); 
    fprintf('%sAverage time used per iteration step is %e sec\n',infoPrefix,mean(tStep)); 
    fprintf('%sEstimated computational order of convergence p=%f\n',infoPrefix,mean(p(5:end)));
    x=xSol(:,new); % solution
    %get fval so we can compare with newton or fsolve:
    F=FEvaluation(x,n,Aphi,Aw,RHSphi,RHSw,R,parameters);
    fval=sqrt(sum(F.^2)); % use L2 norm for now
    fprintf('%sres=%e, fval=%e, tol=%e\n',infoPrefix,res,fval,tol);
    fprintf('%s---------------------------------------\n',infoPrefix); 
end


end