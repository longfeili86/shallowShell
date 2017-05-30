function [x,fval,exitflag,output]=picardSolve(n,x0,W0,Index,mtx,parameters,myGrid,RHSphi,RHSw,R,xOld,dx)
% solve the coupled problem using picard iteration methods
% Input: Wi is the initial guess
% -- Longfei Li

infoPrefix = '--picardSolve--: '; % all info displayed by this function includes this prefix


bcType=parameters.bcType;
maxIter=parameters.maxIter;
tol=parameters.tol;
useLU=parameters.useLU;
solver=parameters.solver;
usePAC=parameters.usePAC; 
% new: we add implicitFactor to combine exPicard and imPicard
implicitFactor=parameters.implicitFactor;

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


Aphi = getMTX_phiEqn(Index,mtx,parameters); % Aphi does not change
if(useLU)
    fprintf('%suse LU decomposition for phi equation\n',infoPrefix);
    [Lphi,Uphi]=lu(Aphi,0.);
end
% Aw does not change for exPicard, so we compute it only once before the
% iteration starts
if(strcmp(solver,'exPicard'))
    Aw=getMTX_wEqn(Index,mtx,parameters,0.*W0);  % we don't need the PHI argument for exPicard, so pass zero 
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
xSol(:,new)=x0;
while(~isConverged && step<=maxIter)
    % time each step:
    tStart=tic;
    step=step+1;
    [prev2,prev,cur,new] = step2IterLevels(step); %shift index
    
    % solve phi equation 
    Rphi=RHSphi(xSol(Nw,cur),xSol(Nphi,cur));
    if(usePAC)
        vxi(1:length(Nphi),1)=xSol(end,cur); % we need add to cur vxi to each of phi eqautions
        bcTypeSaved=parameters.bcType; % save bcType
        if(parameters.bcType==3 || parameters.bcType==5)
            % for free bc and CF bc, the phi eqn uses the clamped bc: phi=dphidn=0
            % so we overwrite the bcType value here to reuse
            % the assignBoundaryConditionsRHS fucntion for the phi eqn as well
            parameters.bcType=2; 
        end
        vxi=assignBoundaryConditionsRHS(vxi,Index,parameters);
        parameters.bcType=bcTypeSaved; % reset bcType to saved value
        Rphi=Rphi-vxi; % add current xi to Rphi 
    end
    if(useLU)
        y=Lphi\RHSphi(xSol(Nw,cur),xSol(Nphi,cur));
        xSol(Nphi,new)=Uphi\y;
    else
        xSol(Nphi,new)=Aphi\Rphi;
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
    
    if(usePAC) % we need to solve the normalizeiton equation
        ds=parameters.ds;
        xSol(2*n+nadd,new)=(ds-dx(1:end-1)'*(xSol(1:end-1,new)-xOld(1:end-1)))/dx(end)+xOld(end);
        %fprintf('--debug-- xSol(end,new)=%e and xSol(end,cur)=%e\n',xSol(end,new),xSol(end,cur));
    end
    
    
    tStep(step)=toc(tStart);
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
%get fval so we can compare with newton or fsolve:
fval=FEvaluation(x,n,Aphi,Aw,RHSphi,RHSw,R,Index,parameters,xOld,dx);


end