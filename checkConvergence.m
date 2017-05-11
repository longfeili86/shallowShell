function [isConverged,res,p]=checkConvergence(step,tol,x,n,tStep)
% check the convergence of iterations
%  stopping criteria: ||xnew-xcur||_{\infty} < tol
%

infoPrefix = '--checkConvergence--: '; % all info displayed by this function includes this prefix

[prev2,prev,cur,new] = step2IterLevels(step); 

Nphi=1:n; % phi solutions
Nw=n+1:2*n;  % w solutions
Nlambda=2*n+1:2*n+3; % lambda solutions

res=max(abs(x([Nphi,Nw],new)-x([Nphi,Nw],cur)));
isConverged=false;
if(res<tol)
   isConverged=true; 
end

%if(res>1e15)
%    error('Iteration is diverging. I am going to abort.');
%end

p=0; % estimated iteration conv rate
if(step>4)
    resCur=max(abs(x([Nphi,Nw],cur)-x([Nphi,Nw],prev)));
    resPrev=max(abs(x([Nphi,Nw],prev)-x([Nphi,Nw],prev2)));

    p = log(res/resCur)/log(resCur/resPrev);
end

%print info:
fprintf('%sStep=%i, res=%e, time used=%e sec, p=%f\n',infoPrefix,step,res,tStep,p);


end