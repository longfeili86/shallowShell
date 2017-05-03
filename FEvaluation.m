function F=FEvaluation(x,n,Aphi,Aw,RHSphi,RHSw,R,parameters)
% evaluate the F function
%
% -- Longfei Li


Nphi=1:n; % phi solutions
Nw=n+1:2*n;  % w solutions
Nlambda=2*n+1:2*n+3; % lambda solutions



phi=x(Nphi);
w=x(Nw);

F(Nphi,1)= Aphi*phi-RHSphi(w,phi); % phi eqn
if(parameters.bcType==3) % w eqn is singular, we need extra equations for lambda
    lambda=x(Nlambda); %  extra unknowns
    F([Nw,Nlambda],1) = Aw*[w;lambda]-[RHSw(w,phi);R]; % w + lambda eqn
else
    F(Nw,1) = Aw*w-RHSw(w,phi); % w eqn
end


end